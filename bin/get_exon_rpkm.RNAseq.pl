#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take a file resulting form running the overlap script and
# extracting the total overlap reads. It will calculate the RPKM for each of
# the exons in the file.

use Bio::Range;
use RNAseq_pipeline3 qw(get_fh parse_gff_line get_exon_coverage_1000nt);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list');

# Declare some variables
my $prefix;
my $exondir;
my $exontable;
my $mapper;
my $readlength;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$exontable=$options{'EXONSCLASSTABLE'};
$mapper=$options{'MAPPER'};
$readlength=$options{'READLENGTH'};

# Connect to the database
my $dbh=get_dbh();
my $commondbh=get_dbh(1);

# Get the exon list;
my %exons=%{get_exon_list($commondbh,
			  $exontable)};

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# Get unique maps for each lane
my $mappingtable=$prefix.'_genome_mapping';
my %unique_maps=%{get_unique_maps($dbh,
				  $mappingtable,
				  \%files,
				  $mapper)};

# Read and process the overlap files
foreach my $lane (keys %lanes) {
    my $type;
    if (keys %{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (keys %{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my %rpkms;

    # Get the exon overlap information
    my $exonoverlap=$exondir.'/'.$lane.'.'.$type.'.unique.gtf.overlap.total';

    if (-r $exonoverlap) {
	print STDERR "Processing $exonoverlap\n";
    } else {
	die "Can't read $exonoverlap\n";
    }

    # Read in the exon coverage
    my %exons_coverage;
    get_exon_coverage_1000nt($exonoverlap,
			    \%exons_coverage,
			     $readlength);

    # For each of the exons in the exon list print the RPKM by normalizing the
    # reads per 1000 bases by the number of uniquely mapped reads
    my $outfile=$exondir.'/'.$lane.'.'.$type.'.exon.rpkm.txt.gz';
    process_exons(\%exons,
		  \%exons_coverage,
		  $outfile,
		  $unique_maps{$lane},
		  $lane);
}

exit;

sub get_unique_maps {
    my $dbh=shift;
    my $maptable=shift;
    my $files=shift;
    my $mapper=shift;

    my %unique_maps;
    my %lane2pair;

    my ($query,$sth,$count);

    print STDERR "Getting unique mappings from $maptable...";
    foreach my $file (keys %{$files}) {
	$lane2pair{$files->{$file}->[1]}=$files->{$file}->[0];
    }

    $query ='SELECT LaneName, uniqueReads ';
    $query.="FROM $maptable";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count == keys %{$files}) {
	die "Wrong number of mapped readsd\n";
    }

    while (my ($lane,$unique)=$sth->fetchrow_array()) {
	my $pair=$lane2pair{$lane};
	$unique_maps{$pair}+=$unique;
    }
    print STDERR "done\n";

    foreach my $lane (keys %unique_maps) {
	print STDERR join("\t",
			  $lane,
			  $unique_maps{$lane}),"\n";
    }

    return(\%unique_maps);
}

sub get_exon_list {
    my $dbh=shift;
    my $table=shift;

    my %exons;

    print STDERR "Retrieving exon list from $table...\n";
    my ($query,$sth,$count);
    $query ='SELECT DISTINCT exon_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count) {
	warn "No records retrieved from $table\n";
    }

    while (my ($exon_id)=$sth->fetchrow_array()) {
	$exons{$exon_id}=1;
    }
    print STDERR $count,"\tExons retrieved\n";

    return(\%exons);
}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}


sub process_exons {
    my $exonlist=shift;
    my $exoncoverage=shift;
    my $output=shift;
    my $uniquemaps=shift;
    my $lane=shift;

    print STDERR "Collecting events...\n";

    my $outfh=get_fh($output,1);

    foreach my $exon_id (keys %{$exonlist}) {
	my $rpkm=0;

	if ($exoncoverage->{$exon_id}) {
	    $rpkm+=$exoncoverage->{$exon_id};
	}

	# Normalize by read number
	$rpkm=sprintf "%.3f",($rpkm * 1000000) / $uniquemaps;

	# Print results only if they are positive
	if ($rpkm != 0) { 
	    print $outfh join("\t",
			      $exon_id,
			      $rpkm,
			      $lane),"\n";
	}
    }
    close($outfh);

    print STDERR "\rDone\n";
}
