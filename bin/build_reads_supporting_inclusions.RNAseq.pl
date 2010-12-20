#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script should take information from the exon_inclusion table, and build
# a table that contains the same information, but with the exon locations
# separated into chr start and end (also in the place of the lane field
# we will have the sample field) Also in the place of the rpkms we will just
# have the read counts

# Load some modules
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file get_dbh read_file_list get_coords_from_exon_sub);

# Declare som variables
# Define some variables
my $project;
my $prefix;
my $exclusionfn;
my $genomedir;
my $junctiondir;
my $exondir;
my $mapper;
my $juncfile;
my $readlength;

# Read the configuration file
my %options=%{read_config_file()};
$exclusionfn=$options{'EXCLUSIONFILE'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$genomedir=$options{'GENOMEDIR'};
$junctiondir=$options{'JUNCTIONSDIR'};
$exondir=$options{'EXONDIR'};
$mapper=$options{'MAPPER'};
$juncfile=$options{'JUNCTIONSFASTA'};
$readlength=$options{'READLENGTH'};

# Check the input
unless ($exclusionfn) {
    die "No exclusion file???\n";
}

# Connect to the database
my $dbh=get_dbh();
my $dbhcommon=get_dbh(1);

# Get some subroutines
my $table=get_junctions_table($dbh,
			      $prefix);

my $mappingtable=$prefix.'_genome_mapping';
*get_inclusion_exons=get_inclusion_junctions_sub($dbhcommon,
						 $table);
*get_exon_coords=get_coords_from_exon_sub();

my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};
my %groups=%{get_groups(\%files)};

# Read and process the files
foreach my $group (keys %groups) {
    my %exons_coverage;
    my %juncs_coverage;
    foreach my $lane (@{$groups{$group}}) {
	my $type;
	if (keys %{$lanes{$lane}} == 1) {
	    $type='single';
	} elsif (keys %{$lanes{$lane}} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type\n";
	}
	my %genes;
	
	# Get the exon information
	my $exonfilename=$exondir.'/'.$lane.'.'.$type.'.unique.gtf.overlap.total';
	
	if (-r $exonfilename) {
	    print STDERR "Processing $exonfilename\n";
	} else {
	    die "Can't read $exonfilename\n";
	}
	
	# Read in the exon coverage
	get_exon_coverage($exonfilename,
			  \%exons_coverage);
	
	# Get the junction information
	my @junctionfns;
	
	foreach my $track (keys %{$lanes{$lane}}) {
	    my $juncfilename=$junctiondir.'/'.$track.'.single.unique.overlap.total';
	    push @junctionfns, $juncfilename;
	}
	
	# Process each of the junction files
	foreach my $juncfilename (@junctionfns) {
	    if (-r $juncfilename) {
		print STDERR "Processing $juncfilename\n";
	    } else {
		die "Can't read $juncfilename\n";
	    }
	    get_junctions_coverage($juncfilename,
				   \%juncs_coverage);
	}
    }

    # Go through the exclusion file and for each of the exons print the
    # Gene, exon, number of inclusion reads and number of exclusion reads
    my $outfile=$exondir.'/'.$group.'.inc.exc.reads.txt';
    process_exons($exclusionfn,
		  \%juncs_coverage,
		  \%exons_coverage,
		  $outfile,
		  $group);
}


exit;

sub process_exons {
    my $filename=shift;
    my $junctions=shift;
    my $exons=shift;
    my $output=shift;
    my $lane=shift;

    print STDERR "Collecting $lane events...\n";

    my $infh=get_fh($filename);
    my $outfh=get_fh($output,1);

    while (my $line=<$infh>) {
	my ($inc_rate_ex,$inc_rate_junc,$ex_rate)=(0,0,0);
	chomp($line);
	my ($gene,$exon_id,$exclusion)=split("\t",$line);
	my @inclusion_exons=get_inclusion_exons($exon_id);
	my @exclusion_exons=split(',',$exclusion);

	if ($exons->{$exon_id}) {
	    $inc_rate_ex+=$exons->{$exon_id};
	}

	foreach my $inc (@inclusion_exons) {
	    if ($junctions->{$inc}) {
		$inc_rate_junc+=$junctions->{$inc};
	    }
	}

	foreach my $exc (@exclusion_exons) {
	    if ($junctions->{$exc}) {
		$ex_rate+=$junctions->{$exc};
	    }
	}

	my $inclusion='\N';
	my $all_inclusion=$inc_rate_junc + $inc_rate_ex;
	my $all_events=$inc_rate_junc + $inc_rate_ex + $ex_rate;
	if ($all_events > 0) {
	    $inclusion= sprintf "%.3f",$all_inclusion / $all_events;
	}

	# Print results
	my @coords=@{get_exon_coords($exon_id)};
	print $outfh join("\t",
			  $exon_id,
			  @coords[0,1,2],
			  $inc_rate_ex,
			  $inc_rate_junc,
			  $ex_rate,
			  $inclusion,
			  $lane),"\n";
    }
    close($outfh);
    close($infh);

    print STDERR "\rDone\n";

}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}

sub get_junctions_coverage {
    my $infn=shift;
    my $features=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage)=split("\t",$line);
	$features->{$feature}+=$coverage;
    }

    print STDERR "done\n";
}

sub get_exon_coverage {
    my $infn=shift;
    my $features=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage,$length,$type)=split("\t",$line);

	if ($coverage < 0) {
	    die "Coverage is below zero ????\n";
	}
	$features->{$feature}+=$coverage;
    }

    print STDERR "done\n";
}

sub get_junctions_table {
    my $dbh=shift;
    my $prefix=shift;

    my $table=$prefix.'_junctions';

    my ($query,$sth,$count);
    $query ='SELECT table_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count == 1) {
	warn "Incorrect number of records retrieved from $table\n";
    }

    my ($junctable)=$sth->fetchrow_array();

    return($junctable);
}

sub get_inclusion_junctions_sub {
    my $dbh=shift;
    my $table=shift;

    my ($query,$sth);
    $query ='SELECT DISTINCT chr, start, end ';
    $query.="FROM $table ";
    $query.='WHERE exon1_id = ? OR exon2_id = ?';
    $sth=$dbh->prepare($query);

    my $no_hits=0;
    my $processed=0;
    my $get_exons=sub {
	my $exon=shift;
	my @junctions;
	my $count=$sth->execute($exon,$exon);

	while (my ($chr,$start,$end)=$sth->fetchrow_array()) {
	    my $junction=join('_',$chr,$start,'splice',$end);
	    push @junctions, $junction;
	}
	return(@junctions);
    };
    return($get_exons);
}

sub get_groups {
    my $files=shift;
    my %groups;
    
    foreach my $file (keys %{$files}) {
	my $group=$files->{$file}->[2] || 'All';
	push @{$groups{$group}},$files->{$file}->[0];
    }

    return(\%groups);
}
