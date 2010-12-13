#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take a file resulting form running the overlap script and
# extracting the total overlap reads and a file resulting from projecting all
# the genes onto the genome. It will calculate the RPKM for each of
# the genes in the file

use Bio::Range;
use RNAseq_pipeline2 qw(get_fh parse_gff_line get_feature_overlap);
use RNAseq_pipeline_settings ('read_config_file','get_dbh','read_file_list');

# Declare some variables
my $prefix;
my $exondir;
my $genomedir;
my $mapper;
my $projfile;
my $tmpdir;
my $stranded;
my $readlength;
my $threshold=10;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$genomedir=$options{'GENOMEDIR'};
$mapper=$options{'MAPPER'};
$projfile=$genomedir.'/'.$prefix.'.proj.gtf';
$tmpdir=$options{'LOCALDIR'};
$stranded=$options{'STRANDED'};
$readlength=$options{'READLENGTH'};

# Connect to the database
my $dbh=get_dbh();

# Get the exon list;
my %genes;
my %exon2gene;
my %exons=%{get_proj_list_from_gtf($projfile,
				   \%genes)};

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

    # Get the proj overlap information
    my $projoverlap=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.proj.overlap.total';

    if (-r $projoverlap) {
	print STDERR "Processing $projoverlap\n";
    } else {
	die "Can't read $projoverlap\n";
    }

    # Read in the proj coverage
    my %gene_coverage;
    get_proj_coverage_1000nt($projoverlap,
			     \%exons,
			     \%gene_coverage,
			     $readlength,
			     $threshold);

    # For each of the exons in the exon list print the RPKM by normalizing the
    # reads per 1000 bases by the number of uniquely mapped reads
    my $outfile=$genomedir.'/'.$lane.'.'.$type.'.proj.rpkm.txt.gz';
    process_exons(\%genes,
		  \%gene_coverage,
		  $outfile,
		  $unique_maps{$lane},
		  $lane);
}

exit;

# This will get all the exons from the annotation file
sub get_proj_list_from_gtf {
    my $projfile=shift;
    my $genes=shift;
    my $repeated_exons='rep.exons.txt';

    my %exons;
    my %remove;

    print STDERR "Extracting projection list lists from $projfile...";
    my $repeatfh=get_fh($repeated_exons,1);
    my $exonfh=get_fh($projfile);
    while (my $line=<$exonfh>) {
	my %line=%{parse_gff_line($line)};

	# Complain if they give us something that is not an exon
	unless ($line{'type'} eq 'proj') {
	    warn "Non projection line: $line";
	}

	# get the strand
	my $strand=$line{'strand'};
	if ($strand eq '+') {
	    $strand=1;
	} elsif ($strand eq '-') {
	    $strand= -1;
	} else {
	    # complain
	    warn "Unstranded exon: Strand $line{'strand'},$line\n";
	}

	# Get the exon_id
	my $exon_id=join('_',
			 $line{'chr'},
			 $line{'start'},
			 $line{'end'},
			 $strand);

	# Get the gene_id
	my $gene_id=$line{'feature'}{'gene_id'};
	$gene_id=~s/"//g;

	# Check we have everything
	unless ($gene_id && $exon_id) {
	    warn "Problem parsing $line\n";
	}

	if ($exons{$exon_id} && 
	    ($exons{$exon_id} ne $gene_id)) {
	    $remove{$exon_id}=1;
	} else {
	    $exons{$exon_id}=$gene_id;
	    $genes->{$gene_id}->{$exon_id}=1;
	}
    }
    close($exonfh);
    close($repeatfh);
    print STDERR "done\n";

    # Remove those exons that map to multiple genes
    foreach my $exon (keys %remove) {
	delete $exons{$exon};
    }

    my $count=keys %remove;
    print STDERR $count,"\tProjections mapping to multiple genes removed\n";

    $count=keys %exons;
    print STDERR $count,"\tProjections obtained\n";
    $count=keys %{$genes};
    print STDERR $count,"\tGenes\n";

    return(\%exons);
}

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

sub get_proj_coverage_1000nt {
    my $infn=shift;
    my $proj_lengths=shift;
    my $features=shift;
    my $readlength=shift;
    my $threshold=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    my $nolength=[0,0];
    while (my $line=<$infh>) {
	chomp($line);
	my ($exon_id,$coverage,$length)=split("\t",$line);

	# Normalize by length as reads per 1000 nt
	if ($coverage >= $threshold) {
	    $coverage=($coverage * 1000) / $length;
	    $features->{$exon_id}+=$coverage;
	} elsif ($coverage) {
	    $features->{$exon_id}=-1;
	} else {
	    $nolength->[0]++;
	    $nolength->[1]+=$coverage;
	}
    }
    print STDERR "done\n";

    print STDERR $nolength->[0],"\tFeatures lacked coverage information\n";
#    print STDERR $nolength->[1],"\tWas their length\n";
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
    my $genelist=shift;
    my $exoncoverage=shift;
    my $output=shift;
    my $uniquemaps=shift;
    my $lane=shift;

    print STDERR "Collecting events...\n";

    my $outfh=get_fh($output,1);

    foreach my $gene (keys %{$genelist}) {
	my $rpkm=0;
	my $exon_id='';
	# Choose the highest rpkm for the gene
	foreach my $exon (keys %{$genelist->{$gene}}) {
#	    if ($exoncoverage->{$exon} &&
#		($exoncoverage->{$exon} > $rpkm)) {
	    if ($exoncoverage->{$exon}) {
		
		if ($exoncoverage->{$exon} > $rpkm) {
		    $rpkm=$exoncoverage->{$exon};
		    $exon_id=$exon;
		}
	    } else {
#		print STDERR join("\t",
#				  $exon,
#				  $gene),"\n";
	    }
	}
	# Normalize by read number
	$rpkm=sprintf "%.3f",($rpkm * 1000000) / $uniquemaps;
	
	# Print results only if they are positive
	if ($rpkm > 0) { 
	    print $outfh join("\t",
			      $gene,
			      $exon_id,
			      $rpkm,
			      $lane),"\n";
	} else {
	    print $outfh join("\t",
			      $gene,
			      '-',
			      '-',
			      $lane),"\n";
	}
    }
    close($outfh);

    print STDERR "\rDone\n";
}
