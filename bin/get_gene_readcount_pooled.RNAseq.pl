#!/soft/bin/perl
# DGK

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
# This script should take a file resulting form running the overlap script and
# extracting the total overlap reads and a file resulting from projecting all
# the genes onto the genome. It will calculate the RPKM for each of
# the genes in the file

use Bio::Range;
use Bio::SeqIO;
use RNAseq_pipeline3 ('get_fh','parse_gff_line','get_feature_overlap',
		      'get_exon_list_from_gtf');
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list',
			      'get_gene_from_short_junc_sub');

# Declare some variables
my $prefix;
my $exondir;
my $genomedir;
my $mapper;
my $exonfile;
my $projfile;
my $tmpdir;
my $stranded;
my $readlength;
my $junctiondir;
my $junctionsfasta;
my $debug=1;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$genomedir=$options{'GENOMEDIR'};
$mapper=$options{'MAPPER'};
$exonfile=$exondir.'/'.$prefix.'.exon.gtf';
$projfile=$genomedir.'/'.$prefix.'.proj.gtf';
$tmpdir=$options{'LOCALDIR'};
$stranded=$options{'STRANDED'};
$readlength=$options{'READLENGTH'};
$junctiondir=$options{'JUNCTIONSDIR'};
$junctionsfasta=$options{'JUNCTIONSFASTA'};

# Connect to the database
my $dbh=get_dbh();
my $dbh_common=get_dbh(1);

# Get some usefull subs
*junc2gene=get_gene_from_short_junc_sub($dbh_common);

# Get the exon list;
my %genes;
my %exons=%{get_exon_list_from_gtf($exonfile,
				   \%genes)};

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};
my %groups=%{get_groups(\%files)};

# Get unique maps for each lane
my $mappingtable=$prefix.'_genome_mapping';
my %unique_maps=%{get_unique_maps($dbh,
				  $mappingtable,
				  \%files,
				  $mapper)};

# Read and process the overlap files
foreach my $group (keys %groups) {
    my %gene_juncs_coverage;
    my %gene_coverage;

    # Check if the file is already present and skip it if so
    my $outfile=$genomedir.'/'.$group.'.gene.readcount.pooled.txt.gz';
    if ($debug && (-r $outfile)) {
	print STDERR $outfile,"\tIs present. Skipping...\n";
	next;
    }

    foreach my $lane (@{$groups{$group}}) {
	my $type;
	if (keys %{$lanes{$lane}} == 1) {
	    $type='single';
	} elsif (keys %{$lanes{$lane}} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type $lane\n";
	}
	
	$unique_maps{$group}+=$unique_maps{$lane};
	
	# Get the projected exon overlap information
	my $exonoverlap=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.proj.overlap.total';
	
	if (-r $exonoverlap) {
	    print STDERR "Processing $exonoverlap\n";
	} else {
	    die "Can't read $exonoverlap\n";
	}
	
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
	    get_feature_coverage_junctions($juncfilename,
					   \%gene_juncs_coverage,
					   $readlength);
	}
	
	
	# Read in the exon coverage
	get_gene_coverage($exonoverlap,
			  \%gene_coverage,
			  $readlength);
	
    }
    # For each of the genes in the exon list print the RPKM by normalizing the
    # reads per 1000 bases by the number of uniquely mapped reads
    process_exons(\%genes,
		  \%gene_coverage,
		  \%gene_juncs_coverage,
		  $outfile,
		  $unique_maps{$group},
		  $group);
}

exit;

sub get_feature_coverage_junctions {
    my $infn=shift;
    my $features=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage)=split("\t",$line);

	my $genes=junc2gene($feature);
	foreach my $gene (@{$genes}) {
	    $features->{$gene}+=$coverage;
	}
    }

    print STDERR "done\n";
}

sub get_projected_length {
    my $genes=shift;
    my $lengths=shift;
    my $projfile=shift;
    my $read_length=shift;

    print STDERR "Calculating gene projections\n";

    # If the projected file exists read the projection lengths from it, and
    # else create it
    if (-r $projfile) {
	print STDERR $projfile,"\tpresent\n";
	my $projfh=get_fh($projfile);
	while (my $line=<$projfh>){
	    my %line=%{parse_gff_line($line)};
	    my $gene=$line{'feature'}{'gene_id'};
	    $gene=~s/"//g;
	    my $length=$line{'end'} - $line{'start'} + 1;
	    $lengths->{$gene}+=$length;
	}
	close($projfh);
    } else {
	print STDERR $projfile,"\tAbsent...building\n";
	my $projfh=get_fh($projfile,1);
	foreach my $gene (keys %{$genes}) {
	    my $chr;
	    my @ranges;
	    foreach my $exon (keys %{$genes->{$gene}}) {
		my @location=split('_',$exon);
		my $exon_chr=shift(@location);
		if ($chr &&
		    $exon_chr ne $chr) {
		    warn "$gene has exons in $chr and $exon_chr\n";
		} else {
		    $chr=$exon_chr;
		}
		my $range=Bio::Range->new(-start => $location[0],
					  -end => $location[1],
					  -strand => $location[2]);
		push @ranges,$range;
	    }
	    my @disc_ranges = Bio::Range->disconnected_ranges(@ranges);
	    foreach my $range (@disc_ranges) {
		my $length=$range->length();
		my $strand;
		if ($range->strand() == 1) {
		    $strand='+';
		} elsif ($range->strand() == -1) {
		    $strand='-';
		} else {
		    warn "Unknown strand for $gene\n";
		}
		my $string='gene_id "'.$gene.'";';
		print $projfh join("\t",
				   $chr,
				   'RNAseqPipe',
				   'proj',
				   $range->start(),
				   $range->end(),
				   '.',
				   $strand,
				   '.',
				   $string),"\n";
		$lengths->{$gene}+=$length;
	    }
	}
	close($projfh);
    }
    print STDERR "done\n";
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

sub get_gene_coverage {
    my $infn=shift;
    my $features=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($gene,$coverage,$length)=split("\t",$line);

	# Normalize reads per 1000 nt (the length will be divided later
	if ($coverage) {
	    $features->{$gene}+=$coverage;
	} 
    }
    print STDERR "done\n";
}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
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

sub process_exons {
    my $genelist=shift;
    my $exoncoverage=shift;
    my $junccoverage=shift;
    my $output=shift;
    my $uniquemaps=shift;
    my $lane=shift;

    print STDERR "Collecting events...\n";

    my $outfh=get_fh($output,1);

    foreach my $gene (keys %{$genelist}) {
	my $readcount=0;
	# add coverage from projected exons
	if ($exoncoverage->{$gene}) {
	    $readcount+=$exoncoverage->{$gene};
	}

	# add coverage from junctions
	if ($junccoverage->{$gene}) {
	    $readcount+=$junccoverage->{$gene};
	}

	# Print results only if they are positive
	if ($readcount > 0) { 
	    print $outfh join("\t",
			      $gene,
			      $readcount,
			      $lane),"\n";
	} elsif ($readcount < 0) {
	    die "negative readcounts should not happen\n";
	}
    }
    close($outfh);

    print STDERR "\rDone\n";
}
