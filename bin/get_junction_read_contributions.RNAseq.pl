#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script should take the reads from the single total files in the junctions
# directory and it will produce a single file with the reads grouped by gene

# Load some modules
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list',
			       'get_gene_from_short_junc_sub');

# Declare some variables
my $prefix;
my $stranded;
my $junctiondir;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$stranded=$options{'STRANDED'};
$junctiondir=$options{'JUNCTIONSDIR'};

# Connect to the database
my $dbh_common=get_dbh(1);

# Get some usefull subs
*junc2gene=get_gene_from_short_junc_sub($dbh_common);

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

my %genes;

# Read and process the overlap files
foreach my $lane (keys %lanes) {

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
				       \%genes);
    }

    my $outfn=$junctiondir.'/'.$lane.'.paired.gene.overlap.total';
    my $outfh=get_fh($outfn,1);
    foreach my $gene (keys %genes) {
	    print $outfh join("\t",
			      $gene,
			      $genes{$gene}),"\n";
    }
    close($outfh);
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

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}
