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
# This script should run overlap on a set of read clusters, coming from the
# genome mapping, the junctions mapping and the split mapping filtered for
# selecting those junctions between clusters supported by a certain threshold
# number of reads and in which the junction is also supported by that number
# of reads.

# Load some modules
use RNAseq_pipeline3 ('get_fh','get_feature_overlap_sub');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list');
use Getopt::Long;

# Declare some variables
my $threshold=1;
my $tmpdir;
my $clusterdir;
my $genomedir;
my $exondir;
my $splitdir;
my $junctionsdir;
my $projectid;
my $prefix;
my $stranded;
my $file_list;
my $parallel='default';
my $paralleltmp;

GetOptions('threshold|t=i' => \$threshold,
    );

# Read the options file
my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$genomedir=$options{'PROJECT'}.'/'.$options{'GENOMEDIR'};
$exondir=$options{'EXONDIR'};
$splitdir=$options{'SPLITMAPDIR'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$projectid=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$stranded=$options{'STRANDED'};
$file_list=$options{'FILELIST'};
$paralleltmp=$options{'PROJECT'}.'/'.$options{'LOCALPARALLEL'};
if ($options{'PARALLEL'}) {
    $parallel='parallel';
}

unless ($paralleltmp) {
    $paralleltmp=$options{'PROJECT'}.'/work';
}

# Get the required sub
*get_feature_overlap=get_feature_overlap_sub($parallel,
					     $paralleltmp);

# First ge the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get for each lane the files corresponding to the genome mapping and the
# junctions mapping
my %lane_files=%{get_lane_files(\%files,
				$genomedir)};

my $annotation=$genomedir.'/'.$prefix.'.gene.gtf';

my %coverage;

foreach my $pair (keys %lane_files) {
    # Find the coverage granted by genome mapping;
    my $outfile=$lane_files{$pair};
    $outfile=~s/.gz$//;
    get_feature_overlap($lane_files{$pair},
			$annotation,
			$stranded,
			$outfile,
			$tmpdir);
}

exit;

sub get_lane_files {
    my $files=shift;
    my $dir=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	$lane_files{$pair}{$lane}=1;
    }

    # Determine if we are looking at paired or single reads
    foreach my $pair (keys %lane_files) {
	my @lanes=keys %{$lane_files{$pair}};
	if (@lanes == 1) {
	    # single files
	    print STDERR "$pair set identified as single reads\n";
	    $lane_files{$pair}=$dir.'/'.$lanes[0].'.single.unique.gtf.gz';
	} elsif (@lanes == 2) {
	    # paired files
	    print STDERR "$pair set identified as paired reads\n";
	    $lane_files{$pair}=$dir.'/'.$pair.'.paired.unique.gtf.gz';
	} else {
	    # This should never happen
	    die "Apparently there are more than two files grouped\n";
	}
    }

    return(\%lane_files);
}
