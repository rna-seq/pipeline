#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take the unique mappings to the genome (initially, and
# eventually those to junctions and splitmaps) and build a set of clusters with
# them separated into two groups, clusters with uniquely mapped reads and those
# with multimapped reads.
# If all reads are in the same strand it will conserve the strand, if not it
# it will set it to '.'

# Load some modules
use Getopt::Long;
use RNAseq_pipeline3 ('get_fh','cluster_gff','get_gff_from_junc_id',
		      'get_sorted_gff_fh');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list');

# Declare some variables
my $tmpdir;
my $file_list;
my $genomedir;
my $junctionsdir;
my $splitdir;
my $projectid;
my $readlength;
my $mismatches;
my $threshold=1;
my $staggered;

# Read command line options
GetOptions(
    'threshold|t=s' => \$threshold,
    'staggered' => \$staggered
    );

# Read the options file
my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$file_list=$options{'FILELIST'};
$genomedir=$options{'GENOMEDIR'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$splitdir=$options{'SPLITMAPDIR'};
$projectid=$options{'PROJECTID'};
$readlength=$options{'READLENGTH'};
$mismatches=$options{'MISMATCHES'};

# Get the files
my %files=%{read_file_list($file_list)};
my %lanes_paired=%{get_lanes_paired(\%files)};
my %lanes_single=%{get_lanes_single(\%files)};

# First read in the gemome files and cluster them
cluster_genome_files(\%lanes_paired,
		     $genomedir,
		     $tmpdir,
		     $projectid,
		     $staggered,
    );

# Read junction files and cluster them
cluster_split_files(\%lanes_single,
		    $junctionsdir,
		    $tmpdir,
		    $projectid,
		    $staggered,
		    'junctions'
    );

# Read split files and cluster them
cluster_split_files(\%lanes_single,
		    $splitdir,
		    $tmpdir,
		    $projectid,
		    $staggered,
		    'split'
    );

# combine the three generated temporary files into the final one

exit;

sub cluster_split_files {
    my $lanes=shift;
    my $splitdir=shift;
    my $tmpdir=shift;
    my $project=shift;
    my $stagger=shift;
    my $splittype=shift;

    my $type;
    my @uniquefiles;
    my @multifiles;

    print STDERR "Clustering split reads...\n";

    # Collect the necessary files
    foreach my $lane (keys %{$lanes}) {
	my $type='single';

	my $infilename1=$splitdir.'/'.$lane.'.'.$type.'.unique.gtf.gz';
	my $infilename2=$splitdir.'/'.$lane.'.'.$type.'.multi.gtf.gz';

	unless (-r $infilename1) {die "$infilename1 is not readable\n";}
	print STDERR "Collecting $infilename1\n";
	unless (-r $infilename2) {die "$infilename2 is not readable\n";}
	print STDERR "Collecting $infilename2\n";

	push @uniquefiles, $infilename1;
	push @multifiles, $infilename2;
    }

    # Cluster the unique reads
    my $uniqueclusterfn=$projectid.'.unique.'.$splittype.'.cluster.gtf.gz';
    my $uniqueinfh=get_sorted_gff_fh(\@uniquefiles,
				     $stagger,
				     $tmpdir);
    cluster_gff($uniqueinfh,
		$uniqueclusterfn,
		$threshold,
	);
    close($uniqueinfh);

    # Cluster the multireads
    my $multiclusterfn=$projectid.'.multi.'.$splittype.'.cluster.gtf.gz';
    my $multiinfh=get_sorted_gff_fh(\@multifiles,
				    $stagger,
				    $tmpdir);
    cluster_gff($multiinfh,
		$multiclusterfn,
		$threshold,
	);
    close($multiinfh);

    print STDERR "done\n";
}

sub cluster_genome_files {
    my $lanes=shift;
    my $genomedir=shift;
    my $tmpdir=shift;
    my $project=shift;
    my $stagger=shift;

    my $type;
    my @uniquefiles;
    my @multifiles;

   print STDERR "Clustering genome reads...\n";

    # Collect the necessary files
    foreach my $lane (keys %{$lanes}) {
	my $type;
	if ($lanes->{$lane} == 1) {
	    $type='single';
	} elsif ($lanes->{$lane} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type\n";
	}

	my $infilename1=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.gz';
	my $infilename2=$genomedir.'/'.$lane.'.'.$type.'.multi.gtf.gz';

	unless (-r $infilename1) {die "$infilename1 is not readable\n";}
	print STDERR "Collecting $infilename1\n";
	unless (-r $infilename2) {die "$infilename2 is not readable\n";}
	print STDERR "Collecting $infilename2\n";

	push @uniquefiles, $infilename1;
	push @multifiles, $infilename2;
    }

    # Cluster the unique reads
    my $uniqueclusterfn=$projectid.'.unique.genome.cluster.gtf.gz';
    my $uniqueinfh=get_sorted_gff_fh(\@uniquefiles,
				     $stagger,
				     $tmpdir);
    cluster_gff($uniqueinfh,
		$uniqueclusterfn,
		$threshold,
		);
    close($uniqueinfh);
    
    # Cluster the multireads
    my $multiclusterfn=$projectid.'.multi.genome.cluster.gtf.gz';
    my $multiinfh=get_sorted_gff_fh(\@multifiles,
				    $stagger,
				    $tmpdir);
    cluster_gff($multiinfh,
		$multiclusterfn,
		$threshold,
		);
    close($multiinfh);

    print STDERR "done\n";
}

sub get_lanes_single {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[1]}++;
    }

    return(\%lanes);
}

sub get_lanes_paired {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}++;
    }

    return(\%lanes);
}
