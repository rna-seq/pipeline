#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take those reads that support the fusion events from the
# genome and the transcriptome unique gff files and it will cluster them
# together.

# Load some modules
use Getopt::Long;
use RNAseq_pipeline3 ('get_fh','parse_gff_line',
		      'cluster_gff','get_gff_from_junc_id',
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
my $transdir;

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
$transdir=$options{'TRANSDIR'};

# Get the files
my %files=%{read_file_list($file_list)};
my %lanes_paired=%{get_lanes_paired(\%files)};
my %lanes_single=%{get_lanes_single(\%files)};

# Get the genome read files
my %groups;
foreach my $lane (keys %lanes_paired) {
    my $type;
    if ($lanes_paired{$lane} == 1) {
	$type='single';
    } elsif ($lanes_paired{$lane} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }
    
    my $fusionfile=$transdir.'/'.$lane.'.'.$type.'.trans.fusions.txt.gz';

    if (-r $fusionfile) {
	print STDERR "Processing $fusionfile";
    } else {
	print STDERR "I can't find $fusionfile\n";
    }

    # Get the necessary reads from the file and build a genome and transcriptome
    # gff file with only those reads
    my %reads=%{get_reads_from_file($fusionfile)};
    my $genomefile=build_genome_file(\%reads,
				     $lane,
				     $genomedir,
				     $type);
    my $transfile=build_genome_file(\%reads,
				    $lane,
				    $transdir,
				    $type);
    my $clusteredgen=cluster_files($genomefile,
				   $threshold,
				   $staggered,
				   'genome',
				   $genomedir);
    my $clusteredtran=cluster_files($transfile,
				    $threshold,
				    $staggered,
				    'trans',
				    $transdir);

    $groups{$lane}=[$genomefile,$transfile,
		    $clusteredgen,$clusteredtran];

}

exit;

sub build_genome_file {
    my $reads=shift;
    my $lane=shift;
    my $genomedir=shift;
    my $type=shift;

    my $genomefile=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.gz';
    my $outfn=$genomedir.'/'.$lane.'.'.$type.'.fusion.supporting.gtf.gz';
    print STDERR "Extracting reads from $genomefile\n";

    my $genfh=get_fh($genomefile);
    my $outfh=get_fh($outfn,1);
    while (my $line=<$genfh>) {
	my %line=%{parse_gff_line($line)};
	my $read_id=$line{'feature'}{'read_id'};
	$read_id=~s/\|?p?[12]$//;
	$read_id=~s/\/$//;
	if ($reads->{$read_id}) {
	    print $outfh $line;
	}
    }
    close($outfh);
    close($genfh);

    return($outfn);
}

sub get_reads_from_file {
    my $file=shift;

    my %reads;
    my $transfh=get_fh($file);

    while (my $line=<$transfh>) {
	chomp($line);
	my @line=split("\t",$line);
	if ($line[7] eq $line[13]) {
	    next;
	}
	my $read_id=$line[0];
	$read_id=~s/\|$//;
	$reads{$read_id}++;
    }
    close($transfh);

    my $count=keys %reads;
    print STDERR $count,"\tFusion supporting reads found in $file\n";

    return(\%reads);
}

sub cluster_files {
    my $infile=shift;
    my $threshold=shift;
    my $stagger=shift;
    my $type=shift;
    my $dir=shift;

    # Cluster the unique reads
    my $uniqueclusterfn=$dir.'/'.$projectid.'.unique.'.$type.'.cluster.gtf.gz';
    print STDERR "Clustering $infile in $uniqueclusterfn\n";

    my $uniqueinfh=get_sorted_gff_fh([$infile],
				     $stagger,
				     $tmpdir);
    cluster_gff($uniqueinfh,
		$uniqueclusterfn,
		$threshold,
		);
    close($uniqueinfh);
    return($uniqueclusterfn);    
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
