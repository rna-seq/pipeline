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
# This script will calculate the md5sum of the fastq files uncompressed
# and compressed. It will also calculate an additional md5sum index which will
# be the one for the pipeline run. This will be md5sum of the sorted
# concatenation of all the md5sums of the unzipped files corresponding to the
# experiment.

# Load some modules
use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'set_exp_info_sub');
use Bio::SeqIO;

# Get some options from the configuration file
my $project;
my $tmpdir;
my $readdir;
my $samdir;
my $file_list;
my $prefix;
my $debug=0;
my $proj_id;
my $exp_id;

my %options=%{read_config_file()};
$project=$options{'PROJECT'};
$tmpdir=$options{'LOCALDIR'};
$readdir=$options{'READDIR'};
$samdir=$options{'PROJECT'}.'/SAM';
$file_list=$options{'FILELIST'};
$prefix=$options{'PREFIX'};
$proj_id=$options{'PROJECTID'};
$exp_id=$options{'EXPID'};

# Get some subs
*set_exp_info=set_exp_info_sub('md5sum');

# Get the files we are going to process
my %files=%{read_file_list()};

# Get a log file
my $log_fh=get_log_fh('build_store_reads.RNAseq.log',
		      $debug);

# Initialize a hash for the checksums
my %sums;
my %sumsunzipped;

foreach my $infile (keys %files) {
    print $log_fh "Processing $infile...";

    my $lanename=$files{$infile}->[1];

    # Process the file
    my $filepath=$tmpdir.'/'.$infile;
    # Calculate the md5sum of the unzipped files
    if (-r $filepath) {
	my ($file,$sum)=get_md5sum($filepath);
	$sums{$infile}{'unzipped'}=[$file,$sum];
	$sumsunzipped{$sum}++;
    } else {
	die "Cannot read $infile\n";
    }

    # Calculate the md5sum of the zipped files
    my $readpath=$readdir.'/'.$infile.'.gz';
    if (-r $readpath) {
	my ($file,$sum)=get_md5sum($readpath);
	$sums{$infile}{'zipped'}=[$file,$sum];
    } else {
	die "Cannot read $infile\n";
    }

    print $log_fh "done\n";
}

# Get the final md5sum
my $finalsum=get_global_md5sum(\%sumsunzipped);
set_exp_info($proj_id,
	     $exp_id,
	     $finalsum);

# print out the results
foreach my $file (keys %sums) {
    my @line;
    push @line, @{$sums{$file}{'unzipped'}};
    push @line, @{$sums{$file}{'zipped'}};
    push @line, $finalsum;
    print join("\t",
	       $file,
	       @line),"\n";
}

exit;

sub get_md5sum {
    my $filepath=shift;
    my $filename=$filepath;
    $filename=~s/.+\///;

    my $sum=`md5sum $filepath`;
    my ($md5sum)=split(' ',$sum);
    
    return($filename,$md5sum);
}

sub get_global_md5sum {
    my $sums=shift;

    my %md5sums;
    my $outfile=$$.'.md5sum.txt';
    my $outfh=get_fh($outfile,1);
    foreach my $sum (sort keys %{$sums}) {
	print $outfh $sum,"\n";
    }
    close($outfh);

    my $sum=`md5sum $outfile`;
    my ($md5sum)=split(' ',$sum);
    my $command="rm $outfile";
    run_system_command($command);

    return($md5sum);
}

