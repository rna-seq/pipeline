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
# This script should read the file_list file, and it will take each of the
# read files and copy it unzipped to the local directory that has been
# specified. This should avoid problems with compressed or uncompressed files

# Load some modules
use RNAseq_pipeline3 qw(get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

# Get some options from the configuration file
my $tmpdir;
my $readdir;
my $file_list;
my $debug;

my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$readdir=$options{'READDIR'};

my $log_fh=get_log_fh('prepare_files.RNAseq.log',
		      $debug);
my %files=%{read_file_list()};

foreach my $readfile (keys %files) {
    print $log_fh "Processing $readfile\n";
    my $command;

    # Set source and target
    my $infile=$readdir.'/'.$readfile;
    my $target=$tmpdir.'/'.$readfile;

    # Check if target is already present and readable
    if (-e $target && 
	-r $target) {
	print $log_fh "$readfile already present in $tmpdir\n";
	next;
    }

    if (-r $infile) {
	print $log_fh $infile,"\tIs unzipped gzipping for storage\n";
	$command="cp $infile $target; gzip -7 $infile";
    } elsif (-r $infile.'.gz') {
	print $log_fh $infile,"\tIs gzipped. Inflating...\n";
	$command="gunzip $infile.gz -c > $target";
    } else {
	die "I can't find $infile(.gz)\n";
    }

    print $log_fh "Executing: $command\n";
    system($command);
}
close($log_fh);

exit;
