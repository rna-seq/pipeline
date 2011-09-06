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
# This script will take preprocessiong script as input
# the  preprocessing script should take as an argument a single 
# fastq/fasta file, and must be in the path, preferably copied to the bin dir.
# The preprocessing script will output the preprocessed file to the STDOUT.
# This script will run the preprocessing script on each of the files in the
# readData directory and copy the zipped results back to the readDirectory
# overwritting the original files (which will be mover to sequences currently,
# but maybe not always)

# Load some modules
use RNAseq_pipeline3 qw(get_fh run_system_command get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use Bio::SeqIO;

# Get some command line options and set their defaults
my $species;
my $project;
my $sizethreshold=2; # This is the maximum allowe filesize if it is higher
                     # the script will use tmp space on the disk
my $debug=1;
my $localdir;
my $projdir;
my $readdir;
my $script;
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};
$species=$options{'SPECIES'};
$localdir=$options{'LOCALDIR'};
$projdir=$options{'PROJECT'};
$readdir=$projdir.'/'.$options{'READDIR'};
$script=shift || $options{'PREPROCESS'};

# If mismatches is set, the reads will be filtered

# Paired will indicate paired end reads. This is only appliable to the initial
# case where the file is a txt file, as both pairs will appear as the first and
# second halves of the read.

# Check the localdir
if ($localdir &&
    -e $localdir) {
    $localdir=~s/\/$//;
    print $log_fh "Using $localdir for tmp files\n";
} else {
    die "Could not find $localdir\n";
}

# Process the input
my ($scriptpath)=split(/\s+/,$script);
unless ($script &&
	-e $scriptpath) {
    die "Could not find input script $script at $scriptpath\n";
}

# Get the files we are going to process
my %files=%{read_file_list()};

# Get a log file
my $log_fh=get_log_fh('preprocess.RNAseq.log',
		      $debug);

foreach my $infile (keys %files) {
    print $log_fh "Processing $infile\n";

    my $outfile=$localdir.'/'.$infile.'.out.gz';
    my $command="$script $readdir/$infile.gz |gzip -9 -c > $outfile";
    run_system_command($command,
		       $log_fh);

    # save the original read files
    my $original=$readdir.'/'.$infile.'.gz';
    my $destination=$projdir.'/sequence/'.$infile.'.gz';
    $command="mv $original $destination";
    run_system_command($command,
		       $log_fh);

    # Put the preprocessed file into the directory
    $command="mv $outfile $readdir/$infile.gz";
    run_system_command($command,
		       $log_fh);
}
close($log_fh);

exit;
