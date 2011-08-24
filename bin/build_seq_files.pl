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
# This script will take as an input a list of bam files and it generate for
# each one of them the equivalent fastq file.

# Load some modules
use RNAseq_pipeline3 qw(get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list print_table_file);
use Bio::SeqIO;
use Tools::Bam qw(bam2sequence);

# Get some options from the configuration file
my $prefix;
my $tmpdir;
my $debug=1;

my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$tmpdir=$options{'LOCALDIR'};

# Get the files we are going to process
my %files=%{read_file_list()};

# Get a log file
my $log_fh=get_log_fh('build_seq_files.RNAseq.log',
		      $debug);

foreach my $infn (keys %files) {
    print $log_fh "Processing $infn\n";

    my $infile=$tmpdir.'/'.$infn;

    # Check the file type
    unless ($infile=~/.bam$/) {
	warn "File $infile does not seem to be bam from the name I'll try to extrac it anyway, but this should not happen\n";
	next;
    }

    # Process the file
    my $outfn=bam2sequence($infile);
    print join("\t",
	       $infn,
	       $outfn),"\n";
}

# Print out the table file
my $outtable=$prefix.'_seq_files';
my $outtablesql="CREATE TABLE $outtable (
    BAM_file varchar(50) not null,
    Seq_file varchar(100) not null
    );";
print_table_file($outtable,
		 $outtablesql);

close($log_fh);

exit;
