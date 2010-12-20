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
# This script will take the necessary commands for running the mapper, and it
# will check if the input files have qualities or not as well as deciding the
# type of qualities at some stage.

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file);
use RNAseq_GEM3 ('check_index','determine_quality_type','split_map',
		 'check_input');

my $index;
my $mapper;
my $infile;
my $outfile;
my $outdir;
my $threads;
my $tempdir;
my $filetype;
my $mismatches;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'infile|i=s' => \$infile,
	   'outdir|o=s' => \$outdir);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$threads=$options{'THREADS'};
$tempdir=$options{'LOCALDIR'};
$mismatches=$options{'MISMATCHES'};

# Make sure we have all parameters
unless ($index &&
	$infile) {
    die "At least an index file and an input file are required\n";
}

unless ($outdir) {
    die "Where should I put the results ???\n";
}
# check for the existence of the index
my $index_ok=check_index($index);
unless ($index_ok) {
    die $index,"\tIs not a valid index\n";
}

# Check the input file
my $input_ok=check_input($infile,
			 \$filetype);
unless ($input_ok) {
    die $infile,"\tIs not readable\n";
}

# Make the outfile
$outfile=$infile;
$outfile=~s/.*\///;
$outfile=$outdir.'/'.$outfile;
if ($mapper eq 'GEM') {
    $outfile=~s/(.fa(stq)*)$/.gem/;
} else {
   $outfile=~s/(.fa(stq)*)$/.mapper/;
}

if (-r "$outfile.map" ||
    -r "$outfile.map.gz") {
    print STDERR "$outfile Exists.. Skipping\n";
    exit;
}

# Run the split mapper
$mapper{$mapper}->($index,
		   $infile,
		   $outfile,
		   $threads,
    		   $qualities,
		   $tempdir,
		   $mismatches,
		   1);

exit;

