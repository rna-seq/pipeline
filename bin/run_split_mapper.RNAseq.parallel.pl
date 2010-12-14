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
# This script will take the necessary commands for running the split mapper,
# This is unfinished
use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file);
use RNAseq_GEM3 ('check_index','get_split_mapper_routines',
		 'check_input','check_index');

my $index;
my $mapper;
my $infile;
my $outfile;
my $tempdir;
my $filetype;
my $mismatches;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'infile|i=s' => \$infile,
	   'outfile|o=s' => \$outfile,
	   'mismatches|m=i' => \$mismatches);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$tempdir=$options{'LOCALPARALLEL'};
unless ($mismatches) {
    $mismatches=$options{'MISMATCHES'};
}

# Make sure we have all parameters
unless ($index &&
	$infile) {
    die "At least an index file and an input file are required\n";
}

unless ($outfile) {
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

if (-r "$outfile.split-map" ||
    -r "$outfile.split-map.gz") {
    print STDERR "$outfile Exists.. Skipping\n";
    exit;
}

if ($filetype eq 'fastq') {
    # The readfile has qualities
    die "Split mapper cannot use qualities only fasta files\n";
}

my %mapper=%{get_split_mapper_routines()};

# Run the mapper
$mapper=uc($mapper);
unless (exists $mapper{$mapper}) {
    die "I am not familiar with $mapper. This split-mapper is not implemented\n";
}

# Map the file
$mapper{$mapper}->($index,
		   $infile,
		   $outfile,
		   $tempdir,
		   $mismatches,
		   1);

exit;

