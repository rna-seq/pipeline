#!/soft/bin/perl

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

#    Author : David Gonzalez, david.gonzalez@crg.eu

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
use RNAseq_GEM3 ('check_index','get_mapper_routines',
		 'check_input');

my $index;
my $mapper;
my $infile;
my $outfile;
my $qualities;
my $threads;
my $ignore_quals;
my $tempdir;
my $filetype;
my $mismatches;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'infile|i=s' => \$infile,
	   'outfile|o=s' => \$outfile,
	   'noqual' => \$ignore_quals,
	   'mismatches=i' => \$mismatches,
	   'threads|t=i' => \$threads);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$tempdir=$options{'LOCALPARALLEL'};
$qualities=$options{'QUALITIES'};
unless ($threads) {
    $threads=$options{'THREADS'};
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

if (-r "$outfile.map" ||
    -r "$outfile.map.gz") {
    print STDERR "$outfile Exists.. Skipping\n";
    exit;
}

if ($ignore_quals) {
    print STDERR "Ignoring qualities\n";
    $qualities='ignore';
}

my %mapper=%{get_mapper_routines()};

# Run the mapper
$mapper=uc($mapper);
unless (exists $mapper{$mapper}) {
    die "I am not familiar with $mapper. The mapper is not implemented\n";
}

# Map the file
$mapper{$mapper}->($index,
		   $infile,
		   $outfile,
		   $threads,
    		   $qualities,
		   $tempdir,
		   $mismatches,
		   1);

exit;

