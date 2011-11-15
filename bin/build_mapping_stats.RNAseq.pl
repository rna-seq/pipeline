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
# This script should build for a mapped file the following stats:
# Total number of reads, number of reads mapped, number of unique mappings
# within the mismatch threshold and number of 1:0:0 matches
# it will also use R to build a graphical summary with this info.

### TO DO
# Modify, so the list of files does not have to be supplied on the command line

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use RNAseq_GEM3 qw(parse_gem_line);

my $mismatches;
my $file_list;

my %options=%{read_config_file()};
$mismatches=$options{'MISMATCHES'};
$file_list=$options{'FILELIST'};

my %files=%{read_file_list()};

unless (@ARGV) {
    die "No input files\n";
}

# This should give us the equivalent files in the same order in both lists
my @input=sort @ARGV;

# We will check if there are repeated files in the input (one zipped and one
# unzipped), if so we will use the unzipped one.
my %input;
foreach my $file (@input) {
    my $unzipped=$file;
    $unzipped=~s/.gz$//;
    if (exists $input{$unzipped}) {
	print STDERR "WARNING: An unzipped version of $file is present so $file will NOT be used\n";
    } else {
	$input{$file}=1;
    }
}

@input=sort keys %input;
my @files=sort keys %files;
unless (@input == @files) {
    die "We have a filenumber problem??\n";
}

for (my $i=0;$i<@input;$i++) {
    my $infile=$input[$i];
    my $shortname=$files{$files[$i]}->[1];
    my $infh=get_fh($infile);
    my $total=0;
    my $mapped=0;
    my $unique=0;
    my $unique100=0;
    
    while (my $line=<$infh>) {
	if ($line=~/^\s*$/) {
	    warn "Empty line in $infile, maybe something is wrong...\n";
	    next;
	}
	my %line=%{parse_gem_line($line)};
	$total++;
	
	if ($line{'matches'}=~/^0(:0)*$/) {
	    # If the read does not map we needn't continue
	    next;
	} else {
	    # This means the read maps somewhere
	    $mapped++;
	}
	
	# Check if the read is unique
	my @hits=split(':',$line{'matches'});

	if (@hits < $mismatches - 1) {
	    $mismatches= @hits - 1;
	    warn "File was mapped with $mismatches mismatches, reducing mismatch threshold accordingly\n";
	}

	for (my $i=0;$i<=$mismatches;$i++) {
	    if ($hits[$i] eq '-') {
		last;
	    } elsif ($hits[$i] eq '!') {
		last;
	    } elsif ($hits[$i] > 0) {
		# This is the best hit
		if ($hits[$i] == 1) {
		    $unique++;
		}
		last;
	    }
	}	 
	
	# Check if the read is 1:0:0
	if ($line{'matches'}=~/^1(:0)+$/) {
	    $unique100++;
	}
	
    }
    close($infh);

    $infile=~s/.*\///;
    print join("\t",
	       $infile,
	       $total,
	       $mapped,
	       $unique,
	       $unique100,
	       $shortname),"\n";
}

exit;

