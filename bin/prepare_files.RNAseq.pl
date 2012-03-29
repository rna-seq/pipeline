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
# This script should read the file_list file, and it will take each of the
# read files and copy it unzipped to the local directory that has been
# specified. This should avoid problems with compressed or uncompressed files

# Load some modules
use RNAseq_pipeline3 qw(get_log_fh run_system_command);
use GRAPE::Logs;
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

# Get some options from the configuration file
my $tmpdir;
my $readdir;
my $file_list;
my $debug;

my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$readdir=$options{'READDIR'};
$debug=$options{'DEBUG'};

my $logobj=GRAPE::Logs->new('prepare_files.RNAseq.log',
			    $debug);
my %files=%{read_file_list()};

foreach my $readfile (keys %files) {
    $logobj->printlog("Processing $readfile");
    my $command;

    # Set source and target
    my $infile=$readdir.'/'.$readfile;
    my $target=$tmpdir.'/'.$readfile;

    # Check if the necessary directory is present
    unless (-e $tmpdir) {
	die "WARNING: $tmpdir is not present. Please check settings or rerun start\n";
    }

    # Check if target is already present and readable
    if (-e $target && 
	-r $target) {
	$logobj->printlog("$readfile already present in $tmpdir");
	next;
    }

    if (-r $infile) {
	if ($infile=~/.bam$/) {
	    $logobj->printlog($infile,"\tIs in BAM format. Copying to $tmpdir");
	    $command="cp $infile $target";
	} else {
	    $logobj->printlog($infile,"\tIs unzipped gzipping for storage");
	    $command="cp $infile $target; gzip -7 $infile";
	}
    } elsif (-r $infile.'.gz') {
	$logobj->printlog($infile,"\tIs gzipped. Inflating...");
	$command="gunzip $infile.gz -c > $target";
    } else {
	die "I can't find $infile(.gz)\n";
    }

    run_system_command($command);
}

exit;
