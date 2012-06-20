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
# This script will run cufflinks on the merged SAM files with the most basic
# settings

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

my $prefix;
my $bindir;
my $tmpdir;
my $samdir;
my $debug=0;

my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$bindir=$options{'BIN'};
$tmpdir=$options{'LOCALDIR'};
$samdir=$options{'SAMDIR'};

# Get a log file
my $log_fh=get_log_fh('run_cufflinks.RNAseq.log',
		      $debug);

# First get a list of the bam files we are going to process
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# Run cufflinks for each of the files that we have and move the output to
# a file with the adequate name
foreach my $lane (keys %lanes) {
    my %reads;

    # Check if the file corresponding to this lane exists already
    my $outfile=$samdir.'/'.$lane.'.transcripts.gtf';
    if (-s $outfile) {
	print STDERR $lane,"\tTranscripts present. Skipping...\n";
	next;
    }

    my $outdir=$tmpdir.'/'.$lane;

    # Get the files corresponding to the both halves of the reads
    my $infilename=$samdir.'/'.$lane.'.merged.bam';
    if (-r $infilename) {
	print $log_fh "Processing $infilename\n";
    } else {
	die "Can't read $infilename\n";
    }

    my %results;
    # Run the flux capacitor
    # This is the parallel step
    my $cufffile=run_cufflinks($infilename,
			       $outfile,
			       $outdir);

    # Move the file to the correct location
    my $command="mv $cufffile $outfile";
    run_system_command($command);
}

exit;

sub run_cufflinks {
    my $infile=shift;
    my $outfile=shift;
    my $outdir=shift;

    my $command='cufflinks --library-type fr-unstranded ';
    $command.="-o $outdir $infile";
    run_system_command($command);

    my $cufffile=$outdir.'/'.'transcripts.gtf';

    return($cufffile);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	push @{$lanes{$files->{$file}->[0]}},$files->{$file}->[1];
    }

    return(\%lanes);
}
