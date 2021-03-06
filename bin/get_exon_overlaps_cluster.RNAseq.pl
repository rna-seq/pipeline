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
# This script should run overlap on a set of annotated exons

# Load some modules
use POSIX qw(:signal_h :errno_h :sys_wait_h);
use RNAseq_pipeline3 ('get_fh','parse_gff_line');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_feature_overlap_sub');
use Getopt::Long;

# Declare some variables
my $threshold=1;
my $tmpdir;
my $clusterdir;
my $genomedir;
my $exondir;
my $splitdir;
my $junctionsdir;
my $projectid;
my $prefix;
my $stranded;
my $file_list;
my $parallel='cluster';
my $paralleltmp;
my $bindir;

GetOptions('threshold|t=i' => \$threshold,
    );

# Read the options file
my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$genomedir=$options{'GENOMEDIR'};
$exondir=$options{'PROJECT'}.'/'.$options{'EXONDIR'};
$splitdir=$options{'SPLITMAPDIR'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$projectid=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$stranded=$options{'STRANDED'};
$file_list=$options{'FILELIST'};
$paralleltmp=$options{'PROJECT'}.'/'.$options{'LOCALPARALLEL'};
$bindir=$options{'BIN'};
unless($options{'CLUSTER'}) {
    print STDERR "Running locally, as no cluster has been defined\n";
    $parallel='default';
}

unless ($paralleltmp) {
    $paralleltmp=$options{'PROJECT'}.'/work';
}

# Get the required sub
*get_feature_overlap=get_feature_overlap_sub($parallel,
					     $paralleltmp,
					     $bindir,
					     $options{'CLUSTER'});

# First get the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get for each lane the files corresponding to the genome mapping
my %lane_files=%{get_lane_files(\%files,
				$genomedir)};

my $annotation=$exondir.'/'.$prefix.'.exon.gtf';

my %coverage;

foreach my $pair (keys %lane_files) {
    # Find the coverage granted by genome mapping;
    my $outfile=$lane_files{$pair};
    $outfile=~s/genome/exons/;
    $outfile=~s/.gz$//;
    
    # Check if file exists already
    my $overfinalname=$outfile;
    $overfinalname=~s/.gz//;
    $overfinalname.='.overlap.gz';
    if (-r $overfinalname) {
	print STDERR $overfinalname,"\tAlready exists. Skipping\n";
	next;
    } else {       
	get_feature_overlap($lane_files{$pair},
			    $annotation,
			    $stranded,
			    $outfile,
			    $tmpdir);
    }
}
    
exit;

sub get_lane_files {
    my $files=shift;
    my $dir=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	$lane_files{$pair}{$lane}=1;
    }

    # Determine if we are looking at paired or single reads
    foreach my $pair (keys %lane_files) {
	my @lanes=keys %{$lane_files{$pair}};
	if (@lanes == 1) {
	    # single files
	    print STDERR "$pair set identified as single reads\n";
	    $lane_files{$pair}=$dir.'/'.$lanes[0].'.single.unique.gtf.gz';
	} elsif (@lanes == 2) {
	    # paired files
	    print STDERR "$pair set identified as paired reads\n";
	    $lane_files{$pair}=$dir.'/'.$pair.'.paired.unique.gtf.gz';
	} else {
	    # This should never happen
	    die "Apparently there are more than two files grouped\n";
	}
    }

    return(\%lane_files);
}
