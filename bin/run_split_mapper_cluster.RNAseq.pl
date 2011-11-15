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

### TO DO
# Make the mapper take a list of files

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list send2cluster);
use RNAseq_GEM3 ('check_index','determine_quality_type',
		 'get_split_mapper_routines','check_input');

my $index;
my $mapper;
my $outdir;
my $threads;
my $tempdir;
my $filetype;
my $mismatches;
my $paralleldir;
my $readdata;
my $projdir;
my $bindir;
my $splitdir;
my $queue;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'outdir|o=s' => \$outdir);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$threads=$options{'THREADS'};
$tempdir=$options{'LOCALDIR'};
$mismatches=$options{'MISMATCHES'};
$paralleldir=$options{'LOCALPARALLEL'};
$readdata=$options{'READDIR'};
$projdir=$options{'PROJECT'};
$bindir=$options{'BIN'};
$splitdir=$options{'SPLITMAPDIR'};
$queue=$options{'CLUSTER'};

# Decide where to put the output
unless ($outdir) {
    die "Where should I put the results ???\n";
}
$outdir=~s/\/$//;

# Make sure we have a valid index
unless ($index) {
    die "An index file is required\n";
}

# Decide if we are on the cluster or not
my $usecluster=0;
if ($options{'CLUSTER'}) {
    print STDERR "Mapping in the cluster\n";
    $usecluster=1;
} else {
    print STDERR "No cluster name provided. Mapping locally\n";
}

# Make sure we have a valid index
unless ($index) {
    die "An index file is required\n";
}

# check for the existence of the index
my $index_ok=check_index($index);
unless ($index_ok) {
    die $index,"\tIs not a valid index\n";
}

# Get the list of files we need to map
my %files=%{unmapped_file_list()};

# Check each of the files and determine if it is present in the paralleldir
# directory if not unzip it and copy it there.
my %locations=%{check_unmapped_files(\%files,
				     $splitdir)};

my @filepairs;
foreach my $file (keys %locations) {
    my $qualities;
    # Check the input file
    my $infile=$locations{$file};
    my $input_ok=check_input($infile,
			     \$filetype);
    unless ($input_ok) {
	die $infile,"\tIs not readable\n";
    }

    # Make the outfile
    my $outfile=$infile;
    $outfile=~s/.*\///;
    $outfile=$outdir.'/'.$outfile;
    $outfile=~s/(.fa(stq)*)$/.gem/;

    if (-r "$outfile.split-map" ||
	-r "$outfile.split-map.gz") {
	print STDERR "$outfile Exists.. Skipping\n";
	next;
    }

    if ($filetype eq 'fastq') {
	# The readfile has qualities
	die "I cannot split-map using qualities\n";
    }

    push @filepairs,[$infile,$outfile,$qualities];

}

unless (@filepairs) {
    print STDOUT "Everything seems to be mapped already\n";
    exit;
}

if ($usecluster) {
    # Build the submission file
    my $jobname='RNAseqSplitMap';
    my $subfile=build_run_mapper_submission(\@filepairs,
					    $bindir,
					    $index,
					    $jobname);
    send2cluster($subfile,
		 $queue,
		 $jobname);

    # clean up
    my $command="rm $subfile";
    run_system_command($command);

    } else {
	my %mapper=%{get_split_mapper_routines()};
	foreach my $pair (@filepairs) {
	    # Map the file
	    $mapper{$mapper}->($index,
			       $pair->[0],
			       $pair->[1],
			       $tempdir,
			       $mismatches,
			       1);
	}
}
exit;

# This subroutine will get a list of the unmapped files that should be present
# based on the contents of the read.list.txt file
sub unmapped_file_list {
    my %files;
    my $file_list=shift;
    # Set a default for the configuration file
    unless ($file_list) {
	my %options=%{read_config_file()};
	$file_list=$options{'FILELIST'};
    }
    my $read_fh=get_fh($file_list);

    while (my $line=<$read_fh>) {
	chomp($line);
	my @line=split("\t",$line);
	if (@line < 3 ) {
	    die "Number of elements in $line is smaller than 3\n";
	}
	my $file=$line[2].'.unmapped.fa';
	$files{$file}=1;
    }

    close($file_list);

    return(\%files);
}

sub build_run_mapper_submission {
    my $pairs=shift;
    my $bidir=shift;
    my $index=shift;
    my $jobname=shift;

    print STDERR 'Building submission file...';
    my $filenum=@{$pairs};
     
    unless(@{$pairs}) {
	die "No input supplied\n";
    }
    
    # Get the input and output files
    my @infiles;
    my @outfiles;
    foreach my $pair (@{$pairs}) {
	push @infiles,$pair->[0];
	push @outfiles,$pair->[1];
    }

    # Print the submission file
    my $subfile="subfile.$$.job";
    my $outfh=get_fh($subfile,1);
    
    print $outfh <<FORMEND;
# Get the job name
#\$ -N $jobname
    
# Set the array jobs
#\$ -t 1-$filenum

# Request 8 cpus this cannot be done, but we can request memmory
#\$ -l h_vmem=16G

# Write in to the current working directory
#\$ -cwd 
export PATH=\$PATH:/soft/bin
# Make sure the sorting order is "a la C"
export LC_ALL=C

infiles=(@infiles)
outfiles=(@outfiles)

export infile=\${infiles[\$SGE_TASK_ID-1]}
export outfile=\${outfiles[\$SGE_TASK_ID-1]}

echo \$HOSTNAME >&2
$bindir/run_split_mapper.RNAseq.parallel.pl -index $index -infile \$infile -outfile \$outfile > \$infile.mapping.log
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile);
}

sub check_unmapped_files {
    my $files=shift;
    my $paralleldir=shift;

    my %locations;

    foreach my $file (keys %{$files}) {
	print STDERR "Checking $file...";
	my $filepath1=$paralleldir.'/'.$file;
	if ( -r $filepath1) {
	    print STDERR "Present in $paralleldir\n";
	    $locations{$file}=$filepath1;
	} else {
	    die "I can't find file $file\n";
	}
    }

    return(\%locations);
}
