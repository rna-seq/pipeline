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
use RNAseq_GEM3 ('check_index','determine_quality_type','get_mapper_routines',
		 'check_input');
use Tools::Bam ('bam2sequence');

my $index;
my $mappermem;
my $mapper;
my $outdir;
my $qualities;
my $threads;
my $ignore_quals;
my $tempdir;
my $filetype;
my $mismatches;
my $paralleldir;
my $readdata;
my $projdir;
my $bindir;
my $zip=1;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'outdir|o=s' => \$outdir,
	   'noqual' => \$ignore_quals,
	   'zip' => \$zip);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$mappermem=$options{'MAPPERMEM'} || '5G';
$threads=$options{'THREADS'};
$tempdir=$options{'LOCALDIR'};
$mismatches=$options{'MISMATCHES'};
$paralleldir=$options{'LOCALPARALLEL'};
$readdata=$options{'READDIR'};
$projdir=$options{'PROJECT'};
$bindir=$options{'BIN'};
$qualities=$options{'QUALITIES'};
$threads=$options{'THREADS'};

# Decide where to put the output
unless ($outdir) {
    die "Where should I put the results ???\n";
}
$outdir=~s/\/$//;
my $usecluster=0;

if ($options{'CLUSTER'}) {
    print STDERR "Mapping in the cluster with $threads threads\n";
    $usecluster=1;
} else {
    print STDERR "I can only use the cluster if I am introduced to it, and I have no idea of its name. Mapping locally\n";
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
my %files=%{read_file_list()};

# Check each of the files and determine if it is present in the paralleldir
# directory if not unzip it and copy it there.
my %locations=%{check_read_files(\%files,
				 $paralleldir,
				 $tempdir,
				 $readdata)};

my @filepairs;
foreach my $file (keys %locations) {
    # Check the input file
    my $infile=$locations{$file};

    # If the input is in BAM format we need to make sure the fastq file with the
    # sequence is available. For this we will create a fastq file changing the
    # name .bam for .fastq of .fa, but we will keep the original name for the
    # purpose of naming the output.

    # Make the outfile
    my $outfile=$infile;
    $outfile=~s/.*\///;
    $outfile=$outdir.'/'.$outfile;
    $outfile=~s/(.fa(stq)*)$/.gem/;

    if ($infile=~/.bam$/) {
	$infile=bam2sequence($infile);
    }

    my $input_ok=check_input($infile,
			     \$filetype);

    unless ($input_ok) {
	die $infile,"\tIs not readable\n";
    }

    if (-r "$outfile.map" ||
	-r "$outfile.map.gz") {
	print STDERR "$outfile Exists.. Skipping\n";
	next;
    }

    if ($ignore_quals) {
	print STDERR "Ignoring qualities\n";
	$qualities='ignore';
    }

    push @filepairs,[$infile,$outfile,$qualities];

}

unless (@filepairs) {
    print STDOUT "Everything seems to be mapped already\n";
    exit;
}

if ($usecluster) {
    # Build the submission file
    my $jobname='RNAseqMap';
    my $subfile=build_run_mapper_submission(\@filepairs,
					    $bindir,
					    $index,
					    $threads,
					    $jobname,
					    $mismatches,
					    $mappermem);
    my $queue=$options{'CLUSTER'};
    send2cluster($subfile,
		 $queue,
		 $jobname);

    # clean up
    my $command="rm $subfile";
    run_system_command($command);

    } else {
	my %mapper=%{get_mapper_routines()};
	foreach my $pair (@filepairs) {
	    # Map the file
	    $mapper{$mapper}->($index,
			       $pair->[0],
			       $pair->[1],
			       $threads,
			       $pair->[2],
			       $tempdir,
			       $mismatches,
			       $zip);
	}
}
exit;

sub build_run_mapper_submission {
    my $pairs=shift;
    my $bidir=shift;
    my $index=shift;
    my $threads=shift || 2;
    my $jobname=shift;
    my $mismatches=shift;
    my $mappermem=shift;

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

# Request 8 cpus and request memmory (if configured)

#\$ -l h_vmem=$mappermem
#\$ -pe smp $threads

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
ypwhich >&2
$bindir/run_mapper.RNAseq.parallel.pl -index $index -infile \$infile -outfile \$outfile -t $threads -mismatches $mismatches > \$infile.mapping.log
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile);
}

sub check_read_files {
    my $files=shift;
    my $paralleldir=shift;
    my $tempdir=shift;
    my $readdir=shift;

    my %locations;

    foreach my $file (keys %{$files}) {
	print STDERR "Checking $file...";
	my $filepath1=$paralleldir.'/'.$file;
	my $filepath2=$tempdir.'/'.$file;
	my $filepath3=$readdir.'/'.$file;
	if ( -r $filepath1) {
	    print STDERR "Present in $paralleldir\n";
	    $locations{$file}=$filepath1;
	} elsif (-r $filepath2) {
	    print STDERR "copying to $paralleldir\n";
	    $locations{$file}=$filepath1;
	    my $command="cp $filepath2 $filepath1";
	    run_system_command($command);
	} elsif (-r $filepath3) {
	    print STDERR "Present at $readdir\n";
	    $locations{$file}=$filepath3;
	} elsif (-r "$filepath3.gz") {
	    print STDERR "Unzipping in paralleldir\n";
	    my $command="gunzip -c $filepath3.gz > $filepath1";
	    run_system_command($command);
	    $locations{$file}=$filepath1;
	} else {
	    die "I can't find file $file\n";
	}
    }

    return(\%locations);
}
