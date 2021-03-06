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
# This script will take the necessary commands for running the recursive mapping
# on the cluster

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_dbh send2cluster);
use RNAseq_GEM3 ('check_index','determine_quality_type','get_mapper_routines',
		 'check_input');

my $index;
my $mapper;
my $outdir;
my $filetype;
my $mismatches;
my $paralleldir;
my $readdata;
my $projdir;
my $bindir;
my $splittable;
my $splitdir;
my $recmapdir;
my $queue='main';
my $threads=2;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'outdir|o=s' => \$outdir);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$mismatches=$options{'MISMATCHES'};
$paralleldir=$options{'LOCALPARALLEL'};
$readdata=$options{'READDIR'};
$projdir=$options{'PROJECT'};
$bindir=$options{'BIN'};
$splittable=$options{'PREFIX'}.'_split_mapping';
$splitdir=$projdir.'/'.$options{'SPLITMAPDIR'};
$recmapdir=$projdir.'/'.$options{'RECMAPDIR'};
$queue=$options{'CLUSTER'};
$threads=$options{'THREADS'};

# Decide where to put the output
unless ($outdir) {
    die "Where should I put the results ???\n";
}
$outdir=~s/\/$//;
my $usecluster=0;

if ($options{'CLUSTER'}) {
    print STDERR "Mapping in the cluster\n";
    $usecluster=1;
} else {
    print STDERR "I can only use a cluster if I know its name. Mapping locally\n";
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
my %files=%{split_file_list($splittable,
			    $splitdir)};

# Check each of the files and determine if it is present in the paralleldir
# directory if not unzip it and copy it there.
my %locations=%{check_read_files(\%files,
				 $paralleldir,
				 $splitdir)};

my @filepairs;
foreach my $file (keys %locations) {
    my $qualities;
    # Check the input file
    my $infile=$locations{$file};

    # Make the outfile
    my $outdir=$recmapdir;

    # Set the name of the final output, and exit if it is already present
    my $basename=$infile;
    $basename=~s/.*\///;
    my $final_mapping=$outdir.'/'.$basename.".recursive.map";
    if (-r $final_mapping) {
	print STDERR "$final_mapping is already present. Skipping...\n";
    } else {
	push @filepairs,[$infile,$outdir];
    }
}

unless (@filepairs) {
    print STDOUT "Everything seems to be mapped already\n";
    exit;
}

if ($usecluster) {
    # Build the submission file
    my $jobname='RNAseqRecMap';
    my $subfile=build_run_mapper_submission(\@filepairs,
					    $bindir,
					    $index,
					    $jobname,
					    $threads);
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
	    # Run the recursive mapper for each of the files
	    my $command="$bindir/run_recursive_mapper.RNAseq.pl ";
	    my $genomeindex=$options{'GENOMEINDEX'};
	    my $recmapdir=$options{'RECMAPDIR'};
	    $command.="-index $genomeindex ";
	    $command.='-i '.$pair->[0].' ';
	    $command.="-o $recmapdir > ";
	    $command.=$pair->[1].'/rec.mapping.log';
	    run_system_command($command);
	}
}
exit;

# This subroutine will get a list of the split-mapped files that should be
# present
sub split_file_list {
    my $table=shift;
    my $splitdir=shift;

    my %files;
    my $dbh=get_dbh();

    my ($query,$sth);
    $query ='SELECT filename ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($file)=$sth->fetchrow_array()) {
	$files{$file}=1;
    }

    return(\%files);
}

sub build_run_mapper_submission {
    my $pairs=shift;
    my $bidir=shift;
    my $index=shift;
    my $jobname=shift;
    my $threads=shift;

    # Threads are hardcoded to 2 here in order to prevent wasting too much CPU
    # time during the split-mapping which is not parallelized.

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
#\$ -pe smp $threads

#\$ -l h_vmem=6G

# Write in to the current working directory
#\$ -cwd 
export PATH=\$PATH:/soft/bin
# Make sure the sorting order is "a la C"
export LC_ALL=C

infiles=(@infiles)
outfiles=(@outfiles)

export infile=\${infiles[\$SGE_TASK_ID-1]}
export outdir=\${outfiles[\$SGE_TASK_ID-1]}

echo \$HOSTNAME >&2
$bindir/run_recursive_mapper.RNAseq.pl -index $index -infile \$infile -outdir \$outdir > \$infile.mapping.log
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile);
}

sub check_read_files {
    my $files=shift;
    my $paralleldir=shift;
    my $readdir=shift;

    my %locations;

    foreach my $file (keys %{$files}) {
	print STDERR "Checking $file...";
	my $filepath1=$paralleldir.'/'.$file;
	my $filepath3=$readdir.'/'.$file;
	if ( -r $filepath1) {
	    print STDERR "Present at $readdir\n";
	    $locations{$file}=$filepath1;
	}  elsif (-r $filepath3) {
	    print STDERR "Present at $readdir\n";
	    $locations{$file}=$filepath3;
	} elsif (-r "$filepath3.gz") {
	    print STDERR "Unzipping in $paralleldir\n";
	    my $command="gunzip -c $filepath3.gz > $filepath1";
	    run_system_command($command);
	    $locations{$file}=$filepath1;
	} else {
	    die "I can't find file $file\n";
	}
    }

    return(\%locations);
}
