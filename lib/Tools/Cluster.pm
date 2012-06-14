package Tools::Cluster;

#  GRAPE
#  Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('send2cluster','build_script_submission');

use strict;
use warnings;

# Load other modules required
use POSIX qw(uname);

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0) || die "Unable to determine absolute path: $!\n";
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# This file contains subroutines used for submitting jobs to the cluster

# Load subs from other modules
use RNAseq_pipeline3 ('MySQL_DB_Connect','get_fh','run_system_command',
		      'parse_gff_line','check_field_existence');
use Cwd;

### Sub to send stuff to the cluster
sub send2cluster {
    my $file=shift;
    my $queue=shift;
    my $jobname=shift;
    my %options=%{read_config_file()};

    if ($queue && 
	($queue ne '-')) {
	print STDERR "Using $queue\n";
    } elsif ($options{'CLUSTER'}  &&
	     ($options{'CLUSTER'} ne '-')) {
	$queue=$options{'CLUSTER'};
	warn "CLUSTER had no value, there may be a problem set to: $queue\n";
    } else {
	die "No valid queue is defined, so I can't submit to the cluster\n";
    }

    my $logs=$options{'LOGS'};

    # Submit to the cluster
    my $command="qsub -q $queue $file";
    sleep(1);
    my $execute=`$command`;

    # monitor the execution
    my ($job_id)=(split(/\s+/,$execute))[2];
    if ($job_id) {
	print STDERR "Running job $job_id on $queue\n";
    } else {
	warn "WARNING: No job ID returned by cluster. Something when wrong\n";
	die "Better luck next time...\n";
    }

    $job_id=~s/\..+$//;    
    $command="qstat -j $job_id 2>&1";
    print STDERR "Waiting for job $job_id";
    while (1) {
	my @running=`$command`;
	my $success=$?;
	my $finished=1;

	while (my $line=shift(@running)) {
	    if ($line=~/^=+$/) {
		next;
	    } elsif ($line=~/^job_number/o) {
		# The job is queued or running
		my @line=split(/\s+/,$line);
		$finished=0;
		print STDERR '.';
		last;
	    } elsif ($line=~/^Following\sjobs\sdo\snot\sexist/o) {
		# The job has finished
		last;
	    } else {
		# There is a problem
		print STDERR $line,"\n";
		print STDERR "\nProblem\n",@running,"\n";
		die;
	    }
	}
	
	if ($finished) {
	    print STDERR "done\n";
	    last;
	} else {
	    sleep(60);
	}
    }

    # collect the output of the cluster node into the log file for the
    # step
    # First indicate which job the output beloiongs to:
    $command="echo $job_id >> $logs/$jobname.$queue.log";
    run_system_command($command);

    # Add the command output
    $command="cat $jobname.*[eo]$job_id* >> $logs/$jobname.$queue.log";
    run_system_command($command);

    # Remove the files
    $command="rm $jobname.*[eo]$job_id*";
    run_system_command($command);

    return($job_id);
}


sub build_script_submission {
    my $script=shift;
    my $bindir=shift;
    my $memmory=shift;
    my $threads=shift || 1;

    unless ($bindir) {
	die "I don't know where to find the binaries\n";
    }
     
    unless($script) {
	die "No input supplied\n";
    }

    my $jobname=$script;
    $jobname=~s/(.pl)?$//;
    $jobname=~s/_//g;
    
    # Print the submission file
    my $subfile="subfile.$$.job";
    my $outfh=get_fh($subfile,1);
    
    print STDERR "Building submission file for $jobname...";

    print $outfh <<FORMEND;
# Get the job name
#\$ -N $jobname
    
#\$ -l h_vmem=$memmory
#\$ -pe smp $threads

# Write in to the current working directory
#\$ -cwd

export PATH=\$PATH:/soft/bin
export LC_ALL='C'

$bindir/$script
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile,$jobname);
}

# For running overlap in the cluster
sub build_run_overlap_submission {
    my $pairs=shift;
    my $bindir=shift;
    my $flags=shift;
    my $annotation=shift;
    my $jobname=shift;

    unless ($jobname) {
	die "No jobname provided\n";
    }

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
infiles=(@infiles)
outfiles=(@outfiles)

export infile=\${infiles[\$SGE_TASK_ID-1]}
export outfile=\${outfiles[\$SGE_TASK_ID-1]}

echo \$HOSTNAME >&2
until ( [[ -s \$outfile ]] );
do $bindir/overlap $annotation \$infile $flags > \$outfile; sleep 1;
if [[ -s \$outfile ]]; then echo done; else echo repeating; fi
done

FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile);
}



1;
