package Tools::FastQC;

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
@EXPORT_OK=('run_fastqc','parse_fastqc_report');

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
use RNAseq_pipeline3 qw(run_system_command get_fh);
use Cwd;

sub parse_fastqc_report {
    my $fastqc_dir=shift;
    my %report;

    # Check that the report file exists
    my $reportfile=$fastqc_dir.'/fastqc_data.txt';
    unless (-s $reportfile) {
	die "I cant find the report file $reportfile\n";
    }

    my $reportfh=get_fh($reportfile);
    while (my $line=<$reportfh>) {
	chomp($line);
	if ($line=~/^\>\>/) {
	    # Begin block
	    $line=~s/^>>//;
	    my ($block,$status)=split("\t",$line);
	    my $key=$block;
	    $key=~s/\s//g;
	    $report{$key}{'status'}=[$block,$status];

	    # Get the header and values
	    while (my $line2=<$reportfh>) {
		chomp($line2);
		if ($line2=~/^>>END_MODULE/) {
		    # end block
		    last;
		} elsif ($line2=~/^#/) {
		    $line2=~s/^#//;
		    my @header=split("\t",$line2);
		    push @{$report{$key}{'header'}},[@header];
		} else {
		    my @values=split("\t",$line2);
		    push @{$report{$key}{'values'}},[@values];
		}
	    }

	} else {
	    next;
	}
    }
    close($reportfh);

    # print the results
#    foreach my $key (keys %report) {
#	foreach my $value (@{$report{$key}{'header'}}) {
#	    print STDERR join("\t",
#			  @{$value}),"\n";
#	}
#	foreach my $value (@{$report{$key}{'values'}}) {
#	    print STDERR join("\t",
#			  @{$value}),"\n";
#	}
#    }

    return(\%report)
}

sub run_fastqc {
    my $infile=shift;
    my $readdir=shift;
    my $outdir=shift;

    my $readfile=$readdir.'/'.$infile;
    unless ($infile=~/.bam$/) {
	$readfile.='.gz';
    }

    my $command='fastqc --nogroup ';
    $command.="-o $outdir ";
    $command.="-t 1 ";
    $command.="$readfile";

    # Build the path for the fastqc directory
    my $fastqcdir=$outdir.'/'.$infile;
    $fastqcdir=~s/\.[^\.]+?$//;
    $fastqcdir.='_fastqc';

    # Check if the directory is already present and if it is skip the creation
    if (-d $fastqcdir) {
        print STDERR $fastqcdir,"\tIs already present skipping fastqc\n";
    } else {
       run_system_command($command);
       # First check if the directory was created correctly
       if (-d $fastqcdir) {
  	   print STDERR $fastqcdir,"\tcreated correctly\n";
       } else {
	   die $fastqcdir,"\tdoes not seem to exist\n";
       }
    }
    return($fastqcdir);
}

1;
