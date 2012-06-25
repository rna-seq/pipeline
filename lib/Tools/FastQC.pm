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
@EXPORT_OK=('run_fastqc');

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
use RNAseq_pipeline3 qw(run_system_command);
use Cwd;

sub parse_fastqc_report {
    my $fastqc_dir=shift;
    my %report;

    return(\%report)
}

sub run_fastqc {
    my $infile=shift;
    my $readdir=shift;
    my $outdir=shift;

    my $readfile=$readdir.'/'.$infile.'.gz';

    my $command='fastqc --nogroup ';
    $command.="-o $outdir ";
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
