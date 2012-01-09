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

# Author : David Gonzalez, david.gonzalez@crg.eu


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
# This script should read a bam file, extract the sam, sort it and produce an
# output file contaiing only the uniquely mapping reads in it

use Getopt::Long;
use Bio::DB::Sam;
use RNAseq_pipeline3 ('get_fh','run_system_command');
use RNAseq_pipeline_settings3 ('read_config_file');
#use Tools::Bam ('generate_sorted_bam');
use Data::Dumper;

# Declare variables & Get command line options
my $infile;
my $tmpdir;
my $add_cuff=1;

my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};

# Get the infile
$infile=shift;
unless($infile) {
    die "No input file provided\n";
}

# Process the file
print STDERR "Processing $infile\n";
print STDERR "Extracting uniquely mapping reads acordig to the alignmnent\n";

# produce a sam file sorted by read_id, pipe it into the script and print out
# only those cases where the read id appears only once.
my $outfile=filterbam($infile,
		      $tmpdir);

# Get the header from the original bam file
my $headerfile=get_header($infile);

# Concatenate both files and create the sorted bam
my $bamfile=build_output_bam($headerfile,
			     $outfile,
			     $infile);

# Clean up
my $command="rm $headerfile $outfile";
#run_system_command($command);

exit;

sub build_output_bam {
    my $header=shift;
    my $samfile=shift;
    my $bamfile=shift;

    $bamfile=~s/bam$/filtered/;
    $bamfile=~s/.*\///;

    my $command="cat $header $samfile | samtools view -S -b - |samtools sort - $bamfile";
    run_system_command($command);

    
}

sub get_header {
    my $infile=shift;

    my $outfile=$infile;

    $outfile=~s/.*\///;
    $outfile=$$.'.header.'.$outfile;

    print STDERR "Extracting header from $infile\n";

    my $command="samtools view -H $infile > $outfile";
    run_system_command($command);

    return($outfile);
}

sub filterbam {
    my $infile=shift;
    my $outfile=$infile;
    my $tmpdir=shift || '/tmp';

    $outfile=~s/.*\///;
    $outfile=$$.'.'.$outfile;

    my $command="samtools view $infile| sort -T $tmpdir -k1,1|";
    print STDERR $command,"\n";

    my $infh;
    open($infh,$command);

    my $outfh=get_fh($outfile,1);

    my $old_read_id='';
    my @lines=();

    while (my $line=<$infh>) {
	my @line=split("\t",$line);

	my $read_id=$line[0];
	my $flag=$line[1];
	my $cigar=$line[5];
#	my $flags=$line[11];
#	chomp($flags);

	# If the last previous read is differnt check if it is unique
	# If it is unique print it
	if ($read_id ne $old_read_id) {
	    if (@lines == 1) {		
		print $outfh @lines;
	    }
	    $old_read_id=$read_id;
	    @lines=();
	}
	push @lines,$line;
    }

    # Print the last record
    if (@lines == 1) {
	print $outfh @lines;
    }
    close($outfh);
    close($infh);

    return($outfile);
}
