#!/soft/bin/perl
# DGK

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
# This script should take a fasta or fastq file and trim it to the specified
# number of nucleotides

use Bio::SeqIO;
use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file');
use Bio::SeqIO::fastq;

my %options=%{read_config_file()};

my $infile;
my $qualities=$options{'QUALITIES'};
my $length=$options{'READLENGTH'};
my $format='fastq';

GetOptions('length|l=i' => \$length);

if ($qualities eq 'ignore') {
    $format='fasta';
}

my %parsers=%{get_parsing_subs($length,
			       $qualities)};
$infile=shift;
unless($infile) {
    die "No input file provided\n";
}

unless ($parsers{$format}) {
    die "Unknown format %format\n";
}
	     
print STDERR "Trimming $infile to $length nt reads\n";
  
# Parse the file
# Deal with zipped files
if ($infile=~/.gz$/) {
    $infile="zcat $infile |";
}
$parsers{$format}->($infile);

exit;

sub get_parsing_subs {
    my $read_length=shift;
    my $format=shift;
    my %parsing_subs;

    $parsing_subs{'fasta'}=sub {
	my $infn=shift;

	my $infh=Bio::SeqIO->new(-file => $infn,
				 -format => 'fasta');
	my $outfh=Bio::SeqIO->new(-fh => \*STDOUT,
				  -format => 'fasta');
	while (my $seq_obj=$infh->next_seq()) {
	    my $trimmed_seq=$seq_obj->trunc(1,$read_length);
	    $outfh->write_seq($trimmed_seq);
	}
	$infh->close();
    };

    $parsing_subs{'fastq'}=sub {
	my $infn=shift;

	if ($format eq 'phred') {
	    $format='fastq';
	} else {
	    $format='fastq-illumina';
	}

	my $in=Bio::SeqIO->new(-format    => $format,
				 -file      => $infn);

	my $out = Bio::SeqIO->new(-format    => $format,
				  -fh => \*STDOUT);

	while (my $seq = $in->next_seq) {
	    my $trimmed=$seq->trunc(1,$read_length);
	    $out->write_seq($trimmed);
	}

	return();
    };

    return(\%parsing_subs);
}
