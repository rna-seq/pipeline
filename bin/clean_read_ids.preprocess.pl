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
# This script should take a fasta or fastq file and clean the read ids removing
# any | or p from the end

use Bio::SeqIO;
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file');

my %options=%{read_config_file()};

my $infile;
my $qualities=$options{'QUALITIES'};
my $format='fastq';

if ($qualities eq 'ignore') {
    $format='fasta';
}

my %parsers=%{get_parsing_subs()};
$infile=shift;
unless($infile) {
    die "No input file provided\n";
}

unless ($parsers{$format}) {
    die "Unknown format %format\n";
}
	     
print STDERR "Cleaning $infile IDs\n";
  
# Parse the file
# Deal with zipped files
if ($infile=~/.gz$/) {
    $infile="zcat $infile |";
}
$parsers{$format}->($infile);

exit;

sub get_parsing_subs {
    my %parsing_subs;

    $parsing_subs{'fasta'}=sub {
	my $infn=shift;

	my $infh=Bio::SeqIO->new(-file => $infn,
				 -format => 'fasta');
	my $outfh=Bio::SeqIO->new(-fh => \*STDOUT,
				  -format => 'fasta');
	while (my $seq_obj=$infh->next_seq()) {
	    my $id=$seq_obj->display_id();
	    $id =~s/\|p1$/\/1/o;
	    $id =~s/\|p2$/\/2/o;
	    $seq_obj->display_id($id);
	    $outfh->write_seq($seq_obj);
	}
	$infh->close();
    };
    $parsing_subs{'fastq'}=sub {
	my $infn=shift;

	my $infh=get_fh($infn);

	my $line_type=0;
	my %sequence;
	while (my $line=<$infh>) {
	    chomp($line);

	    # Decide which line we are in
	    if ($line=~s/^@//) {
		my $id=$line;
		$id =~s/\|p1$/\/1/o;
		$id =~s/\|p2$/\/2/o;
		%sequence=('id' => $id);
		$line_type=1;
		next;
	    } elsif ($line=~/^\+/) {
		$line_type=2;
		next;
	    }
	    
	    if ($line_type == 1) {
		$sequence{'seq'}.=$line;
		next;
	    } elsif ($line_type == 2) {
		$sequence{'qual'}.=$line;
	    } else {
		warn "Shouldn't be here\n";
	    }
	    
	    my $seq=$sequence{'seq'};
	    my $qual=$sequence{'qual'};

	    print '@',$sequence{'id'},"\n";
	    print "$seq\n";
	    print "+\n";
	    print "$qual\n";
	}
	close($infh);

	return();
    };

    return(\%parsing_subs);
}
