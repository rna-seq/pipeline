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
# This script will take as an input the overlap result of an experiment and it
# will extract from it the genes that have at least one RNA seq hit that
# overlaps completely with the gene 

# Load some modules
use RNAseq_pipeline3 qw(get_fh parse_gff_line);

# Declare some variables
my $threshold=1;

my $infile=shift;
my %genes;

# Check that we have the correct input files
unless($infile) {
    die "Input sequence file required\n";
}

# First read in the file
my $infh=get_fh($infile);
while (my $line=<$infh>) {
    my %line=%{parse_gff_line($line)};
    my $gene_id=$line{'feature'}{'gene_id'};
    my $hits=$line{'feature'}{'list_feat2:'};
    unless ($gene_id) {
	warn "No gene ID found in $line\n";	
    }
    $gene_id=~s/("|;)//g;
    # Check if the hit is within the target
    my @targets=split(',',$hits);
    my $start=$line{'start'};
    my $end=$line{'end'};
    my $length=$end - $start + 1;
    $genes{$gene_id}->[2]=$line{'type'};
    $genes{$gene_id}->[1]=$length;
    my @valid;

    if ($line{'feature'}{'nb_ov_feat2:'} > 0) {
	foreach my $hit (@targets) {
	    # This is because sometimes the scaffold has an undescore
	    my @parts=split('_',$hit);
	    my $h_strand=pop(@parts);
	    my $h_end=pop(@parts);
	    my $h_start=pop(@parts);
	    my $h_chr=join('_',@parts);
	    unless ($h_chr eq $line{'chr'}) {
		warn "Wrong chromosome\n";
	    }
	    if (($h_start >= $start) &&
		($h_end <= $end)) {
		push @valid, $hit;
	    }
	}

	if (@valid) {
	    $line{'feature'}{'nb_ov_feat2'}=@valid;
	    $line{'feature'}{'list_feat2:'}=join(',',@valid);
	    $genes{$gene_id}->[0]=@valid;
	} else {
	    $genes{$gene_id}->[0]=0;
	}
    } else {
	$genes{$gene_id}->[0]=0
    }
}
close($infh);
# print results
foreach my $gene (keys %genes) {
    unless ($genes{$gene}->[0]) {
	next;
    }
    my $coverage=sprintf "%.2f",$genes{$gene}->[0]/$genes{$gene}->[1];
    print join("\t",
	       $gene,
	       @{$genes{$gene}},
	       $coverage),"\n";
}

exit;
