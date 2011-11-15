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

# Converts a a file in fasta format to tbl format
# Usage: FastaToTbl.pl fasta_file

use Bio::SeqIO;
use RNAseq_pipeline_settings3 ('read_config_file');
use RNAseq_pipeline3 qw(get_fh);

my $transfile;

my %options=%{read_config_file()};
$transfile=$options{'TRANSCRIPTOMEFASTA'};

my $in  = Bio::SeqIO->new(-file => $transfile ,
			  -format => 'Fasta');

my $outfn=$options{'TRANSDIR'}.'/transcript.lengths';
my $outfh=get_fh($outfn,1);
while (my $seq = $in->next_seq() ) {
    print $outfh join("\t",
		      $seq->display_id(),
		      $seq->length()),"\n";
}
close($outfh);
$in->close();

exit;
