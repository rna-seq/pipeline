package GRAPE::Formats::Fasta;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export subroutines to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
push @EXPORT_OK,('check_fasta_file');

use strict;
use warnings;

# Bio::Perl modules
use Bio::SeqIO;

# Check if the fasta file identifiers contain any spaces and warn if this is
# the case
sub check_fasta_file {
    my $file=shift;
    my $fileok=1;

    print STDERR "Checking the genome file...";
    my $in  = Bio::SeqIO->new(-file => $file ,
			      -format => 'Fasta');

    while ( my $seq = $in->next_seq() ) {
	my $seqid=$seq->display_id();
	my $desc=$seq->desc();
	if ($desc=~/[=~:;]/o) {
	    print STDERR "Presence of special characters in the header may cause parsing problems after mapping\n";
	    print STDERR "Check: $desc\n";
	    $fileok=0;
	} elsif ($seqid=~/^\s*\w+\s*$/o) {
	    next;
	} else {
	    print STDERR "WARNING: Spacing in the fasta identifiers may cause problems: $seqid\n";
	}
    }
    $in->close();
    if ($fileok) {
	print STDERR "Fasta file seems fine\n";
    }
    return($fileok);
}

1;
