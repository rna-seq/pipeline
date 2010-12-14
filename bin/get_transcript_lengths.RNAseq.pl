#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
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
