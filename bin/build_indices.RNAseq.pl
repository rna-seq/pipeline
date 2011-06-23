#!/soft/bin/perl
# DGK

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
# This script should take a list of files and check if there are indices
# available for them. If not it will create them.

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file);
use RNAseq_GEM3 ('check_index');
use Bio::SeqIO;

my $exonsfile;
my $junctionsfile;
my $genomefile;
my $transcriptomefile;
my $database;
my $genomeindex;
my $junctionsindex;
my $transcriptomeindex;
my $exonsindex;

my %options=%{read_config_file()};
$exonsfile=$options{'EXONSFASTA'};
$junctionsfile=$options{'JUNCTIONSFASTA'};
$genomefile=$options{'GENOMESEQ'};
$transcriptomefile=$options{'TRANSCRIPTOMEFASTA'};

$database=$options{'DB'};

$genomeindex=$options{'GENOMEINDEX'};
$junctionsindex=$options{'JUNCTIONSINDEX'};
$transcriptomeindex=$options{'TRANSCRIPTOMEINDEX'};

unless ($exonsfile && $junctionsfile && $genomefile && $transcriptomefile) {
    die "ERROR:Files for exons, transcript, junctions and genome are needed\n";
}

unless ($junctionsindex && $genomeindex && $transcriptomeindex) {
    die "ERROR:Location of transcript, junctions and genome indices is needed\n";
}

print STDERR "Determining which indices have to be created\n";

# Check if the indices exist already
my %present;
foreach my $index ($junctionsindex,$genomeindex,$transcriptomeindex) {
    $present{$index}=check_index($index);
}

# For each of the missing indices create them:
my %indices=($genomeindex => $genomefile,
	     $transcriptomeindex => $transcriptomefile,
	     $junctionsindex => $junctionsfile);

my %type=($genomeindex => 'genome',
	  $transcriptomeindex => 'transcriptome',
	  $junctionsindex => 'junctions');

my $out_table=$options{'PREFIX'}.'_indices.txt';
my $outfh=get_fh($out_table,1);
foreach my $index (keys %present) {
    if ($present{$index}) {
	print $outfh join("\t",
			  $index,
			  $type{$index},
			  'Present'),"\n";
    } else {
	build_index($index,
		    $indices{$index});
	print $outfh join("\t",
			  $index,
			  $type{$index},
			  'Generated'),"\n";
    }
}
close($outfh);

exit;

sub build_index {
    my $index=shift;
    my $file=shift;

    $index=~s/(.*\/)//;
    my $indexdir=$1;

    # Check if the headers in the fasta file are OK if not copy a "good" version
    # of the file into the work directory and build the index with that.


    my $command="cd $indexdir; gem-do-index -i $file -o $index;";
    print STDERR "Executing $command\n";
    system($command);
}
