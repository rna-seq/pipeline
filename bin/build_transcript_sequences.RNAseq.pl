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

# Objective
# This script should take a gtf/gff file containing the annotation
# and a file containing the sequences for each exon, and it will build
# every annotated transcript for each gene

use RNAseq_pipeline3 ('get_fh','get_annotation_from_gtf','get_exon_sequences');
use RNAseq_pipeline_settings3 qw(read_config_file);
use Bio::Range;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqIO;

my $annotation_file;
my $transcriptsfile;
my $exonsfile;

my %options=%{read_config_file()};
$annotation_file=$options{'ANNOTATION'};
$transcriptsfile=$options{'TRANSCRIPTOMEFASTA'};
$exonsfile=$options{'EXONSFASTA'};

print STDERR "Extracting all transcript sequences from $annotation_file\n";

unless ($exonsfile &&
	$annotation_file) {
    die "ERROR:An exons file and an annotation file are required\n";
}

# Check if the file exists already
my $present=check_file($transcriptsfile);

if ($present) {
    print STDERR "$transcriptsfile is present. Skipping\n";
    exit;
}

print STDERR "Building all transcripts for each gene\n";
print STDERR "Processing $annotation_file\n";

# Get the sequence of the exons
my %exons=%{get_exon_sequences($exonsfile)};

# Extract all the trnscript coordinates belonging to each of the genes
my %genes=%{get_annotation_from_gtf($annotation_file)};


# get all the annotated transcripts
my %transcripts=%{build_transcriptome(\%genes,
				      \%exons)};

# Print out all transcripts
print STDERR "Printing transcript sequences in $transcriptsfile\n";
my $transfh=get_fh($transcriptsfile,1);
foreach my $transcript (keys %transcripts) {
    my $seq=$transcripts{$transcript};
    print $transfh ">$transcript\n";
    for (my $pos=0;$pos < length($seq) ; $pos += 60) {
	print $transfh substr($seq,$pos,60),"\n";
    }
}
close($transfh);
print STDERR "done\n";

exit;

sub build_transcriptome {
    my $genes=shift;
    my $exons=shift;

    my %transcripts;

    print STDERR "Extracting all annotated transcripts\n";
    foreach my $gene (keys %{$genes}) {
	# we cannot define the strand here using the gene because there may
	# (althought there should not) be cases where different annotated
	# transcripts belong to diffferent strands. If this is the case the gene
	# strand will be set to 0 and the if we use it to get the exon id this
	# id will belong to an unitialized element
	foreach my $trans ($genes->{$gene}->{'gene'}->transcripts()) {
	    my @exons=$trans->exons_ordered();
	    my @exon_seqs;
	    my $transcript_id=$trans->display_name();

	    for(my $i=0;$i < @exons; $i++) {
		my $exon=$exons[$i];
		my $exon_id=$exon->display_name();

		my $exon_seq=$exons->{$exon_id};
		push @exon_seqs,$exon_seq;
	    }
	    $transcripts{$transcript_id}=join('',
					      @exon_seqs);
	}
    }
    my $count=keys %transcripts;
    print STDERR $count,"\tTranscripts retrieved\n";
    return(\%transcripts);
}

sub check_file {
    my $file=shift;

    my $present=0;

    my $tablefh=get_fh($options{'PREFIX'}.'_transcript_seqs.txt',1);
    if (-r $file) {
	$present=1;
	# Print the table location for the junctions of this experiment
	print $tablefh join("\t",
			    $file,
			    "Present"),"\n";
    } else {
	# Continue, as the table must be created
	print STDERR $file,"\tIs not present\n";
	print $tablefh join("\t",
			    $file,
			    "Generated"),"\n";
    }
    close ($tablefh);

    return($present);
}
