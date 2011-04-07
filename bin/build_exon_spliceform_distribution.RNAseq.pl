#!/soft/bin/perl

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
# This script should take a gtf/gff or genpred file containing exons with
# gene names and it will build a set of all annotated exon-exon
# junctions for each of the genes. It will also output the sequence for all
# the annotated exons
# The introns, regardless of the strand will alway be ordered from the 5'
# coordinate to the 3'

use RNAseq_pipeline3 ('get_fh','parse_gff_line',
		      'get_annotation_from_gtf');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','check_db',
			       'get_feature_overlap_sub');
use Bio::Range;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqIO;

# Declare variables
my $annotation;
my $exonsfile;
my $reported;
my $overlap;
my $stranded;
my $exondir;
my $exontable;
my $junctable;
my $database;
my $command;
my $prefix;
my $parallel='cluster';
my $paralleltmp;
my $bindir;

# Read the options file
my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$exondir=$options{'EXONDIR'};
$prefix=$options{'PREFIX'};
$exonsfile=$exondir.'/'.$options{'PREFIX'}.'.exon.gtf';
$overlap=$options{'PREFIX'}.'.exon';
$stranded=$options{'STRANDED'};
$exontable=$options{'EXONSCLASSTABLE'};
$junctable=$options{'JUNCTIONSCLASSTABLE'};
$database=$options{'COMMONDB'};
$paralleltmp=$options{'PROJECT'}.'/'.$options{'LOCALPARALLEL'};
$bindir=$options{'BIN'};
unless($options{'CLUSTER'}) {
    print STDERR "Running locally\n";
    $parallel='default';
}

unless ($paralleltmp) {
    $paralleltmp=$options{'PROJECT'}.'/work';
}

# Get some subroutines
*get_feature_overlap=get_feature_overlap_sub($parallel,
					     $paralleltmp,
					     $bindir);

# First check if the tables need to be created or not
my $exon_present=check_db($exontable,
			  $prefix.'_exon_classification.txt');
my $junction_present=check_db($junctable,
			      $prefix.'_junction_classification.txt');

if ($exon_present) {
    print STDERR "$exontable is present. Skipping\n";
} else {
    # Prepare for table generation
    make_exon_table($exontable)
}

if ($junction_present) {
    print STDERR "$junctable is present. Skipping\n";
} else {
    # Prepare for table generation
    make_junctions_table($junctable)
}

if ($exon_present && $junction_present) {
    exit;
}

# Run overlap on the annotation files
print STDERR "Finding overlapping exons\n";
get_feature_overlap($exonsfile,
		    $exonsfile,
		    $stranded,
		    $overlap);

print STDERR "Building junctions\n";

print STDERR "Processing $annotation\n";

# get the overlap between all the exons
my %overlap;
my $overlapfile=$overlap.'.overlap.gz';
%overlap=%{get_overlap($overlapfile)};

# Extract all the exons coordinates belonging to each of the genes
my %exons=%{get_annotation_from_gtf($annotation)};

# get all the annotated junctions
my ($annotated_junctions,$annotated_exons)=build_annotated_junctions2(\%exons);

# Print out all exon sequences
unless($exon_present) {
    print STDERR "Generating $exontable\n";
    my $exonfh=get_fh("$exontable.txt",1);
    foreach my $gene (keys %{$annotated_exons}) {
	foreach my $trans (keys %{$annotated_exons->{$gene}}) {
	    foreach my $exon (keys %{$annotated_exons->{$gene}->{$trans}}) {
		my ($start,$end)=(split('_',$exon))[-3,-2];
		my $overlapping=$overlap{$exon};
		unless (defined($overlapping)) {
		    print STDERR join("\t",
				      $exon,
				      $start,
				      $end,
				      $gene,
				      $trans),"\n";
		}
		print $exonfh join("\t",
				   $exon,
				   $start,
				   $end,
				   $gene,
				   $trans,
				   $overlapping),"\n";
	    }
	}
    }
    close($exonfh);

    # Add to the database
    $command="mysql $database < $exontable.sql";
    print STDERR "Executing:\t",$command,"\n";
    system($command);
    $command="mysqlimport -L $database $exontable.txt";
    print STDERR "Executing:\t",$command,"\n";
    system($command);
    $command="rm $exontable.sql $exontable.txt";
    print STDERR "Executing:\t",$command,"\n";
    system($command);
}

# Print out all exon junction sequences
unless($junction_present) {
    print STDERR "Generating $junctable\n";
    my $junctionfh=get_fh("$junctable.txt",1);
    foreach my $gene (keys %{$annotated_junctions}) {
	foreach my $trans (keys %{$annotated_junctions->{$gene}}) {
	    foreach my $junction (keys %{$annotated_junctions->{$gene}->{$trans}}) {
		print $junctionfh join("\t",
				       $junction,
				       $gene,
				       $trans),"\n";
	    }
	}
    }
    close($junctionfh);
    # Add to the database
    $command="mysql $database < $junctable.sql";
    print STDERR "Executing:\t",$command,"\n";
    system($command);
    $command="mysqlimport -L $database $junctable.txt";
    print STDERR "Executing:\t",$command,"\n";
    system($command);
    $command="rm $junctable.sql $junctable.txt";
    print STDERR "Executing:\t",$command,"\n";
    system($command);
}

exit;

sub make_exon_table {
    my $table=shift;
    my $tablefh=get_fh("$table.sql",1);
    print $tablefh "DROP TABLE IF EXISTS $table;
CREATE TABLE $table (
       exon_id varchar(50) NOT NULL,
       start int unsigned not null,
       end int unsigned not null,
       gene_id varchar(50) NOT NULL,
       transcript_id varchar(50) NOT NULL,
       overlaps mediumint unsigned NOT NULL,
       index idx_exon (exon_id),
       index idx_start (start),
       index idx_end (end),
       index idx_gene (gene_id),
       index idx_transcript (transcript_id)
);\n";
    close($tablefh);
}

sub make_junctions_table {
    my $table=shift;
    my $tablefh=get_fh("$table.sql",1);
    print $tablefh "DROP TABLE IF EXISTS $table;
CREATE TABLE $table (
       junction_id varchar(100) NOT NULL,
       gene_id varchar(50) NOT NULL,
       transcript_id varchar(50) NOT NULL,
       index idx_junc (junction_id),
       index idx_gene (gene_id),
       index idx_transcript (transcript_id)
);\n";
    close($tablefh);
}

sub get_overlap {
    my $filename=shift;
    my $fh=get_fh($filename);
    my %overlaps;

    print STDERR "Extracting the overlap information from $filename...";
    while (my $line=<$fh>) {
	my %line=%{parse_gff_line($line)};
	my $strand;

	if ($line{'strand'} eq '-') {
	    $strand= -1;
	} elsif ($line{'strand'} eq '+') {
	    $strand=1;
	} else {
	    warn "unknown strand $line{'strand'}\n";
	}

	# Skip the same stuff we are skipping for the parsing of the annotation
	my $chr=$line{'chr'};
#	if ($chr=~/random/io) {
#	    next;
#	}
	if ($chr=~/hap/o) {
	    next;
	} elsif ($chr=~/^chrU/o) {
	    next;
	} elsif ($chr=~/^Un\./o) {
	    # This is for EnsEMBL cow
	    next;
	} elsif ($chr=~/^(chr)?HSCHR/o) {
	    next;
	} elsif ($chr=~/^AAFC03011182/o) {
	    # EnsEMBL cow
	    next;
	}



#	if ($chr=~/chrMT/) {
#	    $chr=~s/chrMT/chrM/;
#	}
	my $exon_id=join('_',
			 $chr,
			 $line{'start'},
			 $line{'end'},
			 $strand);
	my $overlaps=$line{'feature'}{'nb_ov_feat2:'};
	if ($overlaps) {
	    $overlaps{$exon_id}=$overlaps - 1;
	} else {
	    warn "No overlaps for $exon_id\n";
	    print STDERR $line;
	}
    }
    print STDERR "done\n";

    return(\%overlaps);
}

sub build_annotated_junctions2 {
    my $exons=shift;
    my %annotated_junctions;
    my %annotated_exons;

    print STDERR "Extracting all annotated junctions\n";
    foreach my $gene (keys %{$exons}) {
	foreach my $trans ($exons->{$gene}->{'gene'}->transcripts()) {
	    my $strand=$trans->strand();
	    my $trans_id=$trans->display_name();
	    my @exons=$trans->exons_ordered();
	    for(my $i=0;$i<@exons - 1; $i++) {
		my $exon1=$exons[$i];
		my $exon2=$exons[$i + 1];
		my $exon1_id=$exon1->display_name();

		my $exon2_id=$exon2->display_name();

		my $junction_id=$exon1_id.'_splice_'.$exon2_id;
		if ($strand == -1) {
		    $junction_id=$exon2_id.'_splice_'.$exon1_id;
		}
		$annotated_junctions{$gene}{$trans_id}{$junction_id}=1;
	    }
	    for(my $i=0;$i<@exons; $i++) {
		my $exon=$exons[$i];
		my $exon_id=$exon->display_name();
		$annotated_exons{$gene}{$trans_id}{$exon_id}=1;
	    }
	}
    }
    my $count=keys %annotated_junctions;
    print STDERR $count,"\tGenes  contain annotated junctions\n";
    return(\%annotated_junctions,\%annotated_exons);
}



