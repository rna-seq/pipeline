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
# This script will take as input an annotation file.
# it will read the gene_RPKM table for the dataset in question and it will get
# a list of the genes that are expressed.
# for those genes that are not expressed it will fill in with zeros
# after it will read the annotation file and it will generate a table with all
# the necessary information for the gene from the annotation, for all those
# genes that are detected in at least one condition
### TO DO join the subroutines from build_top(genes,trans,exons) into one

use RNAseq_pipeline3 qw(get_fh parse_gff_line);
use RNAseq_pipeline_settings3 ('read_config_file','get_dataset_id',
			       'get_dbh');

# Declare some variables
my %genes;
my @files;
my %expressed;

my $annotation;
my $prefix;
my $gene_rpkm_table;
my $exonclasstab;
my $top_genes_table;
my $datasets_table;
my $db;
my $threshold=1;

my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$prefix=$options{'PREFIX'};
$gene_rpkm_table=$prefix.'_gene_RPKM';
$exonclasstab=$options{'EXONSCLASSTABLE'};
$top_genes_table=$prefix.'_top_genes_expressed';
$datasets_table=$prefix.'_dataset';
$db=$options{'DB'};

my $dbh=get_dbh();
my $dbhcommon=get_dbh(1);

# get some subs
*get_exon_number=get_exon_number_sub($dbhcommon,
				     $exonclasstab);
*get_trans_number=get_trans_number_sub($dbhcommon,
				       $exonclasstab);

my %distribution;

my %lanes=%{get_dataset_id($dbh,
			   $datasets_table)};

# Get the number of detected features in each lane
my @lanes=sort keys %lanes;

# Get the gene hits
for (my $i=0;$i<@lanes;$i++) {
    get_gene_expression(\%genes,
			$gene_rpkm_table,
			$dbh,
			$i,
			$lanes[$i],
			$threshold);
}

# Build the expression table
foreach my $gene (keys %genes) {
    for (my $i=0;$i<@lanes;$i++) {
	$expressed{$gene}->[$i]=$genes{$gene}->[$i] || 0;
    }
}

# Read the annotation and print out the table
print STDERR "Reading $annotation...\n";
my $annotfh=get_fh($annotation);
my %features;
my %geneinfo;
while (my $line=<$annotfh>) {
    if ($line=~/^#/o) {
	next;
    }
    my %line=%{parse_gff_line($line)};

    my $gene_id=$line{'feature'}{'gene_id'};

    if (defined $expressed{$gene_id}) {
	if ($geneinfo{$gene_id}{'chr'}) {
	    if ($geneinfo{$gene_id}{'start'} > $line{'start'}) {
		$geneinfo{$gene_id}{'start'}=$line{'start'};
	    }
	    
	    if ($geneinfo{$gene_id}{'end'} < $line{'end'}) {
		$geneinfo{$gene_id}{'end'}=$line{'end'};
	    }

	    if ($geneinfo{$gene_id}{'strand'} ne $line{'strand'}) {
		$geneinfo{$gene_id}{'strand'}='.';
	    }
	    if ($geneinfo{$gene_id}{'chr'} ne $line{'chr'}) {
		die "$features{$gene_id}{'chr'} ne $line{'chr'}\n";
	    }
	} else {
	    $geneinfo{$gene_id}{'start'}=$line{'start'};
	    $geneinfo{$gene_id}{'end'}=$line{'end'};
	    $geneinfo{$gene_id}{'strand'}=$line{'strand'};
	    $geneinfo{$gene_id}{'chr'}=$line{'chr'};
	}
    }
}
close($annotfh);

my $count=keys %geneinfo;
print STDERR $count,"\tFeatures expressed in total\n";

# get for each of the features the number of exons and the number of transcripts
# from the database
foreach my $gene (keys %geneinfo) {
#    print STDERR $gene,"\n";
    $features{$gene}->[0]=$geneinfo{$gene}{'end'} - $geneinfo{$gene}{'start'} + 1;
    $features{$gene}->[1]=$geneinfo{$gene}{'strand'};
    $features{$gene}->[2]=$geneinfo{$gene}{'chr'};
    $features{$gene}->[3]=get_exon_number($gene);
    $features{$gene}->[4]=get_trans_number($gene);
}

# print out the results
my $taboutfh=get_fh($top_genes_table.'.txt',1);
foreach my $gene (keys %expressed) {
    my $total=0;
    foreach my $value (@{$expressed{$gene}}) {
	$total+=$value;
    }
    print $taboutfh join("\t",
			 $gene,
			 @{$features{$gene}},
			 $total,
			 @{$expressed{$gene}}),"\n";
}
close($taboutfh);

# Build the table in the database
build_db_table($dbh,
	       $top_genes_table,
	       \@lanes,
	       \%lanes);

# Insert into the database
my $command="mysql $db < $top_genes_table.sql";
print STDERR "Executing: $command\n";
system($command);
$command="mysqlimport -L $db $top_genes_table.txt";
print STDERR "Executing: $command\n";
system($command);
$command="rm $top_genes_table.sql $top_genes_table.txt";
print STDERR "Executing: $command\n";
system($command);

exit;

sub build_db_table {
    my $dbh=shift;
    my $table=shift;
    my $lanes=shift;
    my $lane2id=shift;
    my $outtable=$table.'.sql';

    my ($query,$sth);

    $query ="DROP TABLE IF EXISTS $table;";
    $query.="CREATE TABLE $table( ";
    $query.='gene_id varchar(50) NOT NULL,';
    $query.='length int unsigned NOT NULL,';
    $query.='strand char(2) NOT NULL,';
    $query.='locus varchar(50) NOT NULL,';
    $query.='no_exons mediumint unsigned NOT NULL,';
    $query.='no_transcripts mediumint unsigned NOT NULL,';
    $query.='total mediumint unsigned NOT NULL,';
    foreach my $lane (@{$lanes}) {
	my $lane_id=$lane2id->{$lane};
	$lane_id=$lane;
	$query.="$lane mediumint unsigned NOT NULL,";
	$query.="INDEX idx_${lane} ($lane),";
    }
    $query.='INDEX idx_gene (gene_id)';
    $query.=');';

    my $outfh=get_fh($outtable,1);
    print $outfh $query,"\n";
    close($outtable);
}

sub get_exon_number_sub {
    my $dbh=shift;
    my $table=shift;

    my %cache;
    my ($query,$sth);

    $query ='SELECT count(DISTINCT exon_id) ';
    $query.="FROM $table ";
    $query.='WHERE gene_id = ?';
    $query.='GROUP BY gene_id';
    $sth=$dbh->prepare($query);

    my $exon_number=sub{
	my $gene_id=shift;

	unless ($cache{$gene_id}) {
	    my $count=$sth->execute($gene_id);
	    if ($count == 1) {
		($cache{$gene_id})=$sth->fetchrow_array();
	    } else {
		warn "Incorrect hits for $gene_id\n";
		$cache{$gene_id}=0;
	    }
	}

	return($cache{$gene_id});
    };

    return($exon_number);
}

sub get_trans_number_sub {
    my $dbh=shift;
    my $table=shift;

    my %cache;
    my ($query,$sth);

    $query ='SELECT count(DISTINCT transcript_id) ';
    $query.="FROM $table ";
    $query.='WHERE gene_id = ?';
    $query.='GROUP BY gene_id';
    $sth=$dbh->prepare($query);

    my $trans_number=sub{
	my $gene_id=shift;

	unless ($cache{$gene_id}) {
	    my $count=$sth->execute($gene_id);
	    if ($count == 1) {
		($cache{$gene_id})=$sth->fetchrow_array();
	    } else {
		warn "Incorrect hits for $gene_id\n";
		$cache{$gene_id}=0;
	    }
	}

	return($cache{$gene_id});
    };

    return($trans_number);
}

sub get_gene_expression {
    my $genes=shift;
    my $table=shift;
    my $dbh=shift;
    my $index=shift;
    my $lane=shift;
    my $threshold=shift || 1;

    print STDERR "Getting gene expression from $lane...\n";

    my ($query,$sth);
    $query ='SELECT gene_id, RPKM ';
    $query.="FROM $table ";
    $query.='WHERE LaneName = ?';
#    $query.=' limit 2';
    $sth=$dbh->prepare($query);
    $sth->execute($lane);

    while (my ($gene,$rpkm)=$sth->fetchrow_array()) {
	my $expression=int($rpkm + 0.5);
	if ($expression >= $threshold) {
	    $genes->{$gene}->[$index]=$expression
	}
    }

    my $count=keys %{$genes};
    print STDERR $count,"\tGenes expressed in $lane\n";
}


