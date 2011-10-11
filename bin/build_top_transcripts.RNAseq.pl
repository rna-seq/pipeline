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
# it will read the transcript_RPKM table for the dataset in question and it will
# get a list of the transcripts that are expressed.
# for those that are not expressed it will fill in with zeros
# after it will read the annotation file and it will generate a table with all
# the necessary information for the transcript from the annotation, for all 
# transcripts detected in at least one condition

use RNAseq_pipeline3 qw(get_fh parse_gff_line check_table_existence);
use RNAseq_pipeline_settings3 ('read_config_file','get_dataset_id',
			       'get_dbh');

# Declare some variables
my %trans;
my @files;
my %expressed;

my $annotation;
my $prefix;
my $trans_rpkm_table;
my $exonclasstab;
my $top_trans_table;
my $datasets_table;
my $db;
my $threshold=1;

my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$prefix=$options{'PREFIX'};
$trans_rpkm_table=$prefix.'_transcript_expression_levels';
$exonclasstab=$options{'EXONSCLASSTABLE'};
$top_trans_table=$prefix.'_top_transcripts_expressed';
$datasets_table=$prefix.'_dataset';
$db=$options{'DB'};

my $dbh=get_dbh();
my $commondbh=get_dbh(1);

# get some subs
*get_exon_number=get_exon_number_sub($commondbh,
				     $exonclasstab);

my %distribution;

my %lanes=%{get_dataset_id($dbh,
			   $datasets_table)};

# Get the number of detected features in each lane
my @lanes=sort keys %lanes;

# Get the gene hits
for (my $i=0;$i<@lanes;$i++) {
    get_trans_expression(\%trans,
			$trans_rpkm_table,
			$dbh,
			$i,
			$lanes[$i],
			$threshold);
}

# Build the expression table
foreach my $trans (keys %trans) {
    for (my $i=0;$i<@lanes;$i++) {
	$expressed{$trans}->[$i]=$trans{$trans}->[$i] || 0;
    }
}

# Read the annotation and print out the table
my $annotfh=get_fh($annotation);
my %features;
my %transinfo;
while (my $line=<$annotfh>) {
    if ($line=~/^#/) {
	next;
    }
    my %line=%{parse_gff_line($line)};
    my $trans_id=$line{'feature'}{'transcript_id'};

    if ($trans_id && defined $expressed{$trans_id}) {
	if ($transinfo{$trans_id}{'chr'}) {
	    if ($transinfo{$trans_id}{'start'} > $line{'start'}) {
		$transinfo{$trans_id}{'start'}=$line{'start'};
	    }
	    
	    if ($transinfo{$trans_id}{'end'} < $line{'end'}) {
		$transinfo{$trans_id}{'end'}=$line{'end'};
	    }

	    if ($transinfo{$trans_id}{'strand'} ne $line{'strand'}) {
		$transinfo{$trans_id}{'strand'}='.';
	    }
	    if ($transinfo{$trans_id}{'chr'} ne $line{'chr'}) {
		die "$transinfo{$trans_id}{'chr'} ne $line{'chr'}\n";
	    }
	} else {
	    $transinfo{$trans_id}{'start'}=$line{'start'};
	    $transinfo{$trans_id}{'end'}=$line{'end'};
	    $transinfo{$trans_id}{'strand'}=$line{'strand'};
	    $transinfo{$trans_id}{'chr'}=$line{'chr'};
	}
    }
}
close($annotfh);

# get for each of the features the number of exons and the number of transcripts
# from the database
foreach my $trans (keys %transinfo) {
#    print STDERR $gene,"\n";
    $features{$trans}->[0]=$transinfo{$trans}{'end'} - $transinfo{$trans}{'start'} + 1;
    $features{$trans}->[1]=$transinfo{$trans}{'strand'};
    $features{$trans}->[2]=$transinfo{$trans}{'chr'};
    $features{$trans}->[3]=get_exon_number($trans);
}

# print out the results
my $taboutfh=get_fh($top_trans_table.'.txt',1);
foreach my $trans (keys %expressed) {
    my $total=0;
    foreach my $value (@{$expressed{$trans}}) {
	$total+=$value;
    }
    print $taboutfh join("\t",
			 $trans,
			 @{$features{$trans}},
			 $total,
			 @{$expressed{$trans}},
			 ),"\n";
}
close($taboutfh);

# Build the table in the database
build_db_table($dbh,
	       $top_trans_table,
	       \@lanes,
	       \%lanes);

# Insert into the database
my $command="mysql $db < $top_trans_table.sql";
print STDERR "Executing: $command\n";
system($command);
$command="mysqlimport -L $db $top_trans_table.txt";
print STDERR "Executing: $command\n";
system($command);
$command="rm $top_trans_table.sql $top_trans_table.txt";
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
    $query.='transcript_id varchar(50) NOT NULL,';
    $query.='length int unsigned NOT NULL,';
    $query.='strand char(2) NOT NULL,';
    $query.='locus varchar(50) NOT NULL,';
    $query.='no_exons mediumint unsigned NOT NULL,';
    $query.='total mediumint unsigned NOT NULL,';
    foreach my $lane (@{$lanes}) {
	my $lane_id=$lane2id->{$lane};
	$lane_id=$lane;
	$query.="$lane mediumint unsigned NOT NULL,";
	$query.="INDEX idx_${lane} ($lane),";
    }
    $query.='INDEX idx_trans (transcript_id)';
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
    $query.='WHERE transcript_id = ?';
    $query.='GROUP BY transcript_id';
    $sth=$dbh->prepare($query);

    my $exon_number=sub{
	my $trans_id=shift;

	unless ($cache{$trans_id}) {
	    my $count=$sth->execute($trans_id);
	    if ($count !=0) {
		($cache{$trans_id})=$sth->fetchrow_array();
	    } else {
		warn "No hits for $trans_id\n";
		$cache{$trans_id}=0;
	    }
	}

	return($cache{$trans_id});
    };

    return($exon_number);
}

sub get_trans_expression {
    my $trans=shift;
    my $table=shift;
    my $dbh=shift;
    my $index=shift;
    my $lane=shift;
    my $threshold=shift || 1;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	
	$query ='SELECT transcript_id, rpkm ';
	$query.="FROM $table ";
	$query.='WHERE lane_id = ?';
	$sth=$dbh->prepare($query);
	$sth->execute($lane);
	
	while (my ($trans_id,$rpkm,$lane)=$sth->fetchrow_array()) {
	    my $expression=int($rpkm + 0.5);
	    if ($expression >= $threshold) {
		$trans->{$trans_id}->[$index]=$expression
	    }
	}
    }
}


