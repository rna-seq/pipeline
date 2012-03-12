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
# This script should take a pair of gene ids and it will extract for each of
# Them the expression values of each of the corresponding transcripts
# Also correlation between these two genes across all the datasets in the
# database will be calculated.
# The value used for calulating the correlation will be provided by using the
# suffix of the table containing it

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file','get_dbh',
			       'get_junc_expression_data',
			       'get_gene_from_short_junc_sub');
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_samples');
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $breakdown;
my $transfile;
my $tabsuffix='all_junctions_class_pooled';
my @juncs_needed;
my $all=0;
my $normalize; # Normalize by gene (get the EJEI basically)

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'breakdown|b' => \$breakdown,
	   'transcriptfile|f=s' => \$transfile,
	   'trans|t=s' => \@juncs_needed,
	   'all|a' => \$all,
	   'norm' => \$normalize);

if ($breakdown) {
    $tabsuffix='transcript_expression_levels';
}

# Get subs
*gene2trans=gene2junc_sub();
*junc2gene=get_gene_from_short_junc_sub();

my %trans2gene;
if ($transfile) {
    my $transfh=get_fh($transfile);
    while (my $line=<$transfh>) {
	chomp($line);
	my ($id,$type)=split("\t",$line);
	if ($type &&
	    $type eq 'gene') {
	    # Get all the junctions from the gene of interest
	    my $junctions=gene2trans($id);
	    my @transcripts;
	    if ($junctions) {
		@transcripts=@{$junctions};
		push @juncs_needed,@transcripts;
	    }
	} else {
	    push @juncs_needed, $id;
	}
    }
    close($transfh);
}

unless (@juncs_needed >= 2) {
    die "I need at least 2 transcripts to compare (provided with the -t option\n";
}

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_ind_trans.RNAseqComp.log',
		      $debug);
print $log_fh "Extracting expression info for the following transcripts: ",
    join(",",@juncs_needed),"\n";

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);

# Get all experiment tables from the database
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			'',
			$all)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# For each of tables extract the RPKMs of interest
my %samples=%{get_samples(\%tables,
			  $dbh,
			  $breakdown)};
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_trans;
my %gene_exp_junc;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment,2);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_junc_expression_data($dbh,
				       $table,
				       \%all_trans,
				       $sample,
				       $breakdown);
    if ($data) {
	push @values, [$experiment,$data];
    } else {
	print STDERR "Skipping $experiment\n";
    }
}

# Get the human readable lables
foreach my $experiment (@experiments) {
    my $label;
    if ($nolabels) {
	$label=$samples{$experiment}->[1];
    } else {
	$label=get_labels($experiment);
    }
    if ($label) {
	$experiment=$label;
    }
}

# Print the expression values for each gene of interest
my $tmpfn="Junc.ReadCounts.subset.txt";
my $tmpfh=get_fh($tmpfn,1);
foreach my $exp (@values) {
    my @row;
    my $dataset=$exp->[0];
    $dataset=~s/_transcript_expression_levels_pooled_sample//;
    my ($project,$group)=(split('_',$dataset))[0,-1];
    foreach my $trans (@juncs_needed) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$trans})) {
	    $value=$exp->[1]->{$trans};
	}
	my @genes=@{junc2gene($trans)};

	foreach my $gene_id (@genes) {
	    print $tmpfh join("\t",
			      $dataset,
			      $project,
			      $group,
			      $gene_id,
			      $trans,
			      $value),"\n";
	}
    }
}
close($tmpfh);

exit;

sub gene2junc_sub {
    my %options=%{read_config_file()};
    my $table=$options{'JUNCTIONSTABLE'} || die "No junctions table defined\n";

    my $dbh=get_dbh(1);
    my %cache;

    my ($query,$sth);
    $query ='SELECT distinct chr,start,end ';
    $query.="FROM $table ";
    $query.='WHERE gene_id = ? and type= "known"';
    $sth=$dbh->prepare($query);
    
    my $gene2trans=sub {
	my $gene=shift;

	unless ($cache{$gene}) {
	    $sth->execute($gene);
	    while (my ($chr,$start,$end)=$sth->fetchrow_array()) {
		my $junc_id=join("_",$chr,$start,$end);
		push @{$cache{$gene}},$junc_id;
	    }
	}
	return($cache{$gene});
    };

    return($gene2trans);
}
