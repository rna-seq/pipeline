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
# This script should take a pair of gene ids and it will calculate the
# correlation between thes two genes accross all the datasetes in the database.
# The value used for calulating the correlation will be provided by using the
# suffix of the table containing it

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file',
			       'get_trans_expression_data');
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
my $tabsuffix='transcript_expression_levels_pooled';
my @trans_needed;

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'breakdown|b' => \$breakdown,
	   'trans|t=s' => \@trans_needed);

if ($breakdown) {
    $tabsuffix='transcript_expression_levels';
}

unless (@trans_needed >= 2) {
    die "I need at least 2 transcripts to compare (provided with the -t option\n";
}

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_ind_trans.RNAseqComp.log',
		      $debug);
print $log_fh "Extracting expression info for the following transcripts: ",
    join(",",@trans_needed),"\n";

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
			1)};

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
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment,2);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_trans_expression_data($dbh,
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
my $tmpfn="Trans.Expression.subset.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",
		  'Dataset',
		  'Project',
		  'Group',
		  @trans_needed),"\n";
foreach my $exp (@values) {
    my @row;
    foreach my $trans (@trans_needed) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$trans})) {
	    $value=$exp->[1]->{$trans};
	}
	push @row,$value;
    }

    my $dataset=$exp->[0];
    $dataset=~s/_gene_RPKM_pooled_sample//;
    my ($project,$group)=(split('_',$dataset))[0,-1];
    print $tmpfh join("\t",
		      $dataset,
		      $project,
		      $group,
		      @row),"\n";
}
close($tmpfh);

exit;
