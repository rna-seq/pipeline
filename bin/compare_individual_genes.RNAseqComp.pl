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
			       'get_gene_RPKM_data');
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_samples');
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $all=0;
my $debug=1;
my $breakdown;
my $tabsuffix='gene_RPKM_pooled';
my $genefile;
my @genes_needed;

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'breakdown|b' => \$breakdown,
	   'genefile|f=s' => \$genefile,
	   'gene|g=s' => \@genes_needed,
	   'all' => \$all); # get info for all experiments regardless of the
                            # project

if ($breakdown) {
    $tabsuffix='gene_RPKM';
}

if ($genefile) {
    my $genefh=get_fh($genefile);
    while (my $line=<$genefh>) {
	chomp($line);
	my $id=$line;
	push @genes_needed, $id;
    }
    close($genefh);
}

unless (@genes_needed >= 2) {
    die "I need at least 2 genes to compare (provided with the -g option\n";
}

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_ind_genes.RNAseqComp.log',
		      $debug);
print $log_fh "Extracting expression info for the following genes: ",
    join(",",@genes_needed),"\n";

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
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment,2);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_gene_RPKM_data($dbh,
				$table,
				\%all_genes,
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
my $tmpfn="Gene.Expression.subset.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",
		  'Dataset',
#		  'Project',
#		  'Group',
		  @genes_needed),"\n";
foreach my $exp (@values) {
    my @row;
    foreach my $gene (@genes_needed) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	}
	push @row,$value;
    }

    my $dataset=$exp->[0];
    $dataset=~s/_gene_[^_]+_pooled_sample//;
    my ($project,$group)=(split('_',$dataset))[0,-1];
    print $tmpfh join("\t",
		      $dataset,
#		      $project,
#		      $group,
		      @row),"\n";
}
close($tmpfh);

my $outfile='Expression.csv';
my $outfh=get_fh($outfile,1);
foreach my $exp (@values) {
    my $dataset=$exp->[0];
    $dataset=~s/_gene_[^_]_pooled_sample//;
    my ($project,$group)=(split('_',$dataset))[0,-1];
    foreach my $gene (@genes_needed) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	}
	print $outfh join("\t",
			  $dataset,
			  $project,
			  $group,
			  $gene,
			  $value),"\n";
    }
}
close($outfh);

exit;
