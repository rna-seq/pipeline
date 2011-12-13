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
# This script should take a project from the database and for each of the
# experiments belonging to it build a comparison table that can be run
# through R
# We have to be able to decide what table to get the data from (pooled or not)
# And also we have to be able to select genes that are expressed 
# In order to get set the descriptions after using the script get_EnsEMBL_gene_info.pl and a list of genes:
# gawk -F"\t" '{print "UPDATE 1_H_sapiens_EnsEMBL_55_parsed_gt_geneclass set description=\""$4"\" WHERE gene_id=\""$1"\";"}' all.genes.desc.txt > add.description.sql


use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file',
			       'get_exon_readcount_data','get_gene_info_sub',
			       'get_gene_from_exon_sub');
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
my $tabsuffix='exon_RPKM_pooled';
my $fraction='';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'limit|l=s' => \$fraction);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_exon_RPKM.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);
*gene2desc=get_gene_info_sub('description');
*gene2type=get_gene_info_sub('type');
*gene2chr=get_gene_info_sub('chr');
*exon2gene=get_gene_from_exon_sub();

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			$fraction)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# For each of tables extract the RPKMs of interest and get for each of the
# tables the different samples present in them
my %samples=%{get_samples(\%tables,
			  $dbh,
			  0)};# currently set to one until we fix the table naming problem
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_exon_readcount_data($dbh,
				     $table,
				     \%all_genes,
				     $sample,
				     1); # currently set to one until we fix the table naming problem
    if ($data) {
	push @values, [$experiment,$data];
    } else {
	print STDERR "Skipping $experiment\n";
    }
}

# Get the human readable lables
foreach my $experiment (@experiments) {
    my $label;
    $label=$samples{$experiment}->[1];
}

# Print the expression values for each gene in each of the tables into a
# temporary file
my $tmpfn="Exon.ReadCount.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",@experiments),"\n";
foreach my $exon (keys %all_genes) {
    my @genes=@{exon2gene($exon)};
    my $gene='';
    if (@genes == 1) {
	($gene)=@genes;
    } else {
	die "More than one gene found for $exon\n";
    }

    my @row;
    my $no_print=0;
    foreach my $exp (@values) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	}
	push @row,$value;
    }

    # Skip mitochondrial and ribosomal genes
    my $desc=gene2desc($gene) || '';
    if (gene2chr($gene)=~/chrM/o) {
	next;
    } elsif ($desc=~/ribosom(e|al)/io) {
	next;
    } elsif (gene2type($gene)=~/^rRNA/o) {
	next;
    } else {
	print $tmpfh join("\t",
			  $gene,
			  @row),"\n";
    }
}
close($tmpfh);

exit;
