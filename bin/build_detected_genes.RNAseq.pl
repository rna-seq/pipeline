#!/soft/bin/perl

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

# Author : David Gonzalez, david.gonzalez@crg.eu

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
# This script will use the database in order to determine what features have
# been detected. Each gene is checked for the type and reliability and a table
# is built from this data summarizing it

use RNAseq_pipeline3 qw(get_fh get_log_fh);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_dbh','get_gene_info_sub','get_gene_RPKM_data');
use RNAseq_pipeline_comp3 ('get_samples');

# Declare some variables
my %genes;
my %expressed;
my $geneclasstable;
my $generpkmtable;
my $project_id;
my $exp_id;
my $threshold=1;
my $debug=1;

my $log_fh=get_log_fh('build_detected_genes.RNAseq.log',
		      $debug);

my %options=%{read_config_file()};
$project_id=$options{'PROJECTID'};
$exp_id=$options{'EXPID'};
$generpkmtable=$project_id.'_'.$exp_id.'_gene_RPKM_pooled';
$geneclasstable=$options{'GENECLASSTABLE'};

my $dbh=get_dbh();
my $common_dbh=get_dbh(1);

# get some subs
*get_gene_info=get_gene_info_sub('type','status');

# Get the samples we will be looking at:
my %tables=($exp_id => $generpkmtable);
my %samples=%{get_samples(\%tables,
			  $dbh)};

my %detected;

# Initialize the values with a total
$genes{'total'}{'total'}{'total'}=0;

# Get the gene information
my $gene_out=$project_id.'_'.$exp_id.'_detected_genes.txt';
foreach my $exp (keys %samples) {
    my ($table,$sample)=split('_sample_',$exp);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_gene_RPKM_data($dbh,
				$table,
				\%detected,
				$sample);

    foreach my $gene_id (keys %{$data}) {
	# Filter by expression
	if ($data->{$gene_id} < $threshold) {
	    next;
	}
	my $type=get_gene_info($gene_id);
	$genes{$sample}{$type->[1]}{$type->[0]}++;
	$genes{'total'}{'total'}{'total'}++;
    }
}

# Print out the results
my $outfh=get_fh($gene_out,1);
foreach my $sample (sort keys %genes) {
    foreach my $status (sort keys %{$genes{$sample}}) {
	foreach my $type (sort keys %{$genes{$sample}{$status}}) {
	    print $outfh join("\t",
				 $type,
				 $status,
				 $sample,
				 $genes{$sample}{$status}{$type}),"\n";
	}
    }
}
close($outfh);
close($log_fh);

exit;

