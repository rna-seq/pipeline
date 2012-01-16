#!/soft/bin/perl
# DGK

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

#    Author : David Gonzalez, david.gonzalez@crg.eu

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
# through R This table will contain for each junction the number of reads
# detected.
# We have to be able to decide what table to get the data from (pooled or not)
# And also we have to be able to select genes that are expressed 

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command get_list);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file','get_gene_info_sub',
			       'get_gene_from_short_junc_sub',
			       'get_gene_from_exon_sub');
use RNAseq_pipeline_stats3 qw(log10);
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_samples','remove_tables');
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $tabsuffix='fusion_transcripts';
my $limit;
my $threshold=1; # Minimum support required
my $subset;

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'subset|s=s' => \$subset,
	   'limit|l=s' => \$limit,
	   'threshold|t=i' => \$threshold);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_fusions.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);
*gene2chr=get_gene_info_sub('chr');
*gene2desc=get_gene_info_sub('description');
*gene2type=get_gene_info_sub('type');

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			$limit)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# If a subset has been provided remove any tables that are not included in the
# subset
if ($subset && -r $subset) {
    my %subset=%{get_list($subset)};
    remove_tables(\%tables,
		  \%subset);
}

# For each of tables extract the splice site support of interest and get for
# each of the tables the different samples present in them
my %samples=%{get_samples(\%tables,
			  $dbh,
			  0)};
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_fusion_data($dbh,
			       $table,
			       \%all_genes,
			       $sample);
    if ($data) {
	push @values, [$experiment,$data];
    } else {
	print STDERR "Skipping $experiment\n";
    }
}

# Get the human readable labels
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

# Print the detected reads for each junction in each of the tables into a
# temporary file
my $tmpfn="Pair.Fusions.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",@experiments),"\n";
foreach my $gene (keys %all_genes) {
    my @row;
    my $no_print=0;
    my $total=0;
    foreach my $exp (@values) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	    $total++;
	}
	push @row,$value;
    }

    my @gene_names=split('_',$gene);
    my $gene_name=join('_',sort @gene_names);
    my $print=0;

    foreach my $gene_id (@gene_names) {
	my $desc=gene2desc($gene_id) || '';
	if (gene2chr($gene_id)=~/chrM/o) {
	    next;
	} elsif ($desc=~/ribosom(e|al)/io) {
	    next;
	} elsif (gene2type($gene_id)=~/^rRNA/o) {
	    next;
	} else {
	    $print=1;
	    last;
	}
    }
    if ($print) {
	print $tmpfh join("\t",
			  $gene_name.'_'.$gene,
			  @row),"\n";
    }
}
close($tmpfh);

exit;

sub get_fusion_data {
    my $dbh=shift;
    my $table=shift;
    my $all=shift;
    my $sample=shift;

    my %expression;

    my ($query,$sth,$count);
    $query ='SELECT gene1, gene2, chr1, chr2, genomic1mapstart, genomic2mapstart ';
    $query.="FROM $table ";
    $query.='WHERE LaneName = ?';
#    $query.=' AND junc_type not like "split%"';
#    $query.=' limit 100';
    $sth=$dbh->prepare($query);
    $count=$sth->execute($sample);
    
    if ($count && ($count > 1)) {
	print STDERR $count,"\tFusions are detected in $table\n";
    } else {
	die "No fusions present in $table\n";
    }

    # get all the necessary tables
    while (my ($gene1,$gene2,$chr1,$chr2,$start1,$start2)=$sth->fetchrow_array()) {
	my $splice_id=join('_',
			   sort($gene1,$gene2));
	$expression{$splice_id}++;
	$all->{$splice_id}=[$chr1,$start1,$chr2,$start2];
    }

    return(\%expression);
}
