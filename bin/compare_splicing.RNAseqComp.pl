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
my $tabsuffix='all_junctions_class_pooled';
my $limit;
my $threshold=1; # Minimum support required
my $subset;
my $samplesthreshold=1; # Minimum samples in which a junction must be present in
                        # order to consider it

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'subset|s=s' => \$subset,
	   'limit|l=s' => \$limit,
	   'threshold|t=i' => \$threshold,
	   'samples=i' => \$samplesthreshold,);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_splicing.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);
*junc2gene=get_gene_from_short_junc_sub();
*ex2gene=get_gene_from_exon_sub();
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
    my $data=get_splicing_data($dbh,
			       $table,
			       \%all_genes,
			       $sample,
			       $threshold);
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
my $tmpfn="Junctions.$project.txt";
my $fusionsfn="Fusions.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
my $fusionsfh=get_fh($fusionsfn,1);
print $tmpfh join("\t",@experiments),"\n";
print $fusionsfh join("\t",@experiments),"\n";
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

    # Skip those cases present in lesss than $threshold samples
    if ($total < $samplesthreshold) {
	next;
    }

    my @junctions=@{junc2gene($gene)};
    my $gene_id=join('_',@junctions);
    if ($gene_id eq 'Unknown') {
	my @exons=@{$all_genes{$gene}};
	my @gene_names;
	foreach my $exon (@exons) {
	    if ($exon=~/-/) {
		push @gene_names, 'Unknown';
	    } else {
		my @exons=split(';',$exon);
		my %genes;
		foreach my $frag (@exons) {
		    my ($gene_name)=@{ex2gene($frag)};
		    $genes{$gene_name}++;
		}
		if (keys %genes > 1) {
		    print $log_fh "More than one gene corresponds to $exon\n";
		    print $log_fh join("\t",keys %genes),"\n";
		} else {
		    push @gene_names, keys %genes;
		}
	    }
	}
	if (@gene_names) {
	    $gene_id=join('_',sort @gene_names);
	    print $fusionsfh join("\t",
				  $gene_id.'_'.$gene,
				  @row),"\n";
	}
    } else {
	my $print=0;
	foreach my $gene_id (@junctions) {
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
			      $gene_id.'_'.$gene,
			      @row),"\n";
	}
    }
}
close($tmpfh);
close($fusionsfh);

exit;

sub get_splicing_data {
    my $dbh=shift;
    my $table=shift;
    my $all=shift;
    my $sample=shift;
    # Default limit set to 2 in order to prevent the huge datasets from eating
    # up RAM
    my $limit=shift || 1;

    my %expression;

    my ($query,$sth,$count);
    $query ='SELECT chr1, start, chr2, end, support, exons1, exons2 ';
    $query.="FROM $table ";
    $query.='WHERE sample = ? and support >= ?';
#    $query.=' AND junc_type not like "split%"';
#    $query.=' limit 100';
    $sth=$dbh->prepare($query);
    $count=$sth->execute($sample,$limit);
    
    if ($count && ($count > 1)) {
	print STDERR $count,"\tJunctions are detected in $table\n";
    } else {
	die "No junctions present in $table\n";
    }

    # get all the necessary tables
    while (my ($chr1,$start,$chr2,$end,$support,$ex1,$ex2)=$sth->fetchrow_array()) {
	my $splice_id=join('_',
			   $chr1,$start,$end);
	if ($chr1 ne $chr2) {
	    $splice_id=join('_',
			    $chr1,$start,$chr2,$end);
	}
	$expression{$splice_id}=$support;
	$all->{$splice_id}=[$ex1,$ex2];
    }

    return(\%expression);
}
