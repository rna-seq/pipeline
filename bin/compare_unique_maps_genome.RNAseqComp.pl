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
# This script should take all the genome mappings from a project and build
# a table containing the number of reads mapped uniquely to each chromosome in
# each of the datasets

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 qw(get_dbh read_config_file);
use RNAseq_pipeline_stats3 qw(log10);
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_filenames');
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $tabsuffix='unique_maps_genome';
my $fraction='';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'limit|l=s' => \$fraction);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_unique_maps_genome.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			$fraction)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# For each of tables extract the number of uniquely mapping reads to each of
# the chromosomes
my %samples=%{get_filenames(\%tables,
			   $dbh)};
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_mapping_data($dbh,
			      $table,
			      \%all_genes,
			      $sample);

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

# Print the expression values for each gene in each of the tables into a
# temporary file
my $tmpfn="Unique.maps.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",@experiments),"\n";
foreach my $gene (keys %all_genes) {
    my @row;
    my $no_print=0;
    foreach my $exp (@values) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	} else {
#	    $no_print=1;
	}
	push @row,$value;
    }

    unless ($no_print) {
	$gene='chr'.$gene;
	print $tmpfh join("\t",
			  $gene,
			  @row),"\n";
    }
}
close($tmpfh);

# Plot the correlation graph if we have more than 2 samples
if (@experiments > 2) {
    my $outfn=$project.".unique.clusters";
    plot_graphs_R($tmpfn,
		  $outfn);
} else {
    print STDERR "Not enough samples to cluster\n";
}

exit;

sub get_mapping_data {
    my $dbh=shift;
    my $table=shift;
    my $all=shift;
    my $sample=shift;

    my %expression;

    my ($query,$sth,$count);
    $query ='SELECT chromosome, uniqueMaps ';
    $query.="FROM $table ";
    $query.='WHERE filename = ?';
    $sth=$dbh->prepare($query);
    $count=$sth->execute($sample);
    
    if ($count && ($count > 1)) {
	print STDERR $count,"\tChromosomes are detected in $table\n";
    } else {
	die "No chromosomes present in $table\n";
    }

    # get all the necessary tables
    while (my ($chr,$maps)=$sth->fetchrow_array()) {
	$chr=~s/MT/M/;
	$expression{$chr}=$maps;
	$all->{$chr}=1;
    }

    return(\%expression);
}

# Build a postscript tree using the  canberra distance and the complete linkage
# for clustering
sub plot_graphs_R {
    my $statsfn=shift;
    my $outfile=shift;

    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;

    # Read the data
    $r_string.="rpkms<-read.table(\"$statsfn\",sep=\"\t\",header=T)\n";

    # Calculate the distance matrix
    $r_string.="genesdist<-dist(t(rpkms),method='canberra')\n";

    # Setup the figure
    $r_string.="pdf(\"$outfile.pdf\")\n";

    # Build the tree
    $r_string.="plot(hclust(genesdist))\n";
    $r_string.="dev.off()\n";
	
    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet < $execution_file";
    run_system_command($command);

    $command="rm $execution_file";
    run_system_command($command);
}

