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
# This script will take the information from the transcript expression levels
# pooled and it will generate the gene RPKM values by adding the RPKM values
# for each of the transcript corresponding to that gene as calculated by the
# flux capacitor

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh');

# Declare some variables
my $prefix;
my $table;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$table=$prefix.'_transcript_expression_levels_pooled';

# Connect to the database
my $dbh=get_dbh();

# Get the gene expression levels;
my %gene_exp=%{get_gene_expression($dbh,
				   $table)};

# Print out the results
foreach my $sample (keys %gene_exp) {
    foreach my $gene (keys %{$gene_exp{$sample}}) {
	my $expression=sprintf "%.2f",$gene_exp{$sample}{$gene};
	print join("\t",
		   $gene,
		   $expression,
		   $sample),"\n";
    }
}

exit;

sub get_gene_expression {
    my $dbh=shift;
    my $table=shift;

    my %expression;

    my ($query,$sth);
    $query ='SELECT gene_id, rpkm, sample ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($gene,$rpkm,$sample)=$sth->fetchrow_array()) {
	$expression{$sample}{$gene}+=$rpkm;
    }

    return(\%expression);
}
