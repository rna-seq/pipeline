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
# This scritp will take information from the tables containing the inclusion
# data and it will build the distribution of the inclsuion levels

use RNAseq_pipeline3 qw(get_fh check_table_existence);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file');

# Declare some variables
my $prefix;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};

# Get the table names we need
my $incl=$prefix.'_exon_inclusion';
my $incl_pooled=$prefix.'_exon_inclusion_pooled';

# Get the information for the histogram
my %hist;
my $dbh=get_dbh();
if (check_table_existence($dbh,$incl_pooled)) {
#    get_info_from_table($dbh,
#			$incl,
#			%hist);
#} elsif (check_table_existence($dbh,$incl_pooled)) {
    get_info_from_table($dbh,
			$incl_pooled,
			%hist);
} else {
    die "No suitable exon inclusion table found\n";
}

# Print out the result
foreach my $value (keys %hist) {
    foreach my $set (keys %{$hist{$value}}) {
	print join("\t",
		   $hist{$value}{$set},
		   $value,
		   $set),"\n";
    }
}

exit;

sub get_info_from_table {
    my $dbh=shift;
    my $table=shift;
    my $hist=shift;

    my ($query,$sth,$count);
    $query ='SELECT * ';
    $query.="FROM $table ";
    $query.='WHERE inc_rate IS NOT NULL';
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless($count && 
	   ($count > 0)) {
	die "No entries found in $table\n";
    }

    while (my ($gene_id,$exon_id,$ExIncl,$JuncInc,
	       $JuncExc,$inc,$set)=$sth->fetchrow_array()) {
	my $value=int($inc * 100);
	if (defined $value) {
	    $hist{$value}{$set}++;
	}
    }
}
