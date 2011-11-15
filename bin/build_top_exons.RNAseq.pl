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
# This script will take as input an annotation file.
# it will read the exon_RPKM table for the dataset in question and it will get
# a list of the exons that are expressed.
# for those exons that are not expressed it will fill in with zeros
# after it will read the annotation file and it will generate a table with all
# the necessary information for the exon from the annotation, for all those
# exons that are detected in at least one condition

use RNAseq_pipeline3 qw(get_fh parse_gff_line check_table_existence);
use RNAseq_pipeline_settings3 ('read_config_file','get_dataset_id',
			       'get_dbh');

# Declare some variables
my %exons;
my @files;
my %expressed;

my $annotation;
my $prefix;
my $exon_rpkm_table;
my $top_exons_table;
my $datasets_table;
my $db;
my $threshold=1;

my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$prefix=$options{'PREFIX'};
$exon_rpkm_table=$prefix.'_exon_RPKM';
$top_exons_table=$prefix.'_top_exons_expressed';
$datasets_table=$prefix.'_dataset';
$db=$options{'DB'};

my $dbh=get_dbh();

my %distribution;

my %lanes=%{get_dataset_id($dbh,
			   $datasets_table)};

# Get the number of detected features in each lane
my @lanes=sort keys %lanes;

# Get the gene hits
for (my $i=0;$i<@lanes;$i++) {
    get_exon_expression(\%exons,
			$exon_rpkm_table,
			$dbh,
			$i,
			$lanes[$i],
			$threshold);
}

# Build the expression table
foreach my $exon (keys %exons) {
    for (my $i=0;$i<@lanes;$i++) {
	$expressed{$exon}->[$i]=$exons{$exon}->[$i] || 0;
    }
}

# Read the annotation and print out the table
my $annotfh=get_fh($annotation);
my %features;
while (my $line=<$annotfh>) {
    if ($line=~/^#/) {
	next;
    }
    my %line=%{parse_gff_line($line)};

    unless  ($line{'type'}=~/exon/) {
	next;
    }
    my $chr=$line{'chr'};
    my $start=$line{'start'};
    my $end=$line{'end'};
    my $strand;
    if ($line{'strand'} eq '+') {
	$strand=1;
    } elsif ($line{'strand'} eq'-') {
	$strand=-1;
    } else {
	die "Incorrect strand\n";
    }

    my $exon_id=join('_',
		     $chr,
		     $start,
		     $end,
		     $strand);

    if (defined $expressed{$exon_id}) {
	$features{$exon_id}->[0]=$line{'end'} - $line{'start'} + 1;
	$features{$exon_id}->[1]=$line{'strand'};
	$features{$exon_id}->[2]=$line{'chr'};
    }
}
close($annotfh);

# print out the results
my $taboutfh=get_fh($top_exons_table.'.txt',1);
foreach my $exon (keys %expressed) {
    my $total=0;
    foreach my $value (@{$expressed{$exon}}) {
	$total+=$value;
    }
    print $taboutfh join("\t",
			 $exon,
			 @{$features{$exon}},
			 $total,
			 @{$expressed{$exon}}),"\n";
}
close($taboutfh);

# Build the table in the database
build_db_table($dbh,
	       $top_exons_table,
	       \@lanes,
	       \%lanes);

# Insert into the database
my $command="mysql $db < $top_exons_table.sql";
print STDERR "Executing: $command\n";
system($command);
$command="mysqlimport -L $db $top_exons_table.txt";
print STDERR "Executing: $command\n";
system($command);
$command="rm $top_exons_table.sql $top_exons_table.txt";
print STDERR "Executing: $command\n";
system($command);

exit;

sub build_db_table {
    my $dbh=shift;
    my $table=shift;
    my $lanes=shift;
    my $lane2id=shift;
    my $outtable=$table.'.sql';

    my ($query,$sth);
    $query ="DROP TABLE IF EXISTS $table;";
    $query.="CREATE TABLE $table( ";
    $query.='exon_id varchar(50) NOT NULL,';
    $query.='length int unsigned NOT NULL,';
    $query.='strand char(2) NOT NULL,';
    $query.='locus varchar(50) NOT NULL,';
    $query.='total mediumint unsigned NOT NULL,';
    foreach my $lane (@{$lanes}) {
	my $lane_id=$lane2id->{$lane};
	$lane_id=$lane;
	$query.="$lane_id mediumint unsigned NOT NULL,";
	$query.="INDEX idx_${lane_id} ($lane_id),";
    }
    $query.='INDEX idx_exon (exon_id)';
    $query.=');';

    my $outfh=get_fh($outtable,1);
    print $outfh $query,"\n";
    close($outtable);
}

sub get_exon_expression {
    my $exons=shift;
    my $table=shift;
    my $dbh=shift;
    my $index=shift;
    my $lane=shift;
    my $threshold=shift || 1;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT exon_id, RPKM ';
	$query.="FROM $table ";
	$query.='WHERE LaneName = ?';
	$sth=$dbh->prepare($query);
	$sth->execute($lane);
	
	while (my ($exon,$rpkm,$lane)=$sth->fetchrow_array()) {
	    my $expression=int($rpkm + 0.5);
	    if ($expression >= $threshold) {
		$exons->{$exon}->[$index]=$expression
	    }
	}
    }
}


