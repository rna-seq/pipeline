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
# This script should get for every combination of lanes the exon_inclusion
# correlation

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_lanes','get_dbh');

# Declare some variables
my %inclusion;
my %expressed;

my $exon_inclusion_table;
my $threshold=1;

my %options=%{read_config_file()};
$exon_inclusion_table=$options{'PREFIX'}.'_exon_inclusion';

my $dbh=get_dbh();

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};

# Get the number of detected features in each lane
my @lanes=sort keys %lanes;

# Get the gene hits
for (my $i=0;$i<@lanes;$i++) {
    get_exon_inclusion(\%inclusion,
		       $exon_inclusion_table,
		       $dbh,
		       $i,
		       $lanes[$i],
		       $threshold);
}

# Build the inclusion table
foreach my $exon (keys %inclusion) {
    for (my $i=0;$i<@lanes;$i++) {
	$expressed{$exon}->[$i]=$inclusion{$exon}->[$i] || 0;
    }
}

# print out the results
foreach my $exon (keys %expressed) {
    for (my $i=0;$i<@lanes;$i++) {
	for (my $j=0;$j<@lanes;$j++) {
	    if ($i==$j) {
		next;
	    }
	    my $inc1_exp=$expressed{$exon}->[$i] || 0;
	    my $inc2_exp=$expressed{$exon}->[$j] || 0;
	    print join("\t",
		       $lanes[$i],
		       $lanes[$j],
		       $inc1_exp,
		       $inc2_exp),"\n";
	}
    }
}

exit;

sub get_exon_inclusion {
    my $exons=shift;
    my $table=shift;
    my $dbh=shift;
    my $index=shift;
    my $lane=shift;
    my $threshold=shift || 1;

    my ($query,$sth);

    $query ='SELECT exon_id, inc_rate ';
    $query.="FROM $table ";
    $query.='WHERE lane_id = ? AND inc_rate IS NOT NULL ';
    $query.='AND (ExIncl + JuncInc + JuncExc) >= ?';
    $sth=$dbh->prepare($query);
    $sth->execute($lane,5);

    while (my ($exon,$inc_rate)=$sth->fetchrow_array()) {
	my $expression=int(($inc_rate * 100) + 0.5);
	$exons->{$exon}->[$index]=$expression
    }
}

