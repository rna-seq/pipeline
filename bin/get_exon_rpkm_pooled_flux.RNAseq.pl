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
# pooled and it will generate the exon RPKM values by adding the RPKM values
# as calculated by the flux capacitor for each of the transcripts that contain
# the eon

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','print_table_file');

# Declare some variables
my $prefix;
my $table;
my $exontab;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$table=$prefix.'_transcript_expression_levels_pooled';
$exontab=$options{'EXONSCLASSTABLE'};

# Connect to the database
my $dbh=get_dbh();
my $dbhcommon=get_dbh(1);

# First get all exons and for each exon all transcripts in which it is
# contained
my $exons=get_all_exon_transcripts($dbhcommon,
				   $exontab);

# Get the transcript expression levels;
my $trans_exp=get_transcript_expression($dbh,
					$table);

# Get the exon expression levels;
my %exon_exp=%{get_exon_expression($exons,
				   $trans_exp)};

# Print out the results
foreach my $sample (keys %exon_exp) {
    foreach my $exon (keys %{$exon_exp{$sample}}) {
	my $expression=sprintf "%.2f",$exon_exp{$sample}{$exon};
	print join("\t",
		   $exon,
		   $expression,
		   $sample),"\n";
    }
}

# Print out the table file
my $outtable=$prefix.'_exon_RPKM_pooled_flux';
my $outtablesql="CREATE TABLE $outtable (
    exon_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_exon (exon_id),
    INDEX idx_sample (sample)
    );";
print_table_file($outtable,
		 $outtablesql);

exit;

sub get_exon_expression {
    my $exons=shift;
    my $transcripts=shift;

    my %expression;

    # Add up for each exon the expression of every transcript that contains it
    foreach my $exon (keys %{$exons}) {
	foreach my $trans (keys %{$exons->{$exon}}) {
	    if ($transcripts->{$trans}) {
		my @samples=keys %{$transcripts->{$trans}};
		foreach my $sample (@samples) {
		    $expression{$sample}{$exon}+=$transcripts->{$trans}->{$sample};
		}
	    } else {
		next;
	    }
	}
    }

    return(\%expression);
}


sub get_transcript_expression {
    my $dbh=shift;
    my $table=shift;

    my %expression;

    my ($query,$sth);
    $query ='SELECT transcript_id, rpkm, sample ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($trans,$rpkm,$sample)=$sth->fetchrow_array()) {
	$expression{$trans}{$sample}+=$rpkm;
    }

    return(\%expression);
}

sub get_all_exon_transcripts {
    my $dbh=shift;
    my $table=shift;

    my %exons;

    my ($query,$sth);
    $query ='SELECT exon_id,transcript_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($exon,$transcript)=$sth->fetchrow_array()) {
	$exons{$exon}{$transcript}=1;
    }

    return(\%exons);
}
