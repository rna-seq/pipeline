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
			       'get_dbh','get_trans_info_sub',
			       'get_trans_expression_data');
use RNAseq_pipeline_comp3 ('get_samples');

# Declare some variables
my %transcripts;
my %expressed;
my $transclasstable;
my $transexptable;
my $project_id;
my $exp_id;
my $threshold=1;
my $debug=1;

my $log_fh=get_log_fh('build_detected_transcripts.RNAseq.log',
		      $debug);

my %options=%{read_config_file()};
$project_id=$options{'PROJECTID'};
$exp_id=$options{'EXPID'};
$transexptable=$project_id.'_'.$exp_id.'_transcript_expression_levels_pooled';
$transclasstable=$options{'TRANSCLASSTABLE'};

my $dbh=get_dbh();
my $common_dbh=get_dbh(1);

# get some subs
*get_trans_info=get_trans_info_sub('type','status');

# Get the samples we will be looking at:
my %tables=($exp_id => $transexptable);
my %samples=%{get_samples(\%tables,
			  $dbh)};

my %detected;

# Get the gene information
my $trans_out=$project_id.'_'.$exp_id.'_detected_transcripts.txt';
foreach my $exp (keys %samples) {
    my ($table,$sample)=split('_sample_',$exp);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_trans_expression_data($dbh,
				       $table,
				       \%detected,
				       $sample);

    foreach my $trans_id (keys %{$data}) {
	my $type=get_trans_info($trans_id);
	$transcripts{$sample}{$type->[1]}{$type->[0]}++;
    }
}

# Print out the results
my $outfh=get_fh($trans_out,1);
foreach my $sample (sort keys %transcripts) {
    foreach my $status (sort keys %{$transcripts{$sample}}) {
	foreach my $type (sort keys %{$transcripts{$sample}{$status}}) {
	    print $outfh join("\t",
			      $type,
			      $status,
			      $sample,
			      $transcripts{$sample}{$status}{$type}),"\n";
	}
    }
}
close($outfh);
close($log_fh);

exit;

