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

# Objective:
# This script should create a gtf file with the requested features (genes,
# transcripts or exons), adding the RPKM values to them

# Load modules
use Getopt::Long;
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file',
			       'get_gene_RPKM_data',
			       'get_trans_expression_data',
			       'get_exon_readcount_data');
use RNAseq_pipeline3 ('get_fh','parse_gff_line');

# Declare some variables:
my $annotation;
my $feature='gene';
my $prefix;
my %detected;
my %rpkms;
my $dbh;

# Get command line options
GetOptions('feature|f=s' => \$feature);

# read the config file
my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$prefix=$options{'PREFIX'};

print STDERR "Extracting $feature features from $annotation\n";

# Get a dbh
$dbh=get_dbh();

if ($feature eq 'gene') {
    my $tabsuffix='_gene_RPKM_pooled';
    my $table=$prefix.$tabsuffix;
    %rpkms=%{get_gene_RPKM_data($dbh,
				$table,
				\%detected)};
    my $annotfh=get_fh($annotation);
    while (my $line=<$annotfh>) {
	chomp($line);
	if ($line=~/^#/) {next};
	my %line=%{parse_gff_line($line)};
	if ($line{'type'} ne 'gene') {next;}
	my $gene_id=$line{'feature'}{'gene_id'};
	if ($rpkms{$gene_id}) {
	    print $line," RPKM \"$rpkms{$gene_id}\"\n";
	}
    }
    close($annotfh);
} elsif ($feature eq 'transcript') {
    my $tabsuffix='_transcript_expression_levels_pooled';
    my $table=$prefix.$tabsuffix;
    %rpkms=%{get_trans_expression_data($dbh,
				       $table,
				       \%detected)};
    my $annotfh=get_fh($annotation);
    while (my $line=<$annotfh>) {
	chomp($line);
	if ($line=~/^#/) {next};
	my %line=%{parse_gff_line($line)};
	if ($line{'type'} ne 'transcript') {next;}
	my $trans_id=$line{'feature'}{'transcript_id'};
	if ($rpkms{$trans_id}) {
	    print $line," RPKM \"$rpkms{$trans_id}\"\n";
	}
    }
    close($annotfh);
} elsif ($feature eq 'exon') {
    my $tabsuffix='_exon_RPKM_pooled';
    my $table=$prefix.$tabsuffix;
    %rpkms=%{get_exon_readcount_data($dbh,
				     $table,
				     \%detected)};
    my $annotfh=get_fh($annotation);
    while (my $line=<$annotfh>) {
	chomp($line);
	if ($line=~/^#/) {next};
	my %line=%{parse_gff_line($line)};
	if ($line{'type'} ne 'exon') {next;}
	my $strand=1;
	if ($line{'strand'} eq '-') {
	    $strand=-1;
	}
	my $exon_id=join('_',
			 $line{'chr'},
			 $line{'start'},
			 $line{'end'},
			 $strand);
	if ($rpkms{$exon_id}) {
	    print $line," RPKM \"$rpkms{$exon_id}\"\n";
	}
    }
    close($annotfh);
} else {
    die "Unknown feature $feature\n";
}

exit
