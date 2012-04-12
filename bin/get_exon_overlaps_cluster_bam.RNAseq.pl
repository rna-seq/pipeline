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
# This script should run overlap on a set reads and exon projections

# Load some modules
use RNAseq_pipeline3 ('get_fh','get_annotation_from_gtf');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list');
use Getopt::Long;
use Bio::DB::Sam;
use Bio::Range;


# Declare some variables
my $threshold=1;
my $tmpdir;
my $clusterdir;
my $samdir;
my $exondir;
my $splitdir;
my $junctionsdir;
my $projectid;
my $prefix;
my $stranded;
my $file_list;
my $parallel='cluster';
my $paralleltmp;
my $bindir;

GetOptions('threshold|t=i' => \$threshold,
    );

# Read the options file
my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$exondir=$options{'PROJECT'}.'/'.$options{'EXONDIR'};
$samdir=$options{'PROJECT'}.'/'.$options{'SAMDIR'};
$splitdir=$options{'SPLITMAPDIR'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$projectid=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$stranded=$options{'STRANDED'};
$file_list=$options{'FILELIST'};
$paralleltmp=$options{'PROJECT'}.'/'.$options{'LOCALPARALLEL'};
$bindir=$options{'BIN'};
unless($options{'CLUSTER'}) {
    print STDERR "Running locally as no cluster has been defined\n";
    $parallel='default';
}

my $annotation=$exondir.'/'.$prefix.'.exon.gtf';
print STDERR $annotation,"\n";


# First get the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get for each lane the BAM files
my %lane_files=%{get_lane_files(\%files)};

my %coverage;

# Extract the genes from the annotation
my %exons=%{get_annotation_from_gtf($annotation,
				    '',
				    'exons')};

foreach my $pair (keys %lane_files) {
    foreach my $lane (keys %{$lane_files{$pair}}) {
	my $bamfile=$samdir.'/'.$pair.'.merged.bam';
	print STDERR "$bamfile\n";
	my $outfile=$exondir.'/'.$lane.'.single.unique.gtf.overlap.total';
	
	if (-r $outfile) {
	    print STDERR $outfile,"\tis present. Skipping\n";
	    next;
	}
	
	# Get the necessary filehandles
	my $sam = Bio::DB::Sam->new(-bam  =>$bamfile,
				    -autoindex => 1);
	my $outfh=get_fh($outfile,1);
	
	# Process the genes
	my %gene_hits;
	foreach my $exon (keys %exons) {
	    my @exon=split('_',$exon);
	    my $strand=pop(@exon);
	    my $end=pop(@exon);
	    my $start=pop(@exon);
	    my $chr=join('_',@exon);


	    # Get the exon information
	    my $exonobj=Bio::Range->new(-start=>$start,
					-end=>$end,
					-strand=>$strand);
	    my $feat_hits=process_feature($exonobj,
					  $chr,
					  $sam,
					  $exon);
	    
	    $gene_hits{'exon'}{$exon}=$feat_hits;

	
	    # Print the output
	    # only print if there are actually reads in the gene

	    if ($gene_hits{'exon'}{$exon}->[2]) {
		my $length=$exonobj->end() - $exonobj->start() + 1;
		my $frac=sprintf "%.2f",$gene_hits{'exon'}{$exon}->[2]/$length;
		print $outfh join("\t",
				  $exon,
				  $gene_hits{'exon'}{$exon}->[2],
				  $length,
				  'exon',
				  $frac),"\n";
	    }
	}
	close($outfh);
    }
}

exit;

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	$lane_files{$pair}{$lane}=1;
    }

    return(\%lane_files);
}

sub process_feature {
    my $feat=shift;
    my $chr=shift;
    my $sam=shift;
    my $gene_id=shift;

    # Make sure the features are in the correct order to avoid errors in the
    # estimation of minus strand features
    my ($rstart,$rend)=sort {$a <=> $b} ($feat->start(),
					 $feat->end());
    my @alignments = $sam->get_features_by_location(-seq_id => $chr,
						    -start  => $rstart,
						    -end    => $rend);
    
    my %hits;
    my $hits=0;
    my $feat_hits=[0,0,0];
    for my $al (@alignments) {	
	my $seqid  = $al->seq_id;
	
	# where does the alignment start in the reference sequence
	my $start  = $al->start;
	my $end    = $al->end;
	my $cigar  = $al->cigar_str;
	my $strand = $al->strand;

	# check mismatches
	my $mismatches=$al->get_tag_values('NM') || 0;
	if ($mismatches > 2) {
	    next;
	}

	my $multimaps = 0;
	my $hits0=$al->get_tag_values('H0') || 0;
	my $hits1=$al->get_tag_values('H1') || 0;
	my $hits2=$al->get_tag_values('H2') || 0;
	if ($hits0 > 1) {
	    $multimaps=1;
	} elsif ($hits1 > 1) {
	    $multimaps=1;
	} elsif ($hits2 > 1) {
	    $multimaps=1;
	}

	if ($multimaps) {
	    next;
	} elsif ($rstart > $start) {
	    # This removes the possibility of counting the junction
	    # reads twice if one end is in one interval and the other in
	    # another
	    # Skip reads that start outside of the feature
	    next;
	} elsif ($rend < $end) {
	    # Skip reads that end outside of the feature
	    next;
	} elsif ($cigar eq '*') {
	    # Skip hits with no cigar line
	    next;
	} elsif ($cigar=~/(\d+)M(\d+)[NI](\d+)M/o) {
	    # Skip junction hits
	    next;
	}
	# where does the alignment start in the query sequence
	my $query_start  = $al->query->start;
	my $query_end    = $al->query->end;
	my $query_length = $al->query->length; # matched length
		
	$hits{$strand}++;
	$hits++;
    }
    my $forward=$hits{1} || 0;
    my $reverse=$hits{-1} || 0;
    if ($hits) {
	$feat_hits=[$forward,$reverse,$hits];
    }
    return($feat_hits);
}
