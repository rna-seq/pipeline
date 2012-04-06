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



# Declare some variables
my $threshold=1;
my $tmpdir;
my $clusterdir;
my $readdir;
my $genomedir;
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
$genomedir=$options{'PROJECT'}.'/'.$options{'GENOMEDIR'};
$readdir=$options{'PROJECT'}.'/'.$options{'READDIR'};
$exondir=$options{'EXONDIR'};
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

my $annotation=$genomedir.'/'.$prefix.'.proj.gtf';
print STDERR $annotation,"\n";


# First get the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get for each lane the BAM files
my %lane_files=%{get_lane_files(\%files,
				$readdir)};

my %coverage;

# Extract the genes from the annotation
my %genes=%{get_annotation_from_gtf($annotation)};

foreach my $pair (keys %lane_files) {
    foreach my $lane (keys %{$lane_files{$pair}}) {
	my $bamfile=$lane_files{$pair}{$lane};
	print STDERR "$bamfile\n";
	my $outfile=$genomedir.'/'.$lane.'.single.unique.gtf.proj.overlap.total';
	
	if (-r $outfile) {
	    print STDERR $outfile,"\tis present. Skipping\n";
	    next;
	}
	
	# Get the necessary filehandles
	my $sam = Bio::DB::Sam->new(-bam  =>$bamfile,
				    -autoindex => 1);
	my $outfh=get_fh($outfile,1);
	
	# Process the genes
	foreach my $gene (keys %genes) {
	    my @exons=$genes{$gene}{'gene'}->exons();
	    my $chr=$genes{$gene}{'chr'};
	    my %gene_hits;
	    # Get the exon information
	    foreach my $exon (@exons) {
		my $feat_hits=process_feature($exon,
					      $chr,
					      $sam,
					      $gene);
		
		$gene_hits{'exon'}{$exon}=$feat_hits;
		$gene_hits{'gene'}+=$feat_hits->[2];
	    }
	    	    
	    # Print the output
	    foreach my $exon (@exons) {
		my $exon_id=join("\t",
				 $chr,
				 $exon->start(),
				 $exon->end(),
				 $exon->strand());
		# only print if there are actually reads in the gene
		if ($gene_hits{'exon'}{$exon}->[2]) {
		    my $length=$exon->end() - $exon->start() + 1;
		    my $frac=sprintf "%.2f",$gene_hits{'exon'}{$exon}->[2]/$length;
		    print $outfh join("\t",
				      $gene,
				      $gene_hits{'exon'}{$exon}->[2],
				      $length,
				      'exon',
				      $frac),"\n";
		}
	    }
	    
	}
	close($outfh);
    }
}

exit;

sub get_lane_files {
    my $files=shift;
    my $dir=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	$lane_files{$pair}{$lane}=$dir.'/'.$file;
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

	if ($rstart > $start) {
	    # This removes the possibility of counting the junction
	    # reads twice if one end is in one interval and the other in
	    # another
	    # Skip reads that start outside of the feature
	    next;
	} elsif ($cigar=~/(\d+)M(\d+)[NI](\d+)M/o) {
	    # Skip junction hits
	    next;
	} elsif ($rend < $end) {
	    # Skip reads that end outside of the feature
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
