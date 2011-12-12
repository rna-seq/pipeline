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
# This script should take a gtf/gff containing exons with
# gene names and it will build a set of all exon sequences

use RNAseq_pipeline3 qw(get_fh get_log_fh parse_gff_line get_annotation_from_gtf);
use RNAseq_pipeline_settings3 qw(read_config_file);
use Bio::SeqIO;
use Bio::DB::Fasta;

my $annotation_file;
my $genome_file;
my $debug=0;

my %options=%{read_config_file()};
$annotation_file=$options{'ANNOTATION'};
$genome_file=$options{'GENOMESEQ'};
my $exonfn=$options{'EXONSFASTA'};

# Get a log file
my $log_fh=get_log_fh('build_exon_seqs.RNAseq.log',
		      $debug);

print $log_fh "Extracting the exon sequences for $annotation_file\n";

unless ($annotation_file && $genome_file) {
    die "ERROR:An annotation file and a genome file are required\n";
}

# Check if the file exists already
my $present=check_file($exonfn);

if ($present) {
    print $log_fh "$exonfn is present. Skipping\n";
    exit;
}

print $log_fh "Processing $annotation_file\n";

# Extract all the exons coordinates belonging to each of the genes
my %exons=%{get_annotation_from_gtf($annotation_file,
				    $log_fh,
				    'exons')};

# Get the sequence of the chromosome as well as an accessor to extract the exon
# sequence
#my %chromosomes=%{get_sequences($genome_file)};
*get_exon_seq=get_exon_seq_sub($genome_file);

# Get all the exon sequences
my %exon_seqs=%{build_exon_sequences(\%exons)};

# Print out all exon sequences
print $log_fh "Printing exon sequences in $exonfn\n";
my $exonfh=get_fh($exonfn,1);
foreach my $exon (keys %exon_seqs) {
    my $seq=$exon_seqs{$exon};
    print $exonfh ">$exon\n";
    for (my $pos=0;$pos < length($seq) ; $pos += 60) {
	print $exonfh substr($seq,$pos,60),"\n";
    }
}
close($exonfh);
print $log_fh "done\n";

exit;

sub check_file {
    my $file=shift;

    my $present=0;

    my $tablefh=get_fh($options{'PREFIX'}.'_exon_seqs.txt',1);
    if (-r $file) {
	$present=1;
	# Print the table location for the junctions of this experiment
	print $tablefh join("\t",
			    $file,
			    "Present"),"\n";
    } else {
	# Continue, as the table must be created
	print STDERR $file,"\tIs not present\n";
	print $tablefh join("\t",
			    $file,
			    "Generated"),"\n";
    }
    close ($tablefh);

    return($present);
}

# This should accelerate the process of obtaining the exon sequences at the
# expense of RAM
sub get_exon_seq_sub {
    my $genomefile=shift;

    # Generate the chromosomes
    print STDERR "Generating $genomefile index...";
    my $db=Bio::DB::Fasta->new($genomefile);
    print STDERR "done\n";

    my %exons;

    my %revcom=('a' => 't',
		'c' => 'g',
		't' => 'a',
		'g' => 'c',
		'A' => 'T',
		'C' => 'G',
		'T' => 'A',
		'G' => 'C');

    # Build the subroutine for the search
    my $get_seq_from_chrom= sub {
	my $exon_id=shift;
	
	unless (exists $exons{$exon_id}) {
	    my $seq;
	    my @line=split('_',$exon_id);
	    my ($chr,$start,$end,$strand);
	    if (@line == 4) {
		($chr,$start,$end,$strand)=@line;
	    } else {
		my @chr=splice(@line,0,@line - 3,());
		$chr=join('_',@chr);
		($start,$end,$strand)=@line;
	    }

	    if ($strand == 1) {
		$seq=$db->subseq($chr,$start,$end);
	    } elsif ($strand == -1) {
		$seq=$db->subseq($chr,$end,$start)
	    } else {
		warn "unknown strand $strand\n";
		return();
	    }

	    # Attempt to retrieve the sequence removing the initial chr
	    # If the initial multifasta for the genome did not have these
	    # tags the script will not find them if not.
	    my $chr2=$chr;
	    unless ($seq) {
#		$chr2=~s/^chr//;
		$seq=$db->subseq($chr2,$end,$start);
	    }
	    # Try again with the mitochondrial sequences
	    unless ($seq) {
#		$chr2=~s/M$/MT/;
		$seq=$db->subseq($chr2,$end,$start);
	    }

	    unless ($seq) {
		warn "No sequence retrieved for $exon_id... skipping\n";
		return();
	    }

	    # There is a problem when extracting sequences of length 1 as the
	    # subseq method cannot determine if they are plus or minus, and
	    # decides they are plus
	    # To fix this any case where the sequence is on base and strand is
	    # minus has to be reverse complemented

	    if (($start == $end) &&
		($strand == -1)) {
		if ($revcom{$seq}) {
		    $seq=$revcom{$seq};
		} else {
		    warn "Unknown sequence $seq\n";
		}
	    }

	    $exons{$exon_id}=$seq;
	}
	return($exons{$exon_id});
    };
    return($get_seq_from_chrom);
}

# Get the exon sequences from the coordinates
sub build_exon_sequences {
    my $exons=shift;

    my %exon_seqs;
    my $count=0;

    print STDERR "Building all exon sequences\n";
    foreach my $exon (keys %{$exons}) {
	my $exon_seq=get_exon_seq($exon);
	if ($exon_seq) {
	    $exon_seqs{$exon}=$exon_seq;
	    $count++;
	}
	unless ($count % 1000) {
	    print STDERR "$count\tRetrieved\r";
	}
    }
    if ($count) {
	print STDERR $count,"\tAnnotated exons retrieved\n";
    } else {
	die "Sorry, I couldn't extract any exons...\n";
    }

    return(\%exon_seqs);
}
