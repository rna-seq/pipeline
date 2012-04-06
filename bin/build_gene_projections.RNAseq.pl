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
# This script should take annotation file and project all the exons in
# order to obtain a set of non-overlapping projections for each gene

use Bio::Range;
use RNAseq_pipeline3 qw(get_fh parse_gff_line get_log_fh);
use RNAseq_pipeline_settings3 ('read_config_file');
use Getopt::Long;

# Declare some variables
my $prefix;
my $exondir;
my $genomedir;
my $exonfile;
my $projectionfile;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$genomedir=$options{'GENOMEDIR'};

$exonfile=$exondir.'/'.$prefix.'.exon.gtf';
$projectionfile=$genomedir.'/'.$prefix.'.proj.gtf';

my $logfh=get_log_fh('build_gene_projections.log');

# Do some checks to see we have all necessary information
unless (-r $exonfile) {die "No exon file found at $exonfile\n";}

# If the file exists already skip creation
if (-r $projectionfile) {
    print $logfh $projectionfile,"\tis present, Skipping creation...\n";
} else {
    print $logfh "Building $projectionfile...\n";
    # First read in the exon file and get the list of exons belonging to
    # each gene
    my %genes;
    my %exons=%{get_exons($exonfile,
			  \%genes,
			  $logfh)};

    # Calculate the projection for each gene
    get_projections(\%genes,
		    $logfh);

    # Print out the results
    my $projfh=get_fh($projectionfile,1);
    foreach my $gene (keys %genes) {
	my $chr=shift(@{$genes{$gene}});
	foreach my $exon (@{$genes{$gene}}) {
	    my $strand='.';
	    if ($exon->strand == 1) {
		$strand='+';
	    } elsif ($exon->strand == -1) {
		$strand='-';
	    } else {
		warn "Unkonwn strand\n";
	    }
	    print $projfh join("\t",
			       $chr,
			       'Proj',
			       'exon',
			       $exon->start,
			       $exon->end,
			       '.',
			       $strand,
			       '.',
			       'gene_id "'.$gene.'"; transcript_id "'.$gene.'"'),"\n";
	}
    }
    close($projfh);
}
my $total=`wc -l $projectionfile`;
chomp($total);
$total=(split(/\s+/,$total))[0];
print join("\t",
	   'Proj',
	   $total),"\n";
close($logfh);

exit;

# TO DO use a get_coords from exon name sub for getting the exon coordinates
sub get_projections {
    my $genes=shift;
    my $logfh=shift;

    print $logfh 'Calculating gene projections...';
    foreach my $gene (keys %{$genes}) {
	my $chr;
	my @ranges;
	my $length=0;
	foreach my $exon (keys %{$genes->{$gene}}) {
	    my @location=split('_',$exon);
	    my $strand=pop(@location);
	    my $end=pop(@location);
	    my $start=pop(@location);
	    my $exon_chr=join('_',@location);
	    if ($chr &&
		$exon_chr ne $chr) {
		warn "$gene has exons in $chr and $exon_chr\n";
	    } else {
		$chr=$exon_chr;
	    }
	    my $range=Bio::Range->new(-start => $start,
				      -end => $end,
				      -strand => $strand);
	    push @ranges,$range;
	}
	my @disc_ranges = Bio::Range->disconnected_ranges(@ranges);
	$genes->{$gene}=[$chr,@disc_ranges];
    }
    print $logfh "done\n";
}

sub get_exons {
    my $exonfile=shift;
    my $genes=shift;
    my $logfh=shift;
    my $repeated_exons='rep.exons.txt';

    my %exons;
    my %genes;
    my %remove;

    print $logfh "Extracting exon lists from $exonfile...";
    my $exonfh=get_fh($exonfile);
    while (my $line=<$exonfh>) {
	my %line=%{parse_gff_line($line)};

	# Complain if they give us something that is not an exon
	unless ($line{'type'} eq 'exon') {
	    warn "Non exon line: $line\n";
	}

	# get the strand
	my $strand=$line{'strand'};
	if ($strand eq '+') {
	    $strand=1;
	} elsif ($strand eq '-') {
	    $strand= -1;
	} else {
	    # complain
	    warn "Unstranded exon: Strand $line{'strand'},$line\n";
	}

	# Get the exon_id
	my $exon_id=join('_',
			 $line{'chr'},
			 $line{'start'},
			 $line{'end'},
			 $strand);

	# Get the gene_id
	my $gene_id=$line{'feature'}{'gene_id'};

	# Check we have everything
	unless ($gene_id && $exon_id) {
	    warn "Problem parsing $line\n";
	}
	if ($exons{$exon_id} && 
	    ($exons{$exon_id} ne $gene_id)) {
	    $remove{$exon_id}=1;
	} else {
	    $exons{$exon_id}=$gene_id;
	    $genes->{$gene_id}->{$exon_id}=1;
	}
    }
    close($exonfh);
    print $logfh "done\n";

    my $repeatfh=get_fh($repeated_exons,1);
    # Remove those exons that map to multiple genes
    foreach my $exon (keys %remove) {
	print $repeatfh $exon,"\n";
	delete $exons{$exon};
    }
    close($repeatfh);

    my $count=keys %remove;
    print $logfh $count,"\tExons mapping to multiple genes removed\n";

    $count=keys %exons;
    print $logfh $count,"\tExons obtained\n";
    $count=keys %{$genes};
    print $logfh $count,"\tGenes\n";

    return(\%exons);
}
