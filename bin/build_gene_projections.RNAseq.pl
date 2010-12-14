#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script should take an eson annotation file and project the exons in
# order to obtain a set of non-overlapping projections for each gene

use Bio::Range;
use RNAseq_pipeline2 qw(get_fh parse_gff_line);
use RNAseq_pipeline_settings ('read_config_file');
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
# Do some checks to see we have all necessary information
unless (-r $exonfile) {die "No exon file found at $exonfile\n";}

# If the file exists already skip creation
if (-r $projectionfile) {
    print STDERR $projectionfile,"\tis present, Skipping creation...\n";
} else {
    print STDERR "Building $projectionfile...\n";
    # First read in the exon file and get the list of exons belonging to
    # each gene
    my %genes;
    my %exons=%{get_exons($exonfile,
			  \%genes)};

    # Calculate the projection for each gene
    get_projections(\%genes);

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
			       'gene_id "'.$gene.'"'),"\n";
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

exit;

sub get_projections {
    my $genes=shift;

    print STDERR 'Calculating gene projections...';
    foreach my $gene (keys %{$genes}) {
	my $chr;
	my @ranges;
	my $length=0;
	foreach my $exon (keys %{$genes->{$gene}}) {
	    my @location=split('_',$exon);
	    my $exon_chr=shift(@location);
	    if ($chr &&
		$exon_chr ne $chr) {
		warn "$gene has exons in $chr and $exon_chr\n";
	    } else {
		$chr=$exon_chr;
	    }
	    my $range=Bio::Range->new(-start => $location[0],
				      -end => $location[1],
				      -strand => $location[2]);
	    push @ranges,$range;
	}
	my @disc_ranges = Bio::Range->disconnected_ranges(@ranges);
	$genes->{$gene}=[$chr,@disc_ranges];
    }
    print STDERR "done\n";
}

sub get_exons {
    my $exonfile=shift;
    my $genes=shift;
    my $repeated_exons='rep.exons.txt';

    my %exons;
    my %genes;
    my %remove;

    print STDERR "Extracting exon lists from $exonfile...";
    my $repeatfh=get_fh($repeated_exons,1);
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
    close($repeatfh);
    print STDERR "done\n";

    # Remove those exons that map to multiple genes
    foreach my $exon (keys %remove) {
	delete $exons{$exon};
    }

    my $count=keys %remove;
    print STDERR $count,"\tExons mapping to multiple genes removed\n";

    $count=keys %exons;
    print STDERR $count,"\tExons obtained\n";
    $count=keys %{$genes};
    print STDERR $count,"\tGenes\n";

    return(\%exons);
}
