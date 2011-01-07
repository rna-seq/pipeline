#!/soft/bin/perl

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
# This script will take as an input the overlap result of an experiment and it
# will extract from it the genes that have at least one RNA seq hit that
# overlaps completely with the gene 

# Load some modules
use RNAseq_pipeline3 qw(get_fh parse_gff_line);

# Declare some variables
my $threshold=1;

my $infile=shift;
my %genes;

# Check that we have the correct input files
unless($infile) {
    die "Input sequence file required\n";
}

# First read in the file
my %exons;
my $infh=get_fh($infile);
while (my $line=<$infh>) {
    my %line=%{parse_gff_line($line)};

    # Build the exon id
    my $strand;
    if ($line{'strand'} eq '+') {
	$strand=1;
    } elsif ($line{'strand'} eq '-') {
	$strand=-1;
    } else {
	print STDERR $line{'strand'},"\n";
	die "Unknown strand in $line\n";
    }

    my $gene_id=$line{'feature'}{'gene_id'};
    my $hits=$line{'feature'}{'list_feat2:'};

    my $exon_id=join('_',
		     $line{'chr'},
		     $line{'start'},
		     $line{'end'},
		     $strand);

    unless ($gene_id) {
	warn "No gene ID found in $line\n";	
    }
    $exons{$exon_id}{$gene_id}=1;

    # Check if the hit is within the target
    my @targets=split(',',$hits);
    my $start=$line{'start'};
    my $end=$line{'end'};
    my $length=$end - $start + 1;
    $genes{$exon_id}->[2]=$line{'type'};
    $genes{$exon_id}->[1]=$length;
    my @valid;
    if ($line{'feature'}{'nb_ov_feat2:'} > 0) {
	foreach my $hit (@targets) {
	    # This is because sometimes the scaffold has an undescore
	    my @parts=split('_',$hit);
	    my $h_strand=pop(@parts);
	    my $h_end=pop(@parts);
	    my $h_start=pop(@parts);
	    my $h_chr=join('_',@parts);
	    unless ($h_chr eq $line{'chr'}) {
		warn "Wrong chromosome\n";
	    }
	    if (($h_start >= $start) &&
		($h_end <= $end)) {
		push @valid, $hit;
	    }
	}

	if (@valid) {
	    $line{'feature'}{'nb_ov_feat2'}=@valid;
	    $line{'feature'}{'list_feat2:'}=join(',',@valid);
	    $genes{$exon_id}->[0]=@valid;
	} else {
	    $genes{$exon_id}->[0]=0;
	}
    } else {
	$genes{$exon_id}->[0]=0
    }
}
close($infh);
# print results
foreach my $exon_id (keys %genes) {
    unless ($genes{$exon_id}->[0]) {
	next;
    }
    foreach my $gene_id (keys %{$exons{$exon_id}}) {
	my $coverage=sprintf "%.2f",$genes{$exon_id}->[0]/$genes{$exon_id}->[1];
	print join("\t",
		   $gene_id,
		   @{$genes{$exon_id}},
		   $coverage),"\n";
    }
}

exit;
