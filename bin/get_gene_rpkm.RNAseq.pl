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
# This script should take a file resulting form running the overlap script and
# extracting the total overlap reads and a file resulting from projecting all
# the genes onto the genome. It will calculate the RPKM for each of
# the genes in the file

use Bio::Range;
use Bio::SeqIO;
use RNAseq_pipeline3 qw(get_fh parse_gff_line get_feature_overlap);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list',
			       'get_gene_from_short_junc_sub','get_lanes',
			       'process_exons','get_gene_coverage_1000nt');

# Declare some variables
my $prefix;
my $exondir;
my $genomedir;
my $mapper;
my $exonfile;
my $projfile;
my $tmpdir;
my $stranded;
my $readlength;
my $junctiondir;
my $junctionsfasta;
my $min_read_length=76;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$genomedir=$options{'GENOMEDIR'};
$mapper=$options{'MAPPER'};
$exonfile=$exondir.'/'.$prefix.'.exon.gtf';
$projfile=$genomedir.'/'.$prefix.'.proj.gtf';
$tmpdir=$options{'LOCALDIR'};
$stranded=$options{'STRANDED'};
$readlength=$options{'READLENGTH'};
$junctiondir=$options{'JUNCTIONSDIR'};
$junctionsfasta=$options{'JUNCTIONSFASTA'};

# Connect to the database
my $dbh=get_dbh();
my $dbh_common=get_dbh(1);

# Get some usefull subs
if ($stranded) {
    *junc2gene=get_gene_from_short_junc_stranded_sub($dbh_common);
} else {
    *junc2gene=get_gene_from_short_junc_sub($dbh_common);
}

# Get the exon list;
my %genes;
my %exons=%{get_exon_list_from_gtf($exonfile,
				   \%genes)};

# Get the projected length and build a file with the gene projections
my %gene_lengths;
get_projected_length(\%genes,
		     \%gene_lengths,
		     $projfile,
		     $readlength);

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes()};

# Get unique maps for each lane
my $mappingtable=$prefix.'_genome_mapping';
my %unique_maps=%{get_unique_maps($dbh,
				  $mappingtable,
				  \%files,
				  $mapper)};

# Read and process the overlap files
foreach my $lane (keys %lanes) {
    my $type;
    if (keys %{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (keys %{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my %rpkms;

    # Get the exon overlap information
    my $exonoverlap=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.proj.overlap.total';

    if (-r $exonoverlap) {
	print STDERR "Processing $exonoverlap\n";
    } else {
	die "Can't read $exonoverlap\n";
    }

    # Get the junction information
    my @junctionfns;
    my %gene_juncs_coverage;

    foreach my $track (keys %{$lanes{$lane}}) {
	my $juncfilename=$junctiondir.'/'.$track.'.single.unique.overlap.total';
	push @junctionfns, $juncfilename;
    }

    # Process each of the junction files
    foreach my $juncfilename (@junctionfns) {
	if (-r $juncfilename) {
	    print STDERR "Processing $juncfilename\n";
	} else {
	    die "Can't read $juncfilename\n";
	}
	get_feature_coverage_junctions($juncfilename,
				       \%gene_juncs_coverage,
				       $readlength);
    }


    # Read in the exon coverage
    my %gene_coverage;
    get_gene_coverage_1000nt($exonoverlap,
			     \%gene_coverage,
			     $readlength,
			     $min_read_length);

    # For each of the exons in the exon list print the RPKM by normalizing the
    # reads per 1000 bases by the number of uniquely mapped reads
    my $outfile=$genomedir.'/'.$lane.'.'.$type.'.gene.rpkm.txt.gz';
    process_exons(\%genes,
		  \%gene_coverage,
		  \%gene_juncs_coverage,
		  $outfile,
		  $unique_maps{$lane},
		  $lane,
		  \%gene_lengths);
}

exit;

sub get_gene_from_short_junc_stranded_sub {
    my %options=%{read_config_file()};
    my $dbh=shift;
    my $table=$options{'JUNCTIONSTABLE'};

    # For saving time, as the junctions table is huge
    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT gene_id, exon1_id, exon2_id ';
    $query.="FROM $table ";
    $query.='WHERE chr = ? AND start = ? AND end = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $junc=shift;
	
	### TO DO fix for problematic chromosomes
	# This should be fixed in the start
	my ($chr,$start,$splice,$end,$strand)=split('_',$junc);

	print STDERR "genestrand= $strand\n";

	unless ($cache{$junc}) {
	    $count=$sth->execute($chr,$start,$end);

	    if ($count == 0) {
		die "No gene in $table corresponds to $junc\n";
	    } else {
		while (my ($gene,$exon1,$exon2)=$sth->fetchrow_array()) {
		    my $strand1=(split('_',$exon1))[-1];
		    my $strand2=(split('_',$exon2))[-1];
		    my $genestrand='+';
		    if ($strand1 eq $strand2) {
			if ($strand1 == -1) {
			    $genestrand='-';
			}
		    } else {
			die "Strand inconsistency between $exon1 and $exon2\n";
		    }
		    if ($genestrand eq $strand) {
			push @{$cache{$junc}}, $gene;
		    } else {
			warn "Problem with $junc,$strand,$exon1,$exon2,$genestrand\n";
		    }
		}
	    }
	}
	return($cache{$junc});
    };
    return($subroutine);
}

sub get_feature_coverage_junctions {
    my $infn=shift;
    my $features=shift;
    my $readlength=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage)=split("\t",$line);

	$coverage=$coverage * 1000;
	my $genes=junc2gene($feature);
	foreach my $gene (@{$genes}) {
	    $features->{$gene}+=$coverage;
	}
    }

    print STDERR "done\n";
}

# This will get all the exons from the annotation file
sub get_exon_list_from_gtf {
    my $exonfile=shift;
    my $genes=shift;
    my $repeated_exons='rep.exons.txt';

    my %exons;
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
	$gene_id=~s/"//g;

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

sub get_projected_length {
    my $genes=shift;
    my $lengths=shift;
    my $projfile=shift;
    my $read_length=shift;

    print STDERR "Calculating gene projections\n";

    # If the projected file exists read the projection lengths from it, and
    # else create it
    if (-r $projfile) {
	print STDERR $projfile,"\tpresent\n";
	my $projfh=get_fh($projfile);
	while (my $line=<$projfh>){
	    my %line=%{parse_gff_line($line)};
	    my $gene=$line{'feature'}{'gene_id'};
	    unless ($gene) {
		die "Unknown gene in $line\n";
	    }
	    $gene=~s/"//g;
	    # substract the read length
	    my $length=$line{'end'} - $line{'start'} + 1;
	    $lengths->{$gene}+=$length;
	}
	close($projfh);
    } else {
	print STDERR $projfile,"\tAbsent...building\n";
	my $projfh=get_fh($projfile,1);
	foreach my $gene (keys %{$genes}) {
	    my $chr;
	    my @ranges;
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
	    foreach my $range (@disc_ranges) {
		my $length=$range->length();
		my $strand;
		if ($range->strand() == 1) {
		    $strand='+';
		} elsif ($range->strand() == -1) {
		    $strand='-';
		} else {
		    warn "Unknown strand for $gene\n";
		}
		my $string='gene_id "'.$gene.'";';
		print $projfh join("\t",
				   $chr,
				   'RNAseqPipe',
				   'proj',
				   $range->start(),
				   $range->end(),
				   '.',
				   $strand,
				   '.',
				   $string),"\n";
		## substract the read length
		$lengths->{$gene}+=$length;
	    }
	}
	close($projfh);
    }
    print STDERR "done\n";
}

sub get_unique_maps {
    my $dbh=shift;
    my $maptable=shift;
    my $files=shift;
    my $mapper=shift;

    my %unique_maps;
    my %lane2pair;

    my ($query,$sth,$count);

    print STDERR "Getting unique mappings from $maptable...";
    foreach my $file (keys %{$files}) {
	$lane2pair{$files->{$file}->[1]}=$files->{$file}->[0];
    }

    $query ='SELECT LaneName, uniqueReads ';
    $query.="FROM $maptable";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count == keys %{$files}) {
	die "Wrong number of mapped reads\n";
    }

    while (my ($lane,$unique)=$sth->fetchrow_array()) {
	my $pair=$lane2pair{$lane};
	$unique_maps{$pair}+=$unique;
    }
    print STDERR "done\n";

    foreach my $lane (keys %unique_maps) {
	print STDERR join("\t",
			  $lane,
			  $unique_maps{$lane}),"\n";
    }

    return(\%unique_maps);
}
