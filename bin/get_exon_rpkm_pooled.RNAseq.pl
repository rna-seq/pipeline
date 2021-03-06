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
# This script should take a bam file and calculate from the aligned reads
# the rpkm at the exon level

use Bio::Range;
use RNAseq_pipeline3 qw(get_fh parse_gff_line get_exon_coverage_1000nt);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list',
			       'get_unique_maps');

# Declare some variables
my $prefix;
my $exondir;
my $exontable;
my $mapper;
my $readlength;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$exontable=$options{'EXONSCLASSTABLE'};
$mapper=$options{'MAPPER'};
$readlength=$options{'READLENGTH'};

# Connect to the database
my $dbh=get_dbh();
my $commondbh=get_dbh(1);

# Get the exon list;
my %exons=%{get_exon_list($commondbh,
			  $exontable)};

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};
my %groups=%{get_groups(\%files)};

# Get unique maps for each lane
my $mappingtable=$prefix.'_unique_maps_genome';
my %unique_maps=%{get_unique_maps($dbh,
				  $mappingtable,
				  \%files,
				  $mapper)};

# Read and process the overlap files
foreach my $group (keys %groups) {
    my %exons_coverage;
    foreach my $lane (@{$groups{$group}}) {
	my $type;
	if (keys %{$lanes{$lane}} == 1) {
	    $type='single';
	} elsif (keys %{$lanes{$lane}} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type\n";
	}

	# Get the exon overlap information
	my $exonoverlap=$exondir.'/'.$lane.'.'.$type.'.unique.gtf.overlap.total';

	if (-r $exonoverlap) {
	    print STDERR "Processing $exonoverlap\n";
	} else {
	    die "Can't read $exonoverlap\n";
	}

	# Read in the exon coverage
	get_exon_coverage_1000nt($exonoverlap,
				 \%exons_coverage,
				 $readlength);
    }
    # For each of the exons in the exon list print the RPKM by normalizing the
    # reads per 1000 bases by the number of uniquely mapped reads
    my $outfile=$exondir.'/'.$group.'.exon.rpkm.pooled.txt.gz';
    process_exons(\%exons,
		  \%exons_coverage,
		  $outfile,
		  $unique_maps{$group},
		  $group);
}

exit;

sub get_exon_list {
    my $dbh=shift;
    my $table=shift;

    my %exons;

    print STDERR "Retrieving exon list from $table...\n";
    my ($query,$sth,$count);
    $query ='SELECT DISTINCT exon_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count) {
	warn "No records retrieved from $table\n";
    }

    while (my ($exon_id)=$sth->fetchrow_array()) {
	$exons{$exon_id}=1;
    }
    print STDERR $count,"\tExons retrieved\n";

    return(\%exons);
}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}


sub process_exons {
    my $exonlist=shift;
    my $exoncoverage=shift;
    my $output=shift;
    my $uniquemaps=shift;
    my $lane=shift;

    print STDERR "Collecting events...\n";

    my $outfh=get_fh($output,1);

    foreach my $exon_id (keys %{$exonlist}) {
	my $rpkm=0;

	if ($exoncoverage->{$exon_id}) {
	    $rpkm+=$exoncoverage->{$exon_id};
	}

	# Normalize by read number
	$rpkm=sprintf "%.3f",($rpkm * 1000000) / $uniquemaps;

	# Print results only if they are positive
	if ($rpkm != 0) { 
	    print $outfh join("\t",
			      $exon_id,
			      $rpkm,
			      $lane),"\n";
	}
    }
    close($outfh);

    print STDERR "\rDone\n";
}

sub get_groups {
    my $files=shift;
    my %groups;
    
    foreach my $file (keys %{$files}) {
	my $group=$files->{$file}->[2] || 'All';
	push @{$groups{$group}},$files->{$file}->[0];
    }

    return(\%groups);
}
