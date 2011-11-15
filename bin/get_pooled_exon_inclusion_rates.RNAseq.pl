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
# This script should take all the files with the exclusion junctions generated
# from the annotation, and using the information from the database
# for the project it
# will extract those junctions that correstpond to inclusion junctions.
# After this it will determine how many reads fall in the inclusion junctions
# and the exon itself and how many fall in the exclusion reads, and it
# will return for each exon these two numbers together with a ratio between them
# It will provide the inclusion reads from the exon mappings separate from
# those of the junctions
# And it will also nomalize the values according to the number of reads in each
# lane.

use RNAseq_pipeline3 qw(get_fh get_exon_coverage_1000nt);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list',
			       'get_lanes','get_groups','get_unique_maps');
use Bio::SeqIO;

# Define some variables
my $project;
my $prefix;
my $exclusionfn;
my $genomedir;
my $junctiondir;
my $exondir;
my $mapper;
my $juncfile;
my $readlength;

# Read the configuration file
my %options=%{read_config_file()};
$exclusionfn=$options{'EXCLUSIONFILE'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$genomedir=$options{'GENOMEDIR'};
$junctiondir=$options{'JUNCTIONSDIR'};
$exondir=$options{'EXONDIR'};
$mapper=$options{'MAPPER'};
$juncfile=$options{'JUNCTIONSFASTA'};
$readlength=$options{'READLENGTH'};

# Check the input
unless ($exclusionfn) {
    die "No exclusion file???\n";
}

# Connect to the database
my $dbh=get_dbh();
my $dbhcommon=get_dbh(1);

# Get some subroutines
my $table=get_junctions_table($dbh,
			      $prefix);
my $mappingtable=$prefix.'_unique_maps_genome';
*get_inclusion_exons=get_inclusion_junctions_sub($dbhcommon,
						 $table);

my %files=%{read_file_list()};
my %lanes=%{get_lanes()};
my %groups=%{get_groups()};

# Get normalization ratios for each lane
my %read_no_norm=%{get_unique_maps($dbh,
				   $mappingtable,
				   \%files,
				   \%groups,
				   $mapper)};

foreach my $pair (keys %read_no_norm) {
    my $read_norm=sprintf "%.3f",$read_no_norm{$pair};
    print STDERR join("\t",
		      $pair,
		      $read_norm),"\n";
}

# Get junction lengths
my %junclengths=%{get_junction_length($juncfile,
				      $readlength)};

# Read and process the files
foreach my $group (keys %groups) {
    my %exons_coverage;
    my %juncs_coverage;
    my $normfactor=1000000 / $read_no_norm{$group};
    foreach my $lane (keys %{$groups{$group}}) {
	my $type;
	if (keys %{$lanes{$lane}} == 1) {
	    $type='single';
	} elsif (keys %{$lanes{$lane}} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type\n";
	}
	my %genes;
	
	# Get the exon information
	my $exonfilename=$exondir.'/'.$lane.'.'.$type.'.unique.gtf.overlap.total';
	
	if (-r $exonfilename) {
	    print STDERR "Processing $exonfilename\n";
	} else {
	die "Can't read $exonfilename\n";
	}
	
	# Read in the exon coverage
	get_exon_coverage_1000nt($exonfilename,
				 \%exons_coverage,
				 $readlength);
	
	# Get the junction information
	my @junctionfns;
	
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
					   \%juncs_coverage,
					   \%junclengths);
	}
    }
    # Go through the exclusion file and for each of the exons print the
    # Gene, exon, number of inclusion reads and number of exclusion reads
    my $outfile=$exondir.'/'.$group.'.inclusion.exclusion.pooled.txt';
    print STDERR $normfactor,"\n";
    process_exons($exclusionfn,
		  \%juncs_coverage,
		  \%exons_coverage,
		  $outfile,
		  $normfactor,
		  $group);
}

exit;

sub get_junction_length {
    my $juncfile=shift;
    my $read_length=shift;

    my %lengths;

    print STDERR 'Getting junction lengths...';
    my $infh=Bio::SeqIO->new(-file => $juncfile,
			     -format => 'fasta');

    while (my $seq=$infh->next_seq()) {
	my $juncid=$seq->display_id();
	my $length=$seq->length();

	if ($length > (2 * $read_length)) {
	    $length=2 * $read_length;
	}
	
	my $junc_short=get_junc_id_from_long_id($juncid);

	if ($lengths{$junc_short}) {
	    if ($lengths{$junc_short} < $length) {
		$lengths{$junc_short}=$length;
	    }
	} else {
	    $lengths{$junc_short}=$length;
	}

	$lengths{$juncid}=$length;
    }
    print STDERR "done\n";

    return(\%lengths);
}


### TO DO
# This sub should be improved for speed using the MySQL DB or a cache
sub get_junc_id_from_long_id {
    my $junc=shift;

    my ($ex1,$ex2)=split('_splice_',$junc);
    my @ex1=split('_',$ex1);
    my @ex2=split('_',$ex2);

    pop(@ex1);
    my $start=pop(@ex1);
    pop(@ex1);
    my $chr=join('_',@ex1);

    pop(@ex2);pop(@ex2);
    my $end=pop(@ex2);

    my $short_id=join('_',$chr,$start,'splice',$end);

    return($short_id);
}

sub get_junctions_table {
    my $dbh=shift;
    my $prefix=shift;

    my $table=$prefix.'_junctions';

    my ($query,$sth,$count);
    $query ='SELECT table_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count == 1) {
	warn "Incorrect number of records retrieved from $table\n";
    }

    my ($junctable)=$sth->fetchrow_array();

    return($junctable);
}
    
sub get_inclusion_junctions_sub {
    my $dbh=shift;
    my $table=shift;

    my ($query,$sth);
    $query ='SELECT DISTINCT chr, start, end ';
    $query.="FROM $table ";
    $query.='WHERE exon1_id = ? OR exon2_id = ?';
    $sth=$dbh->prepare($query);

    my $no_hits=0;
    my $processed=0;
    my $get_exons=sub {
	my $exon=shift;
	my @junctions;
	my $count=$sth->execute($exon,$exon);

	while (my ($chr,$start,$end)=$sth->fetchrow_array()) {
	    my $junction=join('_',$chr,$start,'splice',$end);
	    push @junctions, $junction;
	}
	return(@junctions);
    };
    return($get_exons);
}

sub process_exons {
    my $filename=shift;
    my $junctions=shift;
    my $exons=shift;
    my $output=shift;
    my $normfactor=shift || 1;
    my $lane=shift;

    print STDERR "Collecting events...\n";
    print STDERR $normfactor,"\n";

    my $infh=get_fh($filename);
    my $outfh=get_fh($output,1);

    while (my $line=<$infh>) {
	my ($inc_rate_ex,$inc_rate_junc,$ex_rate)=(0,0,0);
	chomp($line);
	my ($gene,$exon_id,$exclusion)=split("\t",$line);
	my @inclusion_exons=get_inclusion_exons($exon_id);
	my @exclusion_exons=split(',',$exclusion);

	if ($exons->{$exon_id}) {
	    $inc_rate_ex+=$exons->{$exon_id};
	}

	foreach my $inc (@inclusion_exons) {
	    if ($junctions->{$inc}) {
		$inc_rate_junc+=$junctions->{$inc};
	    }
	}

	foreach my $exc (@exclusion_exons) {
	    if ($junctions->{$exc}) {
		$ex_rate+=$junctions->{$exc};
	    }
	}

	# Normalize by read number
	$inc_rate_ex=sprintf "%.3f",$inc_rate_ex * $normfactor;
	$inc_rate_junc=sprintf "%.3f",$inc_rate_junc * $normfactor;
	$ex_rate=sprintf "%.3f",$ex_rate * $normfactor;


	my $inclusion='\N';
	my $all_inclusion=$inc_rate_junc + $inc_rate_ex;
	my $all_events=$inc_rate_junc + $inc_rate_ex + $ex_rate;
	if ($all_events > 0) {
	    $inclusion= sprintf "%.3f",$all_inclusion / $all_events;
	}

	# Print results
	print $outfh join("\t",
			  $gene,
			  $exon_id,
			  $inc_rate_ex,
			  $inc_rate_junc,
			  $ex_rate,
			  $inclusion,
			  $lane),"\n";
    }
    close($outfh);
    close($infh);

    print STDERR "\rDone\n";

}

sub get_feature_coverage_junctions {
    my $infn=shift;
    my $features=shift;
    my $lengths=shift;

    print STDERR "Getting mappings from $infn...\n";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage)=split("\t",$line);
	# Normalize to reads per 1000 nt
	my $junclength=$lengths->{$feature};
	unless ($junclength) {
	    # Skip those junctions that are not in the annotation and do not
	    # have information on their length
	    warn "No length for $feature\n";
	    next;
	}
	$coverage=$coverage * (1000 /$junclength);
	$features->{$feature}+=$coverage;
    }

    print STDERR "done\n";
}
