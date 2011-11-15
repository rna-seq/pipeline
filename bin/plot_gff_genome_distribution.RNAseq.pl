#!/soft/bin/perl
# DGK

#    GRAPE
#    Copyright (C) 2009-2011 Centre for Genomic Regulation (CRG)
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
# This script should take the gtf file resulting from the mapping and selection
# of the unique reads and plot the resulting summaries using R
# These summaries will be the distribution along the chromosomes of the reads

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

my $species;
my $project;
my $paired;
my $genomedir;

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$genomedir=$options{'GENOMEDIR'};

# Get the outfile
my $outfile=$options{'PREFIX'}.'_unique_maps_genome';

my %distribution;
unless ($outfile) {
    die "No output file name\n";
}

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};

foreach my $lane (keys %lanes) {
    my $type;
    if ($lanes{$lane} == 1) {
	$type='single';
    } elsif ($lanes{$lane} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my $infilename=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.gz';
    print STDERR "Processing $infilename\n";

    unless (-r $infilename) {
	die "$infilename is not readable\n";
    }

    build_distribution($infilename,
		       \%distribution,
		       $lane);

    print STDERR "done\n";
    
}

build_graph($outfile,
	    $species,
	    $project,
	    \%distribution);

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}++;
    }

    return(\%lanes);
}

sub build_distribution {
    my $file=shift;
    my $dist=shift;
    my $lane=shift;

    my $infh=get_fh($file);

    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);
	$dist->{$lane}->{$line[0]}++;
    }

    close($infh);
}

# Because $a and $b belong to the main namespace (I think) this sub cannot
# be imported easily (or I couldn't) from a library.
sub sort_by_chromosome {
    my @chra=split('_',$a);
    my @chrb=split('_',$b);

    $chra[0]=~s/^chr//i;
    $chrb[0]=~s/^chr//i;

    my $ret;
    if (($chra[0]=~/^\d+$/) && ($chrb[0]=~/^\d+$/)) {
	$ret=$chra[0] <=> $chrb[0];
    } elsif ($chra[0]=~/^\d+$/) {
	$ret=-1;
    } elsif ($chrb[0]=~/^\d+$/) {
	$ret=1;
    }

    if ($ret) {
	return($ret);
    } else {
	return($a cmp $b);
    }
}

sub build_graph {
    my $graph=shift;
    my $species=shift;
    my $project=shift;
    my $dist=shift;

    my $stats_file=$graph;
    unless ($stats_file=~/\.txt$/) {
	$stats_file.='.txt';
    }
    my $tmpfh=get_fh($stats_file,1);
    foreach my $filename (keys %{$dist}) {
	foreach my $chr (sort sort_by_chromosome keys %{$dist->{$filename}}) {
	    my $chr_id=$chr;
	    $chr_id=~s/^chr//;
	    print $tmpfh join("\t",
			      $filename,
			      $chr_id,
			      $dist->{$filename}->{$chr}),"\n";
	}
    }
    # Close the file to make sure buffer is flushed
    close($tmpfh);

    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;
    my $lanes=keys %{$dist};
    $r_string ="cols=rainbow($lanes)\n";
    $r_string.="stats1<-read.table(\"$graph.txt\",sep=\"\t\")\n";
    $r_string.="stats<-tapply(stats1\$V3,stats1[,1:2],sum)\n";
    $r_string.="postscript(\"$graph.ps\")\n";
    $r_string.='barplot(stats[,order(as.numeric(colnames(stats)))],names.arg=attributes(stats)$row.names,beside=T,';
    $r_string.="main=\"$species $project unique mappings distribution\",";
    $r_string.="xlab=\"Chromosomes\",ylab=\"Unique reads\",legend=T,";
    $r_string.="col=cols)\n";
    $r_string.="dev.off()\n";

    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet < $execution_file";
    system($command);
    $command="rm $execution_file";
    system($command);
}    
