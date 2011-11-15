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
# This script should take a g[tf]f file and plot the number of mappings in known
# and in novel junctions in each of the files supplied

use RNAseq_pipeline3 qw(get_fh get_log_fh);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_dbh','get_junction_type_sub');

my $species;
my $project;
my $outfile;
my $table;
my $prefix;
my $debug=0;

my %options=%{read_config_file()};
$table=$options{'JUNCTIONSTABLE'};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};

# Get a log file
my $log_fh=get_log_fh('build_read_stats.RNAseq.log',
		      $debug);

# Get subroutines
*junction_type=get_junction_type_sub($table);

# Set the outfile
$outfile=$options{'PREFIX'}.'_unique_maps_junctions';

my %distribution;
unless ($outfile) {
    die "No output file name\n";
}

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};

my $junc_classfn=$prefix.'_junction_maps_class.txt';
my $junc_classfh=get_fh($junc_classfn,1);
foreach my $lane (keys %lanes) {
    my $type='single';

    my %classification;
    my $infilename=$lane.'.'.$type.'.unique.gtf.gz';
    print $log_fh "Processing $infilename\n";

    unless (-r $infilename) {
	die "$infilename is not readable\n";
    }

    build_distribution($infilename,
		       \%distribution,
		       \%classification,
		       $lane);

    foreach my $junction (keys %classification) {
	my $type=junction_type($junction);
	print $junc_classfh join("\t",
				 $junction,
				 $type,
				 $classification{$junction},
				 $lane),"\n";
    }

    print $log_fh "done\n";
}
close($junc_classfh);

build_graph($outfile,
	    $species,
	    $project,
	    \%distribution);

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[1]}++;
    }

    return(\%lanes);
}

sub build_distribution {
    my $file=shift;
    my $dist=shift;
    my $class=shift;
    my $lane=shift;
    my $infh=get_fh($file);

    while (my $line1=<$infh>) {
	chomp($line1);
	my @line1=split("\t",$line1);
	my $line2=<$infh>;
	chomp($line2);
	my @line2=split("\t",$line2);
	my $junction_id=join('_',@line1[0,4],$line2[3]);
	my $type=junction_type($junction_id);
	unless($type) {
	    die "$junction_id has no type\n";
	}
	$class->{$junction_id}++;
	$dist->{$lane}->{$type}++;
    }

    close($infh);
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
	foreach my $chr (sort keys %{$dist->{$filename}}) {
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
    $r_string.="stats1<-read.table(\"$stats_file\",sep=\"\t\")\n";
    $r_string.="stats<-tapply(stats1\$V3,stats1[,1:2],sum)\n";
    $r_string.="postscript(\"$graph.ps\")\n";
    $r_string.='barplot(stats,names.arg=attributes(stats)$row.names,beside=T,';
    $r_string.="main=\"$species $project junctions unique mappings\",";
    $r_string.="xlab=\"Junction Type\",ylab=\"Unique reads\",legend=T,";
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
