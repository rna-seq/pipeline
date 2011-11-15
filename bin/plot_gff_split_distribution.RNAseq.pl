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
# This script should take a g[tf]f file containing split-mapping information
# Which means there is a single ID for the two halves of the match, 
# and plot the number of mappings in the same chromosome and strand, in
# same chromosome different strand and in different chromosomes

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_dbh);

my $species;
my $project;
my $outfile;

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};

# Set the outfile
$outfile=$options{'PREFIX'}.'_unique_maps_split';

my %distribution;

unless ($outfile) {
    die "No output file name\n";
}

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};


foreach my $lane (keys %lanes) {
    my $type='single';

    my $infilename=$lane.'.'.$type.'.unique.gtf.gz';
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
	$lanes{$files->{$file}->[1]}++;
    }

    return(\%lanes);
}

sub build_distribution {
    my $file=shift;
    my $dist=shift;
    my $lane=shift;
    my $infh=get_fh($file);
    my %splitmaps;

    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);
	my @info=split(/\s+/,$line[8]);
	my $split_id;
	for (my $i=0;$i<@info;$i++) {
	    if ($info[$i] eq 'ID') {
		$split_id=$info[$i + 1];
		last;
	    }
	}
	unless ($split_id) {
	    warn "No ID for $line\n";
	}
	push @{$splitmaps{$split_id}},[@line[0,3,4,6]]
    }
    close($infh);

    foreach my $junction (keys %splitmaps) {
	my @frags=@{$splitmaps{$junction}};
	unless (@frags == 2) {
	    warn "Wrong number of fragments for $junction\n";
	}

	# Count the occurrences
	if ($frags[0][0] eq $frags[1][0]){
	    if ($frags[0][3] eq $frags[1][3]) {
		# determine if they are in the right orientation or not
		if ($frags[0][3] eq '+') {
		    if ($frags[0][1] > $frags[1][2]) {
			$dist->{$lane}->{'SameChr&StrandRev'}++;
		    } elsif ($frags[0][2] > $frags[1][1]) {
			$dist->{$lane}->{'SameChr&StrandOver'}++;
		    } else {
			$dist->{$lane}->{'SameChr&StrandOK'}++;
		    }
		} elsif ($frags[0][3] eq '-') {
		    if ($frags[0][2] < $frags[1][1]) {
			$dist->{$lane}->{'SameChr&StrandRev'}++;
		    } elsif ($frags[0][1] < $frags[1][2]) {
			$dist->{$lane}->{'SameChr&StrandOver'}++;
		    } else {
			$dist->{$lane}->{'SameChr&StrandOK'}++;
		    }
		} else {
		    warn "Problem\n";
		}
	    } else {
		$dist->{$lane}->{'SameChrDifStrand'}++;
	    }
	} else {
	    $dist->{$lane}->{'DifChrom'}++;
	}
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
	foreach my $type (sort keys %{$dist->{$filename}}) {
	    print $tmpfh join("\t",
			      $filename,
			      $type,
			      $dist->{$filename}->{$type}),"\n";
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
