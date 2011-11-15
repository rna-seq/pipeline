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
# This script will take the output of the flux capacitor and it will create a
# table to load in the database with this information

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh parse_gff_line);
use RNAseq_pipeline_settings3 qw(read_file_list get_gene_from_trans_sub);

# get subroutines
*trans2gene=get_gene_from_trans_sub();

# First get a list of the bed files we are going to process
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};
my %groups=%{get_groups(\%files)};

# read the results from the files named after the lane ids
foreach my $group (keys %groups) {
    my $filename=$group.'.flux.gtf';

    if (-r $filename) {
	process_file($filename,
		     $group);
	my $command="rm $filename";
	print STDERR "Executing:\t$command\n";
	system($command);
    } else {
	die "I can't find $filename\n";
    }
}

exit;

sub process_file {
    my $infile=shift;
    my $lane=shift;

    my $infh=get_fh($infile);

    while (my $line=<$infh>) {
	my %line=%{parse_gff_line($line)};

	my $type=$line{'type'};

	if ($type=~/transcript/) {
	    my $flux_id=$line{'feature'}{'gene_id'};
	    my $trans_id=$line{'feature'}{'transcript_id'};
	    my $rpkm=$line{'feature'}{'RPKM'};
	    my $gene_id=trans2gene($trans_id);
	    if ($gene_id &&
		($rpkm > 0)) {
		print join("\t",
			   $gene_id,
			   $trans_id,
			   $flux_id,
			   $rpkm,
			   $lane),"\n";
	    }
	}
    }
    close($infh);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	push @{$lanes{$files->{$file}->[0]}},$files->{$file}->[1];
    }

    return(\%lanes);
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
