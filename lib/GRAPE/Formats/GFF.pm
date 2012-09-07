package GRAPE::Formats::GFF;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export subroutines to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
push @EXPORT_OK,('check_gff_file','parse_gff_line');

use strict;
use warnings;

# Load base modules
use GRAPE::Base ('get_fh');

# Parse strictly and return and die if any problems are encountered 
sub parse_gff_line {
    my $line=shift;
    my $log_fh=shift;

    # If no log_fh is provided redirect to STDERR
    unless ($log_fh) {
	$log_fh=*STDERR;
    }

    chomp($line);
    my %line;

    my @line=split("\t",$line);

    unless (@line >=8) {
	die "Insufficient fields in $line\n";
    }

    $line{'chr'}=$line[0];
    $line{'start'}=$line[3];
    $line{'end'}=$line[4];
    $line{'strand'}=$line[6];
    $line{'frame'}=$line[7];

    $line{'source'}=$line[1];
    $line{'type'}=$line[2];
    $line{'score'}=$line[5];
    
    # Parse the info fields
    if ($line[8]) {
	$line[8]=~s/^\s*//o;
	my @info;

	if ($line[8]=~/;/o) {
	    $line[8]=~s/;$//o;
	    @info=split(/; /,$line[8]);
	} elsif ($line[8]=~/^gene_id/o) {
	    # Check is there is only one key-value pair as in this case there
	    # is not necessarily a ';' If this is the case the attribute should
	    # be gene_id
		@info=($line[8]);
	} else {
	    die "Info field does not seem correct in:\n\t$line\n";
	}


	# If the info field is present print it
	if (@info) {
	    for (my $i=0;$i<@info;$i+=1) {
		my ($key,$value)=split(' ',$info[$i],2);

		if ($value) {
		    $value=~s/"//og;
		    # Add something in order to be compatible with the old
		    # version of overlap

		    # Allow for multiple values of the same tag
		    if ((exists $line{'feature'})&&
			(exists $line{'feature'}{$key})) {
			if ($key eq 'gene_id') {
			    die "gene_id is present more than once in $line\n";
			} elsif ($key eq 'transcript_id') {
			    die "transcript_id is present more than once in $line\n";
			} else {
			    $line{'feature'}{$key}.=",$value";
			}
		    } else {
			$line{'feature'}{$key}=$value;
		    }
		} else {
		    if ($line{'type'}!~/(exon|transcript|gene)/o) {
			print $log_fh "Unknown key $key in $line{'type'} line $line\n";
		    } else {
			die "No value found for key $key in $line{'type'} line $line\n";
		    }
		}
	    }
	} else {
	    die "Info field is missing from gtf in:\n$line\n";
	}
    }

    return(\%line);
}

# Check if a file conformas to the UCSC specifications of GTF
sub check_gff_file {
    my $file=shift;

    print STDERR "Checking annotation file...";
    my %features;
    my $infh=get_fh($file);
    while (my $line=<$infh>) {
	if ($line=~/^#/o) {
	    next;
	}
	my %line=%{parse_gff_line($line)};
	$features{$line{'type'}}++;
    }

    # Issue some warnings if the file looks strange
    if ($features{'exon'} < 1) {
	die "No exons found in $infh\n";
    } elsif ($features{'exon'} &&
	     $features{'exon'} < 10000) {
	print STDERR "WARNING: Only $features{'exon'} exons found in $file, maybe it is incomplete\n";
    }
    if ($features{'gene'} &&
	$features{'gene'} < 1) {
	print STDERR "WARNING: No genes found in $file\n";
    }
    if ($features{'transcript'} &&
	$features{'transcript'} < 1) {
	print STDERR "WARNING: No transcripts found in $file\n";
    }
    print STDERR "Annotation seems fine\n";

    close($infh);
}

1;
