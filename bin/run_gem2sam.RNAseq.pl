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
# This script should take the gem-2-sam program and run it in the pipeline.
# As the program is in test phase currently, it will generate the header first
# also rename all the reads so paired reads have the same name.
# Finally after running it it will check the output to remove negative cases
# from the cigar string

# First of all it will get the files that we need to map. For each lane this
# will be the genome files and the split mapping files (we will leave the
# split mapping ones for the moment, as the bug is not necessarily fixed yet

# First read in the files we want to compare and write them to a temporary
# location after filtering

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub run_system_command);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list');
use Bio::SeqIO;
use Tools::Bam ('generate_bam_index','generate_sorted_bam');

my $genomefile;
my $prefix;
my $tmpdir;
my $samdir;
my $file_list;
my $bindir;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$tmpdir=$options{'LOCALDIR'};
$genomefile=$options{'GENOMESEQ'};
$samdir=$options{'SAMDIR'};
$file_list=$options{'FILELIST'};
$bindir=$options{'BIN'};

# Add the bindir to the path
$ENV{'PATH'}.=":$bindir";

# Connect to the database
my $dbh=get_dbh();

# First ge the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get some subroutines
my $merged_map_table=$options{'PREFIX'}.'_merged_mapping';
*get_merged_file=get_files_from_table_sub($dbh,
					  $merged_map_table);

# Get for each lane the files corresponding to the genome mapping and the
# junctions mapping
my %lane_files=%{get_lane_files(\%files)};


# Process the files
foreach my $pair (keys %lane_files) {
    # Add a check to make sure the file is actually present and die if not
    

    my $samfn=$tmpdir.'/'.$pair.'.merged.sam';
    my $bamfn=$samdir.'/'.$pair.'.merged';
    print STDERR "Processing $pair\n";

    # check if the bam file is already present and skip if it is
    if (-r $bamfn.'.bam') {
	print STDERR "$bamfn.bam exists already. Skipping...\n";
	next;
    }

    # Remove the pair info in the case of the paired end reads, this is done
    # because the gem-2-sam does not recognize the p1 and p2 in the ids
    # Also run the gem-2 sam to generate a sam file with no header
    process_files($lane_files{$pair},
		  $pair,
		  $samfn,
		  $tmpdir);
    # This has to be modified, because if not the file is not filtered and it
    # will crash due to the hard clippings from GEM
    print STDERR "$genomefile.fai not present\n";
    generate_sam_header($genomefile,
			$samfn);
    generate_sorted_bam($samfn,
			$bamfn);
    generate_bam_index($bamfn);
    print join("\t",
	       $pair,
	       $pair.'.merged.bam'),"\n";
}


exit;

sub process_files {
    my $set=shift;
    my $pair=shift;
    my $outfn=shift;
    my $tmpdir=shift;

    my @paired=sort keys %{$set};

    if (@paired == 1) {
	print STDERR "Set identified as single reads\n";
	my $infile=$set->{$paired[0]}[0];
	print STDERR $infile,"\n";
	process_single_reads($infile,
			     "$outfn.$$");
    } elsif (@paired == 2) {
	print STDERR "Set identified as paired reads\n";
	my $infile1=$set->{$paired[0]}[0];
	my $infile2=$set->{$paired[1]}[0];
	print STDERR $infile1,"\n";
	print STDERR $infile2,"\n";
	process_paired_reads($infile1,
			     $infile2,
			     "$outfn.$$",
			     $tmpdir);
    } else {
	# This should never happen
	die "Apparently there are more than two files grouped\n";
    }

    return($outfn);
}

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $genomefile=get_merged_file($lane);
	$lane_files{$pair}{$lane}->[0]=$options{'SAMDIR'}.'/'.$genomefile;
	if (-r $lane_files{$pair}{$lane}->[0].'.gz') {
	    $lane_files{$pair}{$lane}->[0].='.gz';
	}
    }

    return(\%lane_files);
}

# Generate the SAM header and filter the file if required
sub generate_sam_header {
    my $genomefile=shift;
    my $outfn=shift;
    my $maxintronlength=shift || 10000;

    my @seqheader;
    my $invalid=0;
    my $incomplete=0;

    print STDERR "Reading $genomefile...";

    my $infh=Bio::SeqIO->new(-format => 'fasta',
			     -file => $genomefile);
    while (my $seq=$infh->next_seq()) {
	my $id=$seq->display_id();
	my $length=$seq->length();
	my $seqline=join("\t",
			 '@SQ',
			 "SN:$id",
			 "LN:$length");
	push @seqheader,$seqline;
    }
    $infh->close();
    print STDERR "done\n";

    my $outfhtmp=get_fh("$outfn.$$");
    my $outfh=get_fh($outfn,1);

    print $outfh join("\n",
		      @seqheader),"\n";

    while (my $line=<$outfhtmp>) {
	chomp($line);
	# skip empty lines
	if ($line=~/^\s*$/o) {
	    next;
	}

	# Filter invalid cigar lines (the indels that samtools does not
	# recognize hard clippings, etc...)
	my @line=split("\t",$line);

	# Filter for negative insert lengths
	if ($line[5]) {
	    if ($line[5]=~/-/o) {
		$invalid++;
		print STDERR join("\t",
				  "Negative insert:",
				  @line),"\n";
		next;
	    } elsif ($line[5]=~/H/o) {
		# Remove hard clippings
		$invalid++;
		print STDERR join("\t",
				  "Hard clipping:",
				  @line),"\n";
		next;
#	    } elsif ($line[5]=~/^\d+M\d+N\d+M$/o) {
#		my @coords=split(/N/,$line[5]);
#		my $gap=$coords[0];
#		$gap=~s/.+[^\d]//o;
#		$gap=~s/N//o;
#		if ($gap && $gap > $maxintronlength) {
#		    $invalid++;
#		    print STDERR join("\t",
#				      "Gap Too large:",
#				      @line),"\n";
#		    next;
#		}
	    }

	    # This probably should be removed as it is a patch to avoid linking
	    # pairs that are too far appart
#	    if ($line[8] > $maxintronlength) {
#		$line[1]=65;
#		$line[6]='*';
#		$line[7]=0;
#		$line[8]=0;
#	    }
	    print $outfh join("\t",
			      @line),"\n";
	} else {
	    next;
	}
	# If the script gets here the lines are valid
	print $outfh join("\t",
			  @line),"\n";
    }
    close($outfh);
    close($outfhtmp);

    if ($invalid) {
	print STDERR $invalid,"\tInvalid cigar lines found\n";
    }
    if ($incomplete) {
	print STDERR $incomplete,"\tIncomplete lines found. Assumed qualities missing and added a *\n";
    }

    my $command="rm $outfn.$$";
    run_system_command($command);
}

# This subroutine should take a gem mapping file and remove from it the pair
# identifiers so both members of the pair have the same id
sub remove_pair_info {
    my $filename=shift;
    my $tmpdir=shift;
    my $order=shift;

    my $tmpfile=$filename.'.tmp';
    $tmpfile=~s/.*\///o;
    $tmpfile=$tmpdir.'/'.$tmpfile;

    my $infh=get_fh($filename);
    my $outfh=get_fh($tmpfile,1);
    my $count=0;
    my $indels=0;

    print STDERR "Checking Pair ID string in $filename...";

    # This assumes all reads in each file are the same and only checks the first
    my $idok=0;
    my $orderok=0;

    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);
	if ($line[0]=~s/\|p?1$/\/1/o) {
	} elsif ($line[0]=~s/\|p?2$/\/2/o) {
#	} elsif ($line[0]=~s/ 1:/_X:/o) {
#	    # Down syndrome reads
#	    $line[0]=~s/$/\/1/o;
#	} elsif ($line[0]=~s/ 2:/_X:/o) {
#	    # Down syndrome reads
#	    $line[0]=~s/$/\/2/o;
	} else {
	    # This should speed things up
	    $idok=1;
	}

	# Check if the lines files are in the correct order
	my $readnumber;
	if ($line[0]=~/1$/o) {
	    $readnumber=1;
	} elsif ($line[0]=~/2$/o) {
	    $readnumber=2;
	} else {
	    die "Unknown read mate (not 1 or 2)??? in line:\n$line\n";
	}

	if (${$order}) { 
	    if (${$order} != $readnumber) {
		die "Incorrect read ordering $line\n";
	    }
	} else {
	    ${$order}=$readnumber;
	}

	# Check if we have enough fields
	unless ($line[4]) {
	    print STDERR $line,"\n";
	    <STDIN>;
	}

	# Count cases with indels
#	if ($line[4]=~/<[+-][0-9]+>/o) {
#	    $indels++;
#	    print STDERR $line,"\n";
#	    next;
#	}

	# Change any chrMT to chrM
#	if ($line[4]=~/chrMT/o) {
#	    $line[4]=~s/chrMT:/chrM:/og;
#	}

	# Keep only the best match for each entry
	my ($best_hit)=split(',',$line[-1]);
	$line[-1]=$best_hit;

	print $outfh join("\t",
			  @line),"\n";
	$count++;
    }

    close($infh);
    close($outfh);
    print STDERR "done\n";
    print STDERR $count,"\tReads processed\n";
    print STDERR $indels,"\tPresented indels\n";

    return($tmpfile);
}

sub process_paired_reads {
    my $infn1=shift;
    my $infn2=shift;
    my $outfn=shift;
    my $tmpdir=shift;
    my ($order1,$order2);

    # First remove the paired information from the reads and check which is the
    # correct order of the files
    my $infn1tmp=remove_pair_info($infn1,
				  $tmpdir,
				  \$order1);
    my $infn2tmp=remove_pair_info($infn2,
				  $tmpdir,
				  \$order2);

    # Run the gem-2-sam command
    my $command;

    $command ='gem-2-sam ';
    if ($order1 &&
	$order2 &&
	($order2 > $order1)) {
	$command.="-i $infn1tmp ";
	$command.="-ii $infn2tmp ";
    } else {
	print STDERR "Read 1 seems to be in $infn1tmp and read 2 in $infn2tmp.\nI'll invert the order of the files before running gem-2-sam\n";
	$command.="-i $infn2tmp ";
	$command.="-ii $infn1tmp ";
    }
    $command.="--max-pairing-distance 10000 ";
    $command.="-o $outfn > $outfn.log";

    run_system_command($command);

    # clean up
    $command="rm $infn1tmp $infn2tmp";
    run_system_command($command);
}

sub process_single_reads {
    my $infn=shift;
    my $outfn=shift;

    my $command;

    # Unzip to a temporary file
    my $tmpfile=$infn.'.tmp';
    $tmpfile=~s/.*\///o;
    $tmpfile=$tmpdir.'/'.$tmpfile;

    $command="gunzip $infn -c > $tmpfile";
    run_system_command($command);

    $command ='gem-2-sam ';
    $command.="-i $tmpfile ";
    $command.="-o $outfn > $outfn.log";
    run_system_command($command);
    
    # clean up
    $command="rm $tmpfile";
    run_system_command($command);
}
