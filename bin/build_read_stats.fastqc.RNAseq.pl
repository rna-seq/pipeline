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
# This script will take as an input a list of fasta, fastq or bam files and it
# will output a summary of the read statistics in the file calculated by fastqc

# Input: Read files
# Output three stats files with information regarding the read qualities

### TO DO
# Add the median to the qualities per position
# Add the avertage to the N's per position
# Add a sub to guess the quality type and rerun the fastq parsing with the
# correct settings if necessary

# Load some modules
use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use Bio::SeqIO;
use Bio::DB::Sam;
use Tools::FastQC qw(run_fastqc);

# Get some options from the configuration file
my $species;
my $project;
my $tmpdir;
my $file_list;
my $prefix;
my $qualities;
my $debug=0;
my $readdir;
my $genomefn;

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECT'};
$tmpdir=$options{'LOCALDIR'};
$file_list=$options{'FILELIST'};
$prefix=$options{'PREFIX'};
$qualities=$options{'QUALITIES'};
$readdir=$options{'READDIR'};
$genomefn=$options{'GENOMESEQ'};

# Get the files we are going to process
my %files=%{read_file_list()};

my %metainfo;
my %stats;
my %quality_stats;
my %ntpos;


# Get a log file
my $log_fh=get_log_fh('build_read_stats.RNAseq.log',
		      $debug);

foreach my $infile (keys %files) {
    print $log_fh "Processing $infile\n";

    my $lanename=$files{$infile}->[1];

    # Identify the file type
    my $filetype;
    if ($infile=~/.fa(sta)?$/) {
	$filetype='fasta';
    } elsif ($infile=~/.fastq$/) {
	$filetype='fastq';
    } elsif ($infile=~/.bam$/) {
	$filetype='bam';
    } else {
	warn "Unknown filetype for $infile.\nFiletype is guessed from the file ending (fa, fasta or fastq), if the file namings are different please fix them.\n";
	next;
    }
    print $log_fh "$infile identified as $filetype\n";

    my ($good,$bad,$total,$length,$unique);

    # Process tha file
    if ($filetype eq 'fasta') {
	($good,$bad,$total,$length,$unique)=process_fasta_file($infile,
							       \%stats,
							       \%quality_stats,
							       $tmpdir,
							       $log_fh,
							       $lanename);
    } elsif ($filetype eq 'fastq') {
	($good,$bad,$total,$length,$unique)=process_fastq_file($infile,
							       \%stats,
							       \%quality_stats,
							       $tmpdir,
							       $log_fh,
							       $lanename,
							       \%ntpos,
							       $qualities);
    } elsif ($filetype eq 'bam') {
	($good,$bad,$total,$length,$unique)=process_bam_file($infile,
							     \%stats,
							     \%quality_stats,
							     $tmpdir,
							     $log_fh,
							     $lanename,
							     \%ntpos,
							     $qualities,
							     $genomefn);
    } else {
	die "How did I get here???\n";
    }

    $stats{$infile}=[$species,$project,
		     $length, $total,
		     $good,$bad,$unique,$lanename];


    my $fastqcdir=run_fastqc($infile,
			     $readdir,
			     $tmpdir);
}

# Print out the read summaries
print $log_fh 'Building the summary file...';
my $summaryfh=get_fh("${prefix}_read_stats.txt",1);
# Sort on the laneid in order to be consistent with the other sortings
foreach my $file (sort {$stats{$a}->[7] cmp $stats{$b}->[7]} keys %stats) {
    my $filename=$file;
    $filename=~s/.*\///;
    print $summaryfh join("\t",
			  $filename,
			  @{$stats{$file}}),"\n";
}
close($summaryfh);
print $log_fh "done\n";

close($log_fh);

exit;

sub process_fasta_file {
    my $infn=shift;
    my $stats=shift;
    my $quals_pos=shift;
    my $tmpdir=shift;
    my $log_fh=shift;
    my $laneid=shift;
    my $ntpos=shift;

    my $unique=0;
    my $read_length=0;
    my $ambiguous_reads=0;
    my $good_reads=0;
    my $total_reads=0;

    # Open the fasta file
    my $infh=Bio::SeqIO->new(-file => $tmpdir.'/'.$infn,
			     -format => 'fasta');

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    # Process the fasta file
    while (my $seq_obj=$infh->next_seq()) {
	$total_reads++;

	my $seq=$seq_obj->seq;
	# Count the number of unique reads
	print $tmpfh $seq,"\n";

	# Get the minimum read length
	my $seq_length=length($seq);
	if (!$read_length ||
	    ($read_length > $seq_length)) {
	    $read_length=$seq_length;
	}

	# Count the number of Ns in the sequence (ambiguous bases)
	if ($seq=~/N/) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}
    }
    close($tmpfh);
    $infh->close();

   # Get the unique reads without going out of the rood using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}

sub process_fastq_file {
    my $infn=shift;
    my $stats=shift;
    my $quals_pos=shift;
    my $tmpdir=shift;
    my $log_fh=shift;
    my $laneid=shift;
    my $ntpos=shift;
    my $qualities=shift;

    my $unique=0;
    my $read_length=0;
    my $ambiguous_reads=0;
    my $good_reads=0;
    my $total_reads=0;
    my %sequence;
    my $line_type=0;

    # determine the quality variant
    my $variant='illumina';
    if ($qualities eq 'phred') {
	$variant='sanger';
    }

    # Open the fastq file. Here we will use our own sub, because some of the
    # files received do not comply whith the quality values they should actually
    # have and this makes them unreadable by the Bio::SeqIO::fastq module which
    # will always crash
    my $infh;
    $infh=get_fh($tmpdir.'/'.$infn);

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    # Process the fastq file
    while (my $line=<$infh>) {
	chomp($line);

	# Decide which line we are in
	if (($line_type == 0) &&
	    ($line=~s/^@//o)) {
	    # Initialize the sequence hash
	    %sequence=('id' => $line);
	    $line_type=1;
	    $total_reads++;
	    $line=<$infh>;
	    chomp($line);
	} elsif (($line_type==1) &&
		 ($line=~/^\+/o)) {
	    $line_type=2;
	    $line=<$infh>;
	    chomp($line);
	}

	if ($line_type == 1) {
	    $sequence{'seq'}.=$line;
	    next;
	} elsif ($line_type == 2) {
	    $sequence{'qual'}.=$line;
	    $line_type=0;
	    if (length($sequence{'qual'}) < length($sequence{'seq'})) {
		warn "Quality string shorter than sequence,skipping???\n";
		next;
	    } elsif (length($sequence{'qual'}) > length($sequence{'seq'})) {
		warn "Quality string longer than sequence,skipping???\n";
		next;
	    }
	} else {
	    warn "Problem with $infn. Shouldn't be here\n";
	}

	my $seq=$sequence{'seq'};
	# Count the number of unique reads
	print $tmpfh $seq,"\n";

	# Get the minimum read length
	my $seq_length=length($seq);
	if (!$read_length ||
	    ($read_length > $seq_length)) {
	    $read_length=$seq_length;
	}

	# Count the number of sequences with N's (ambiguous bases)
	if ($seq=~/N/o) {
	    $ambiguous_reads++;
	} elsif ($seq=~/\./o) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}
    }
    close($tmpfh);
    $infh->close();

    # Get the unique reads without going out of the roof using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}

sub process_bam_file {
    my $infn=shift;
    my $stats=shift;
    my $quals_pos=shift;
    my $tmpdir=shift;
    my $log_fh=shift;
    my $laneid=shift;
    my $ntpos=shift;
    my $qualities=shift;
    my $genomefn=shift;

    my $unique=0;
    my $read_length=0;
    my $ambiguous_reads=0;
    my $good_reads=0;
    my $total_reads=0;
    my %sequence;
    my $line_type=0;

    # The qualities should be phred in the BAM format
    unless ($qualities eq 'phred') {
	warn "Qualities should be phred in BAM file but are specified as $qualities";
    }

    # Open the BAM file.
    my $sam = Bio::DB::Sam->new(-bam  => $tmpdir.'/'.$infn,
				-fasta=> $genomefn,
	);

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    # Process the BAM file
    # get all the alignments
    my $all_alignments=$sam->features(-iterator => 1);

    while (my $a=$all_alignments->next_seq()) {
	$total_reads++;
	$sequence{'seq'}=$a->query->dna;
	$sequence{'qual'}=$a->qscore;
	if (@{$sequence{'qual'}} < length($sequence{'seq'})) {
		warn "Quality string shorter than sequence,skipping???\n";
		next;
	} elsif (@{$sequence{'qual'}} > length($sequence{'seq'})) {
		warn "Quality string longer than sequence,skipping???\n";
		next;
	}

	my $seq=$sequence{'seq'};
	# Count the number of unique reads
	print $tmpfh $seq,"\n";

	# Get the minimum read length
	my $seq_length=length($seq);
	if (!$read_length ||
	    ($read_length > $seq_length)) {
	    $read_length=$seq_length;
	}

	# Count the number of sequences with N's (ambiguous bases)
	if ($seq=~/N/) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}
    }
    close($tmpfh);

    # Get the unique reads without going out of the roof using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}
