package Tools::Bam;

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('process_bam_file');

use strict;
use warnings;

# This package uses the Bioper interface to the samtools in order to query Bam
# and SAM files

use RNAseq_pipeline3 qw(get_fh run_system_command);

sub process_bam_file {
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

    # Open the bam file. 
    my $infh;
    $infh=get_fh($tmpdir.'/'.$infn);

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    my %ntindex=('A' => 1,
		 'C' => 2,
		 'G' => 3,
		 'T' => 4);

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
	if ($seq=~/N/) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}

	# Build the distribution of the N's and qualities per position
	my @nucleotides=split('',$seq);
	my @qualities=split('',$sequence{'qual'});
	unless (scalar(@nucleotides) == scalar(@qualities)) {
	    die "Quality and nucleotide lengths differe\n";
	}

	# Using a single loop to go through both nt and qualities should save
	# processing time
	for (my $i=0;$i<@nucleotides;$i++) {
	    my $pos=$i + 1;

	    # Set the qualities
	    my $qual=ord($qualities[$i]);
	    if ($variant eq 'sanger') {
		$qual-=33;
	    } else {
		$qual-=64;
	    }
	    $quals_pos->{$laneid}->{$pos}+=$qual;

	    # Set the nucleotides
	    if ($nucleotides[$i]=~/[^ACTG]/o) {
		$ntpos->{$laneid}->{$pos}->[0]++;
	    } else {
		$ntpos->{$laneid}->{$pos}->[$ntindex{$nucleotides[$i]}]++;
	    }
	}
    }
    close($tmpfh);
    $infh->close();

    ### TO DO
    # Calculate the average qualities.
    # This is only useful if we know the type
    # of qualities We should guess the quality type here and after that print
    # the actual meaning of the qualities
    foreach my $pos (keys %{$quals_pos->{$laneid}}) {
	$quals_pos->{$laneid}->{$pos}=sprintf "%.2f",($quals_pos->{$laneid}->{$pos} / $total_reads);
    }

    # Get the unique reads without going out of the roof using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}

1;
