package Tools::Bam;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('process_bam_file','bam2sequence','bam2coords','check_mate_pair',
	    'generate_bam_index','generate_sorted_bam','add_tag');

use strict;
use warnings;

# This package uses the Bioper interface to the samtools in order to query Bam
# and SAM files

use RNAseq_pipeline3 qw(get_fh run_system_command);
use Bio::DB::Sam;
use Bio::SeqIO;
use Bio::Seq::Quality;

sub generate_sorted_bam {
    my $samfile=shift;
    my $bamfile=shift;

    my $tmpbam=$bamfile;
    $tmpbam=~s/.*\///o;
    $tmpbam=$$.'.'.$tmpbam;

    print STDERR "Building sorted BAM file from $samfile\n";

    # Sam tools outputs many times some of the reads and we do not want that
    my $command="uniq $samfile | samtools view ";
    $command.='-b -S ';
    $command.="- |samtools sort - $bamfile";
    run_system_command($command);
}

sub add_tag {
    my $line=shift;
    my $tag=shift;
    my @tag=split(':',$tag);
    my $tag_prefix=join(':',@tag[0,1]);

    chomp($line);
    my @line=split("\t",$line);
    my $flags=$line[11];
    if ($flags) {
	unless ($flags=~/$tag_prefix/) {
	    $line[11]=join(' ',$flags,$tag);
	}
    } else {
	$line[11]=$tag;
    }
    $line=join("\t",@line);

    return($line);
}

sub generate_bam_index {
    my $bamfile=shift;

    my $command='samtools index ';
    $command.="$bamfile.bam";
    run_system_command($command);
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

    my %ntindex=('A' => 1,
		 'C' => 2,
		 'G' => 3,
		 'T' => 4);

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

	# Build the distribution of the N's and qualities per position
	my @nucleotides=split('',$seq);
	my @qualities=@{$sequence{'qual'}};
	unless (scalar(@nucleotides) == scalar(@qualities)) {
	    die "Quality and nucleotide lengths differ\n";
	}

	# Using a single loop to go through both nt and qualities should save
	# processing time
	for (my $i=0;$i<@nucleotides;$i++) {
	    my $pos=$i + 1;

	    # Set the qualities
	    my $qual=$qualities[$i];
	    $quals_pos->{$laneid}->{$pos}+=$qual;

	    # Set the nucleotides
	    if ($nucleotides[$i]=~/[^ACTG]/o) {
		$ntpos->{$laneid}->{$pos}->[0]++;
	    } else {
		$ntpos->{$laneid}->{$pos}->[$ntindex{$nucleotides[$i]}]++;
		$ntpos->{$laneid}->{$pos}->[0]+=0;
	    }
	}
    }
    close($tmpfh);

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

sub bam2sequence {
    my $infn=shift;

    my %sequence;

    # Open the BAM file.
    my $sam = Bio::DB::Sam->new(-bam  => $infn);

    my $outfn=$infn;

    # Process the BAM file
    # get all the alignments
    my $all_alignments=$sam->features(-iterator => 1);

    # Check the first line to see if the file has qualities or not.
    my $a1=$all_alignments->next_seq();

    my $seq= Bio::Seq::Quality->new();
	
    $seq->id($a1->name);
    $seq->seq($a1->query->dna);
    $seq->qual([$a1->qscore]);

    # Get the quality string from the values
    my $outfh;
    my $filetype;
    if (@{$seq->qual} != length($seq->seq())) {
	# The file should be extracted as fasta
	print STDERR "Extracting file as fasta, as it is missing qualities\n";
	$outfn=~s/.bam/.fa/;
	$filetype='fasta';
    } else {
	# The file should be extracted as fastq
	print STDERR "Extracting file as fastq\n";
	$outfn=~s/.bam$/.fastq/;
	$filetype='fastq';
    }

    if (-r $outfn) {
	print STDERR $outfn,"\tIs present already...Skipping\n";
    } else {
	$outfh=Bio::SeqIO->new(-format    => $filetype,
			       -file      => ">$outfn");
	# Print the first entry:
	$outfh->write_seq($seq);
	
	while (my $a=$all_alignments->next_seq()) {
	    my $seq= Bio::Seq::Quality->new();
	    
	    $seq->id($a->name());
	    $seq->seq($a->query->dna());
	    $seq->qual([$a->qscore()]);
	    
	    # Print the entry:
	    $outfh->write_seq($seq);
	}
	$outfh->close();
    }
    return($outfn);
}

sub process_aln {
    my $aln=shift;
    my $entry={};

    $entry->{'chr'}=$aln->seq_id();
    $entry->{'start'}=$aln->start();
    $entry->{'end'}=$aln->end();
    $entry->{'strand'}=$aln->strand();
    # Correct strand format
    if ($entry->{'strand'} == 1) {
	$entry->{'strand'}='+';
    } elsif ($entry->{'strand'} == -1) {
	$entry->{'strand'}='-';
    } else {
	$entry->{'strand'}='.';
    }
    return($entry);
}

# Use a different iterator from the main one when using this, as if trying to
# the same one it does not work (also this new one can be called without
# expansion of flags or splices
sub bam2coords {
    my $aln=shift; # This should be a Bio::DB::BAM::Alignment object
    my $global=shift; # This should be a different iterator from the main one
    my $entry={};
    
    my @coords;
    my $paired=$aln->get_tag_values('PAIRED');

    my @subaligns=$aln->get_SeqFeatures();

    # Check if the alignment is spliced
    if (@subaligns) {
	while (my $a=shift(@subaligns)) {
	    my $entry=process_aln($a);
	    $entry->{'id'}=$aln->name();
	    
	    # Add the paired end info to the IDs
	    if ($paired) {
		if ($aln->get_tag_values('SECOND_MATE')) {
		    $entry->{'id'}.='/2';
		} else {
		    $entry->{'id'}.='/1';
		}
	    }
	    # Set the spliced
	    $entry->{'spliced'}=1;
	    
	    # set the unique
	    $entry->{'unique'}=0;
	    $entry->{'matches'}=join(':',
				     $aln->get_tag_values('H0'),
				     $aln->get_tag_values('H1'),
				     $aln->get_tag_values('H2'));
	    # Determine if it is a unique map.
	    if ($entry->{'matches'}=~/^(0:)*1(:.)*/) {
		$entry->{'unique'}=1;
	    } elsif ($entry->{'matches'}) {
		# This is a multimap
		$entry->{'unique'}=0;
	    } elsif ($global) {
		$entry->{'unique'}=check_unique($aln,
						$global);
	    } else {
		# The tags are not defined so we don't know what this is
		my $read_id=$entry->{'id'};
		print STDERR "No H? tag found for $read_id. Setting as unique, but beware\n";
		$entry->{'unique'}=1;
	    }
	    
	    # Add the entry if it has a match
	    $entry->{'cigar'}=$aln->cigar_str();
	    if ($entry->{'cigar'} &&
		$entry->{'cigar'}!~/\*/) {
		push @coords,$entry;
	    }
	}
    } else {
	# The alignment is not spliced
	my $entry=process_aln($aln);
	$entry->{'id'}=$aln->name();
	    
	# Add the paired end info to the IDs
	if ($paired) {
	    if ($aln->get_tag_values('SECOND_MATE')) {
		$entry->{'id'}.='/2';
	    } else {
		$entry->{'id'}.='/1';
	    }
	}
	$entry->{'spliced'}=0;
	    
	# set the unique
	$entry->{'unique'}=0;
	$entry->{'matches'}=join(':',
				 $aln->get_tag_values('H0'),
				 $aln->get_tag_values('H1'),
				 $aln->get_tag_values('H2'));
	# Determine if it is a unique map.
	if ($entry->{'matches'}=~/^(0:)*1(:.)*/) {
	    $entry->{'unique'}=1;
	} elsif ($entry->{'matches'}) {
	    # This is a multimap
	    $entry->{'unique'}=0;
	} elsif ($global) {
	    $entry->{'unique'}=check_unique($aln,
					    $global);
	} else {
	    # The tags are not defined so we don't know what this is
	    my $read_id=$entry->{'id'};
	    print STDERR "No H? tag found for $read_id and I have no reference alignment. Setting as unique, but beware\n";
	    $entry->{'unique'}=1;
	}
	    
	# Determine if it is spliced
	$entry->{'cigar'}=$aln->cigar_str();
	if ($entry->{'cigar'} &&
	    $entry->{'cigar'}!~/\*/) {
	    push @coords,$entry;
	}
    }
    return(\@coords);
}

sub check_unique {
    my $aln=shift;
    my $global=shift;

    my $seq=$aln->query->dna();
    my $seq_id=$aln->name();
    my $iterator=$global->features(-name => $seq_id,
				   -iterator => 1);
    my %hits;
    $hits{$seq}=0;
    my $unique=1;
    while (my $align=$iterator->next_seq()) {
	my $sequence=$align->query->dna();
	$hits{$sequence}++;
	if ($hits{$seq} > 1) {
	    $unique=0;
	    last;
	}
    }
    return($unique)
}

# This should help determine if a read is unique and what we see are the two
# mates with the same ID
sub check_mate_pair {
    my $mate1=shift;
    my $mate2=shift;
    my @mate1=split("\t",$mate1);
    my @mate2=split("\t",$mate2);

    my $pairing_ok=0;

    if (($mate1[3] == $mate2[7]) &&
	($mate1[7] == $mate2[3]) &&
	(abs($mate1[8])== abs($mate2[8]))) {
	$pairing_ok=1;
    }

    return($pairing_ok);
}

1;
