#!/soft/bin/perl

use strict;
use warnings;

BEGIN {
    unshift @INC, '/users/rg/dgonzalez/lib/Perl';
}

# Objective
# This script will take an input file that may be in one of the following
# formats:
# Solexa _seq.txt
# Solexa _qseq.txt
# fasta
# fastq
# Also it will admit paired ends for each one of these cases
# The file will be processed in the case of the first two to generate a fasta
# or fastq file respectively
# It will not filter ambiguous cases, by default as these will be taken care
# of by mapping. However, if a number of mismatches is given it will fiter
# the reads
# It will print out a summary with the name of the file the number of reads,
# the length of the reads and the number that have no ambiguouos nucleotides
# This will only work correctly if the naming of the file is according to what
# the script expects,in order to guess the format.
# In the case of paired ends it expects both reads to be fused, so if they are
# in separate files they must be trated ad single

### CHANGES NEEDED: Allow to select what directory to use as temporary
###                 Modify the parsers so that all can use tmp (currently not
###                 all can

# Load some modules
use RNAseq_pipeline3 qw(get_fh);
use Getopt::Long;
use Bio::SeqIO;

# Get some command line options and set their defaults
my $species='Unknown species';
my $project='';
my $sizethreshold=2; # This is the maximum allowe filesize if it is higher
                     # the script will use tmp space on the disk
my $mismatches;
my $paired;
my $localdir='';
GetOptions('project|p=s' => \$project,
	   'species|s=s' => \$species,
	   'paired=i' => \$paired,
	   'localdir=s' => \$localdir,
	   'mismatches|m=i' => \$mismatches);

# If mismatches is set, the reads will be filtered

# Paired will indicate paired end reads. This is only appliable to the initial
# case where the file is a txt file, as both pairs will appear as the first and
# second halves of the read.

unless ($species) {
    $species='unknown';
    print STDERR "NO species name provided\n";
}

# Check the localdir
if ($localdir) {
    $localdir=~s/\/$//;
    print STDERR "Using $localdir for tmp files\n";
}

# Process the input
unless (@ARGV) {
    die "Missing input file\n";
}
my %stats;

my %parsing_subs=%{get_parsing_subs(\%stats,
				    $species,
				    $project,
				    $mismatches,
				    $paired)};

while (my $infile=shift) {
    print STDERR "Processing $infile\n";

    my $filetype=guess_file_type($infile);

    # If we cannot identify the file type we will not continue
    unless ($filetype) {
	print STDERR "Please make sure file endings are recognizable, I can only handle _qseq.txt or _seq.txt (and the respective paired end versions)\n";
    }

    my $outfiles=build_outfile_names($infile,
				     $filetype,
				     $paired);

    # Determine if it has to be mapped as paired or not
    if ($paired) {
	$parsing_subs{'paired'}{$filetype}->($infile,
					     @{$outfiles},
					     $paired);
    } else {
	$parsing_subs{'single'}{$filetype}->($infile,
					     $outfiles->[0]);
    }
}

foreach my $file (keys %stats) {
    my $filename=$file;
    $filename=~s/.*\///;
    print join("\t",
	       $filename,
	       @{$stats{$file}}),"\n";
}

exit;

sub build_outfile_names {
    my $infile=shift;
    my $type=shift;
    my $paired=shift;

    my $outfilename=$infile;
    $outfilename=~s/.+\///;
    # This si repeated because it can appear as fastq.txt or similar
    $outfilename=~s/.(txt|fa(stq)*)//;
    $outfilename=~s/.(txt|fa(stq)*)//;
    $outfilename=~s/^((lane(\d)*)(\.|_))*//;
    my $lane=$2;
    if ($lane) {
	$outfilename=$lane.'.'.$outfilename;
    } else {
	$outfilename='lane.'.$outfilename;
    }

    my @outfiles;

    # Build file ending:
    my $ending='.filtered';
    if (($type eq 'solexaSeq') ||
	($type eq 'fasta')){
	$ending.='.fa';
    } elsif (($type eq 'solexaQseq') ||
	     ($type eq 'fastq')) {
	$ending.='.fastq';
    }

    if ($paired) {
	$outfiles[0]=$outfilename.'.1'.$ending;
	$outfiles[1]=$outfilename.'.2'.$ending;
    } else {
	$outfiles[0]=$outfilename.$ending;
    }

    return(\@outfiles);
}

# This subroutine will generate a hash with the subroutines neede to parse the
# different types of files
sub get_parsing_subs {
    my $stats=shift;
    my $species=shift;
    my $project=shift;
    my $mismatches=shift;
    my $paired=shift;

    my %parsing_subs;

    # Build the species name
    my @species=split(/\s+/,$species);
    $species=join('_',@species);

    $parsing_subs{'single'}{'fasta'}=sub {
	my $infn=shift;
	my $outfn=shift;

	my %unique;
	my $unique=0;
	my $read_length=0;
	my $ambiguous_reads=0;
	my $good_reads=0;
	my $total_reads=0;
	
	my $infh=Bio::SeqIO->new(-file => $infn,
				 -format => 'fasta');
	my $outfh=Bio::SeqIO->new(-file => ">$outfn",
				  -format => 'fasta');
	while (my $seq_obj=$infh->next_seq()) {
	    $total_reads++;
	    
	    my $seq=$seq_obj->seq;
	    $unique{$seq}++;
	    unless ($read_length) {
		$read_length=length($seq);
	    }

	    # Count the number of Ns in the sequence (ambiguous bases)
	    my $ambiguous=$seq=~s/N/N/g;
	    if ($ambiguous) {
		$ambiguous_reads++;
		if ($mismatches &&
		    ($ambiguous > $mismatches)) {
		    next;
		}
	    } else {
		$good_reads++;
	    }
	    $outfh->write_seq($seq_obj);
	}
	$infh->close();
	$outfh->close();
	
	$unique=keys %unique;
	$stats->{$infn}=[$species,
			 $project,
			 $read_length,
			 $total_reads,
			 $good_reads,
			 $ambiguous_reads,
			 $unique];
    };
    $parsing_subs{'single'}{'fastq'}=sub {
	my $infn=shift;
	my $outfn=shift;

	my %unique;
	my $unique=0;
	my $read_length=0;
	my $ambiguous_reads=0;
	my $good_reads=0;
	my $total_reads=0;

	my $infh=get_fh($infn);
	my $outfh=get_fh("$outfn.gz",1);

	my $line_type=0;
	my %sequence;
	while (my $line=<$infh>) {
	    chomp($line);

	    # Decide which line we are in
	    if ($line=~s/^@//) {
		%sequence=('id' => $line);
		$line_type=1;
		$total_reads++;
		next;
	    } elsif ($line=~/^\+/) {
		$line_type=2;
		next;
	    }
	    
	    if ($line_type == 1) {
		$sequence{'seq'}.=$line;
		next;
	    } elsif ($line_type == 2) {
		$sequence{'qual'}.=$line;
	    } else {
		warn "Shouldn't be here\n";
	    }
	    
	    my $seq=$sequence{'seq'};
	    $unique{$seq}++;
	    unless ($read_length) {
		$read_length=length($seq);
	    }

	    # Count the number of Ns in the sequence (ambiguous bases)
	    my $ambiguous=$seq=~s/N/N/g;
	    if ($ambiguous) {
		$ambiguous_reads++;
		if ($mismatches &&
		    ($ambiguous > 2 * $mismatches)) {
		    next;
		}
	    } else {
		$good_reads++;
	    }
	    print $outfh '@',$sequence{'id'},"\n";
	    print $outfh "$seq\n";
	    print $outfh "+\n";
	    print $outfh "$sequence{'qual'}\n";
	}
	close($infh);
	close($outfh);

	$unique=keys %unique;
	$stats->{$infn}=[$species,
			 $project,
			 $read_length,
			 $total_reads,
			 $good_reads,
			 $ambiguous_reads,
			 $unique];
	return();
    };
    $parsing_subs{'single'}{'solexaSeq'}=sub {
	my $infn=shift;
	my $outfn=shift;

	my %unique;
	my $unique=0;
	my $read_length=0;
	my $ambiguous_reads=0;
	my $good_reads=0;
	my $total_reads=0;

	my $infh=get_fh($infn);
	my $outfh=get_fh("$outfn.gz",1);
	my $tmpsortfile=$localdir.'/'.$$.'_preproces.RNAseq.txt';
	my $tmpfh;
	if ($localdir) {
	    $tmpfh=get_fh($tmpsortfile,1);
	}
	while (my $line=<$infh>) {
	    $total_reads++;
	    chomp($line);
	    
	    my @line=split(/\t/,$line);
	    
	    my $id=join('_',
			$species,
			@line[0..3]);
	    my $seq=$line[4];

	    if ($localdir) {
		print $tmpfh $seq,"\n";
	    } else {
#		$unique{$seq}++;
	    }

	    unless ($read_length) {
		$read_length=length($seq);
	    }
	    
	    my $ambiguous=$seq=~s/\./N/g;
	    if ($ambiguous) {
		$ambiguous_reads++;
		if ($mismatches &&
		    ($ambiguous > $mismatches)) {
		    next;
		}
	    } else {
		$good_reads++;
	    }
	    print $outfh ">$id\n";
	    print $outfh "$seq\n";
	}
	close($infh);
	close($outfh);

	if ($localdir) {
	    close($tmpfh);
	    $unique=`sort $tmpsortfile | uniq | wc -l`;
	    my $command="rm $tmpsortfile";
	    system($command);
	} else {
	    $unique=keys %unique;
	    chomp($unique);
	}

	$stats->{$infn}=[$species,
			 $project,
			 $read_length,
			 $total_reads,
			 $good_reads,
			 $ambiguous_reads,
			 $unique];

	return();
    };
    $parsing_subs{'single'}{'solexaQseq'}=sub {
	my $infn=shift;
	my $outfn=shift;

	my %unique;
	my $unique=0;
	my $read_length=0;
	my $ambiguous_reads=0;
	my $good_reads=0;
	my $total_reads=0;
	
	my $infh=get_fh($infn);
	my $outfh=get_fh("$outfn.gz",1);
	my $tmpsortfile=$localdir.'/'.$$.'_preproces.RNAseq.txt';
	my $tmpfh;
	if ($localdir) {
	    $tmpfh=get_fh($tmpsortfile,1);
	}
	while (my $line=<$infh>) {
	    $total_reads++;
	    chomp($line);
	    
	    my @line=split(/\t/,$line);
	    
	    my $id=join('_',
			$species,
			@line[2..5]);
	    $id.='|p'.$line[7];
	    my $seq=$line[8];
	    my $qual=$line[9];

	    if ($localdir) {
		print $tmpfh $seq,"\n";
	    } else {
		$unique{$seq}++;
	    }

	    unless ($read_length) {
		$read_length=length($seq);
	    }
	    
	    my $ambiguous=$seq=~s/\./N/g;
	    if ($ambiguous) {
		$ambiguous_reads++;
		if ($mismatches &&
		    ($ambiguous > 2 * $mismatches)) {
		    next;
		}
	    } else {
		$good_reads++;
	    }
	    print $outfh '@',"$id\n";
	    print $outfh "$seq\n";
	    print $outfh '+',"$id\n";
	    print $outfh "$qual\n";
	}
	close($infh);
	close($outfh);

	if ($localdir) {
	    close($tmpfh);
	    $unique=`sort $tmpsortfile | uniq | wc -l`;
	    my $command="rm $tmpsortfile";
	    system($command);
	} else {
	    $unique=keys %unique;
	    chomp($unique);
	}

	$stats->{$infn}=[$species,
			 $project,
			 $read_length,
			 $total_reads,
			 $good_reads,
			 $ambiguous_reads,
			 $unique];

	return();
    };
    $parsing_subs{'paired'}{'fasta'}=sub {die "No sub\n";};
    $parsing_subs{'paired'}{'fastq'}=sub {
	my $infn=shift;
	my $outfn1=shift;
	my $outfn2=shift;
	my $read_length=shift;

	my %unique1;
	my %unique2;
	my $ambiguous_reads1=0;
	my $ambiguous_reads2=0;
	my $good_reads1=0;
	my $good_reads2=0;
	my $total_reads=0;
	
	my $infh=get_fh($infn);
	my $outfh1=get_fh("$outfn1.gz",1);
	my $outfh2=get_fh("$outfn2.gz",1);

	my $line_type=0;
	my %sequence;
	while (my $line=<$infh>) {
	    chomp($line);

	    # Decide which line we are in
	    if ($line=~s/^@//) {
		%sequence=('id' => $line);
		$line_type=1;
		$total_reads++;
		next;
	    } elsif ($line=~/^\+/) {
		$line_type=2;
		next;
	    }
	    
	    if ($line_type == 1) {
		$sequence{'seq'}.=$line;
		next;
	    } elsif ($line_type == 2) {
		$sequence{'qual'}.=$line;
	    } else {
		warn "Shouldn't be here\n";
	    }

	    my $id=$sequence{'id'};
	    
	    my $seq2=$sequence{'seq'};
	    my $qual2=$sequence{'qual'};

	    # Now parse independently the two halves
	    my $seq1=substr($seq2,0,$read_length,'');
	    my $qual1=substr($qual2,0,$read_length,'');
	    $unique1{$seq1}++;
	    $unique2{$seq2}++;

	    my $ambiguous1=$seq1=~s/\./N/g;
	    my $ambiguous2=$seq2=~s/\./N/g;
	    my $keep1=1;
	    my $keep2=1;
	    if ($ambiguous1) {
		$ambiguous_reads1++;
		if ($mismatches &&
		    ($ambiguous1 > 2 * $mismatches)) {
		    $keep1=0;
		}
	    } else {
		$good_reads1++;
	    }
	    
	    if ($ambiguous2) {
		$ambiguous_reads2++;
		if ($mismatches &&
		    ($ambiguous2 > $mismatches)) {
		    $keep2=0;
		}
	    } else {
		$good_reads2++;
	    }

	    if ($keep1) {
		print $outfh1 '@'."$id|1\n";
		print $outfh1 "$seq1\n";
		print $outfh1 '+'."$id|1\n";
		print $outfh1 "$qual1\n";
	    }
	    if ($keep2) {
		print $outfh2 '@'."$id|2\n";
		print $outfh2 "$seq2\n";
		print $outfh2 '+'."$id|2\n";
		print $outfh2 "$qual2\n";
	    }
	}
	close($infh);
	close($outfh1);
	close($outfh2);
    
	my $unique1=keys %unique1;
	my $unique2=keys %unique2;
	$stats->{$infn.'.1'}=[$species,
			      $project,
			      $read_length,
			      $total_reads,
			      $good_reads1,
			      $ambiguous_reads1,
			      $unique1];
	$stats->{$infn.'.2'}=[$species,
			      $project,
			      $read_length,
			      $total_reads,
			      $good_reads2,
			      $ambiguous_reads2,
			      $unique2];
	return();
    };
    $parsing_subs{'paired'}{'solexaSeq'}=sub {
	my $infn=shift;
	my $outfn1=shift;
	my $outfn2=shift;
	my $read_length=shift;

	my %unique1;
	my %unique2;
	my $ambiguous_reads1=0;
	my $ambiguous_reads2=0;
	my $good_reads1=0;
	my $good_reads2=0;
	my $total_reads=0;

	my $infh=get_fh($infn);
	my $outfh1=get_fh("$outfn1.gz",1);
	my $outfh2=get_fh("$outfn2.gz",1);
	while (my $line=<$infh>) {
	    $total_reads++;
	    chomp($line);
	    
	    my @line=split(/\t/,$line);

	    # Prevent downstream parsing problems due to whitespace in the
	    # species name
	    my $id=join('_',
			$species,
			@line[0..3]);
	    my $seq2=$line[4];
	
	    # Now parse independently the two halves
	    my $seq1=substr($seq2,0,$read_length,'');
	    $unique1{$seq1}++;
	    $unique2{$seq2}++;

	    my $ambiguous1=$seq1=~s/\./N/g;
	    my $ambiguous2=$seq2=~s/\./N/g;
	    my $keep1=1;
	    my $keep2=1;
	    if ($ambiguous1) {
		$ambiguous_reads1++;
		if ($mismatches &&
		    ($ambiguous1 > $mismatches)) {
		    $keep1=0;
		}
	    } else {
		$good_reads1++;
	    }
	    
	    if ($ambiguous2) {
		$ambiguous_reads2++;
		if ($mismatches &&
		    ($ambiguous2 > $mismatches)) {
		    $keep2=0;
		}
	    } else {
		$good_reads2++;
	    }

	    if ($keep1) {
		print $outfh1 ">$id/1\n";
		print $outfh1 "$seq1\n";
	    }
	    if ($keep2) {
		print $outfh2 ">$id/2\n";
		print $outfh2 "$seq2\n";
	    }
	}
	close($infh);
	close($outfh1);
	close($outfh2);
    
	my $unique1=keys %unique1;
	my $unique2=keys %unique2;
	$stats->{$infn.'.1'}=[$species,
			      $project,
			      $read_length,
			      $total_reads,
			      $good_reads1,
			      $ambiguous_reads1,
			      $unique1];
	$stats->{$infn.'.2'}=[$species,
			      $project,
			      $read_length,
			      $total_reads,
			      $good_reads2,
			      $ambiguous_reads2,
			      $unique2];
	return();
    };
    $parsing_subs{'paired'}{'solexaQseq'}=sub {die "No sub\n";};

    return(\%parsing_subs);
}

sub guess_file_type {
    my $infile=shift;
    my $filetype;

    print STDERR 'Guessing file type from name...';
    # Identify the file type
    if ($infile=~/_seq.txt(.gz)?$/) {
	$filetype='solexaSeq';
    } elsif ($infile=~/_qseq.txt(.gz)?$/) {
	$filetype='solexaQseq';
    } elsif ($infile=~/.fa(.gz)$/) {
	$filetype='fasta';
    } elsif ($infile=~/.fastq(.gz)$/) {
	$filetype='fastq';
    } else {
	warn "Unknown filetype for $infile\n";
	next;
    }

    print STDERR "done\n";
    print STDERR "Type determined as $filetype\n";

    return($filetype);
}
