#!/soft/bin/perl

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
# This script should take a gtf/gff containing exons with
# gene names and it will build a set of all exon sequences

use RNAseq_pipeline3 qw(get_fh get_log_fh parse_gff_line);
use RNAseq_pipeline_settings3 qw(read_config_file);
use Bio::SeqIO;
use Bio::DB::Fasta;

my $annotation_file;
my $genome_file;
my $debug=0;

my %options=%{read_config_file()};
$annotation_file=$options{'ANNOTATION'};
$genome_file=$options{'GENOMESEQ'};
my $exonfn=$options{'EXONSFASTA'};

# Get a log file
my $log_fh=get_log_fh('build_exon_seqs.RNAseq.log',
		      $debug);

print $log_fh "Extracting the exon sequences for $annotation_file\n";

unless ($annotation_file && $genome_file) {
    die "ERROR:An annotation file and a genome file are required\n";
}

# Check if the file exists already
my $present=check_file($exonfn);

if ($present) {
    print $log_fh "$exonfn is present. Skipping\n";
    exit;
}

print $log_fh "Processing $annotation_file\n";

# Extract all the exons coordinates belonging to each of the genes
my %exons=%{get_exons_from_gtf($annotation_file)};

# Get the sequence of the chromosome as well as an accessor to extract the exon
# sequence
#my %chromosomes=%{get_sequences($genome_file)};
*get_exon_seq=get_exon_seq_sub($genome_file);

# Get all the exon sequences
my %exon_seqs=%{build_exon_sequences(\%exons)};

# Print out all exon sequences
print $log_fh "Printing exon sequences in $exonfn\n";
my $exonfh=get_fh($exonfn,1);
foreach my $exon (keys %exon_seqs) {
    my $seq=$exon_seqs{$exon};
    print $exonfh ">$exon\n";
    for (my $pos=0;$pos < length($seq) ; $pos += 60) {
	print $exonfh substr($seq,$pos,60),"\n";
    }
}
close($exonfh);
print $log_fh "done\n";

exit;

sub check_file {
    my $file=shift;

    my $present=0;

    my $tablefh=get_fh($options{'PREFIX'}.'_exon_seqs.txt',1);
    if (-r $file) {
	$present=1;
	# Print the table location for the junctions of this experiment
	print $tablefh join("\t",
			    $file,
			    "Present"),"\n";
    } else {
	# Continue, as the table must be created
	print STDERR $file,"\tIs not present\n";
	print $tablefh join("\t",
			    $file,
			    "Generated"),"\n";
    }
    close ($tablefh);

    return($present);
}

# This should accelerate the process of obtaining the exon sequences at the
# expense of RAM
sub get_exon_seq_sub {
    my $genomefile=shift;

    # Generate the chromosomes
    print STDERR "Generating $genomefile index...";
    my $db=Bio::DB::Fasta->new($genomefile);
    print STDERR "done\n";

    my %exons;

    my %revcom=('a' => 't',
		'c' => 'g',
		't' => 'a',
		'g' => 'c',
		'A' => 'T',
		'C' => 'G',
		'T' => 'A',
		'G' => 'C');

    # Build the subroutine for the search
    my $get_seq_from_chrom= sub {
	my $exon_id=shift;
	
	unless (exists $exons{$exon_id}) {
	    my $seq;
	    my @line=split('_',$exon_id);
	    my ($chr,$start,$end,$strand);
	    if (@line == 4) {
		($chr,$start,$end,$strand)=@line;
	    } else {
		my @chr=splice(@line,0,@line - 3,());
		$chr=join('_',@chr);
		($start,$end,$strand)=@line;
	    }

	    if ($strand == 1) {
		$seq=$db->subseq($chr,$start,$end);
	    } elsif ($strand == -1) {
		$seq=$db->subseq($chr,$end,$start)
	    } else {
		warn "unknown strand $strand\n";
		return();
	    }

	    # Attemt to retrieve the sequence removing the initial chr
	    # If the initial multifasta for the genome did not have these
	    # tags the script will not find them if not.
	    my $chr2=$chr;
	    unless ($seq) {
		$chr2=~s/^chr//;
		$seq=$db->subseq($chr2,$end,$start);
	    }
	    # Try again with the mitochondrial sequences
	    unless ($seq) {
		$chr2=~s/M$/MT/;
		$seq=$db->subseq($chr2,$end,$start);
	    }

	    unless ($seq) {
		warn "No sequence retrieved for $exon_id... skipping\n";
		return();
	    }

	    # There is a problem when extracting sequences of length 1 as the
	    # subseq method cannot determine if they are plus or minus, and
	    # decides they are plus
	    # To fix this any case where the sequence is on base and strand is
	    # minus has to be reverse complemented

	    if (($start == $end) &&
		($strand == -1)) {
		if ($revcom{$seq}) {
		    $seq=$revcom{$seq};
		} else {
		    warn "Unknown sequence $seq\n";
		}
	    }

	    $exons{$exon_id}=$seq;
	}
	return($exons{$exon_id});
    };
    return($get_seq_from_chrom);
}

# Get the exon sequences from the coordinates
sub build_exon_sequences {
    my $exons=shift;

    my %exon_seqs;
    my $count=0;

    print STDERR "Building all exon sequences\n";
    foreach my $exon (keys %{$exons}) {
	my $exon_seq=get_exon_seq($exon);
	if ($exon_seq) {
	    $exon_seqs{$exon}=$exon_seq;
	    $count++;
	}
	unless ($count % 1000) {
	    print STDERR "$count\tRetrieved\r";
	}
    }
    if ($count) {
	print STDERR $count,"\tAnnotated exons retrieved\n";
    } else {
	die "Sorry, I couldn't extract any exons...\n";
    }

    return(\%exon_seqs);
}

# This has to be combined with the get_annotation_from gtf for portability
sub get_exons_from_gtf {
    my $file=shift;
    my $log_fh=shift;

    my $fh=get_fh($file);
    my %exons;
    my $count=0;

    # If we have no $log_fh redirect to STDERR
    unless ($log_fh) {
	$log_fh=*STDERR;
    }

    print $log_fh "Reading $file\n";
    print $log_fh "WARNING: Skipping entries located in chr random, chr U(nknown), haplotypes and EnsEMBL assembly exceptions\n";

    while (my $line=<$fh>) {
	$count++;
	chomp($line);

	# Skip possible comments
	if ($line=~/^#/) {
	    next;
	}

	# use the parse_gff  subroutine to parse the lines
	my ($chr,$type,$start,$end,$strand,$frame,$info);
	my %line=%{parse_gff_line($line)};
	$chr=$line{'chr'};
	$type=$line{'type'};
	$start=$line{'start'};
	$end=$line{'end'};
	$frame=$line{'frame'};

	# Skip entries we are not interested in
	# Skip non-exon entries
	unless ($type=~/^exon$/) {
	    next;
	}
	# Skip random unknown and haplotype chromosomes
	if ($chr=~/random/i) {
	    next;
	} elsif ($chr=~/hap/) {
	    next;
	} elsif ($chr=~/^chrU/) {
	    next;
	} elsif ($chr=~/^Un\./) {
	    # This is for EnsEMBL cow
	    next;
	} elsif ($chr=~/^(chr)?HSCHR/) {
	    # This is for human
	    next;
	} elsif ($chr=~/^AAFC03011182/) {
	    # EnsEMBL cow
	    next;
	}

	### TO DO Fix some naming issues that may occurr
	# If the chromosomes are not named as chr in the file name them so this
	# may cause some problems if we are looking at contigs etc... but
	# it should only activate if the chromosomens are named as the humans
	# but with no chr
	if ($chr!~/^(chr|contig|scaffold|supercontig)/) {
	    $chr=~s/^/chr/;
	}
	# To prevent problems with the naming of the chromosomes we will change
	# the chrMT to chrM
	$chr=~s/chrMT/chrM/;

	# Check the strand
	if ($line{'strand'} eq '+') {
	    $strand=1;
	} elsif ($line{'strand'} eq '-') {
	    $strand=-1;
	} else {
	    warn "Unknown strand $strand\n";
	}

	unless ($start && $end && $strand) {
	    die "Missing required information\n";
	}

	my $exon_id=join('_',
			 $chr,
			 $start,
			 $end,
			 $strand);
	$exons{$exon_id}='';
    }
    close($fh);

    $count=keys %exons;
    print STDERR $count, "\tExon entries obtained\n";

    return(\%exons);
}
