#!/soft/bin/perl
# DGK

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
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list');
use Bio::SeqIO;

my $genomefile;
my $prefix;
my $tmpdir;
my $samdir;
my $file_list;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$tmpdir=$options{'LOCALDIR'};
$genomefile=$options{'GENOMESEQ'};
$samdir=$options{'SAMDIR'};
$file_list=$options{'FILELIST'};

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
    my $samfn=$tmpdir.'/'.$pair.'.merged.sam';
    my $bamfn=$samdir.'/'.$pair.'.merged';
    print STDERR "Processing $pair\n";

    # check if the bam file is already present
    if (-r $bamfn.'.bam') {
	print STDERR "$bamfn.bam exists already. Skipping...\n";
	next;
    }
    process_files($lane_files{$pair},
		  $pair,
		  $samfn,
		  $tmpdir);
    if (-r $genomefile.'.fai') {
	print STDERR "genomefile.fai is present\n";
	generate_merged_sorted_bam($samfn,
				   $bamfn,
				   $genomefile);
    } else {
	print STDERR "genomefile.fai not present\n";
	generate_sam_header($genomefile,
			    $samfn);
	generate_sorted_bam($samfn,
			    $bamfn);
    }
    generate_bam_index($bamfn);
    print join("\t",
	       $pair,
	       $pair.'.merged.bam'),"\n";
}


exit;

sub generate_bam_index {
    my $bamfile=shift;

    my $command='samtools index ';
    $command.="$bamfile.bam";
    print STDERR "Executing: $command\n";
    system($command);
}

sub generate_sorted_bam {
    my $samfile=shift;
    my $bamfile=shift;

    my $tmpbam=$bamfile;
    $tmpbam=~s/.*\///;
    $tmpbam=$$.'.'.$tmpbam;

    print STDERR "Building sorted BAM file from $samfile\n";

    my $command='samtools view ';
    $command.='-b -S ';
    $command.="$samfile |samtools sort - $bamfile";
    print STDERR "Executing: $command\n";
    system($command);
}

sub generate_merged_sorted_bam {
    my $samfile=shift;
    my $bamfile=shift;
    my $indexfile=shift;

    my $tmpbam=$bamfile;
    $tmpbam=~s/.*\///;
    $tmpbam=$$.'.'.$tmpbam;
    my $tmpsam=$$.'.tmpsam';

    print STDERR "Building sorted BAM file from $samfile\n";

    my $command='samtools import ';
    $command.="$indexfile.fai ";
    $command.="$samfile.$$ $tmpsam;samtools sort $tmpsam $bamfile; rm $tmpsam $samfile.$$";
    print STDERR "Executing: $command\n";
    system($command);
}

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
	if ($line=~/^\s*$/) {
	    next;
	}

	# Filter invalid cigar lines
	my @line=split("\t",$line);

	# Make sure quality field is present
#	if ($line[0]!~/^@/ &&
#	    (!$line[10] ||
#	     $line[10]=~/^\s*$/)) {
#	    $line[10] ='*';
#	    $incomplete++;
#	}

	if ($line[5] &&
	    $line[5]=~/-/) {
	    $invalid++;
	    print STDERR join("\t",
			      @line),"\n";
	} else {
	    print $outfh join("\t",
			      @line),"\n";
	}
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
    system($command);
}

# This subroutine should take a gem mapping file and remove from it the pair
# identifiers so both members of the pair have the same id
sub remove_pair_info {
    my $filename=shift;
    my $tmpdir=shift;

    my $tmpfile=$filename.'.tmp';
    $tmpfile=~s/.*\///;
    $tmpfile=$tmpdir.'/'.$tmpfile;

    my $infh=get_fh($filename);
    my $outfh=get_fh($tmpfile,1);
    my $count=0;
    my $indels=0;

    print STDERR "Checking Pair ID string in $filename...";

    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);
#	$line[0]=~s/(\||\/)p?[12]//;
	$line[0]=~s/\|p?1/\/1/;
	$line[0]=~s/\|p?2/\/2/;

	unless ($line[4]) {
	    print STDERR $line,"\n";
	    <STDIN>;
	}

	# Count cases with indels
	if ($line[4]=~/<[+-][0-9]+>/) {
	    $indels++;
#	    next;
	}

	# Change any chrMT to chrM
	if ($line[4]=~/chrMT/) {
	    $line[4]=~s/chrMT:/chrM:/g;
	}

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

    # First remove the paired information from the reads
    my $infn1tmp=remove_pair_info($infn1,
				  $tmpdir);
    my $infn2tmp=remove_pair_info($infn2,
				  $tmpdir);

    # Run the gem-2-sam command
    my $command;

    $command ='gem-2-sam ';
    $command.="-i $infn1tmp ";
    $command.="-ii $infn2tmp ";
    $command.="-o $outfn > $outfn.log";

    print STDERR "Executing: $command\n";
    system($command);

    # clean up
    $command="rm $infn1tmp $infn2tmp";
    print STDERR "Executing: $command\n";
    system($command);
}

sub process_single_reads {
    my $infn=shift;
    my $outfn=shift;

    my $command;

    $command ='gem-2-sam ';
    $command.="-i $infn ";
    $command.="-o $outfn > $outfn.log";

    print STDERR "Executing: $command\n";
    system($command);

    # clean up
    $command="rm $infn";
    print STDERR "Executing: $command\n";
    system($command);
}
