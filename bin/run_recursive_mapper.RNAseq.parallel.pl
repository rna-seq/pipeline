#!/soft/bin/perl

use strict;
use warnings;

BEGIN {
    unshift @INC, '/users/rg/dgonzalez/lib/Perl';
}

# Objective
# This script will take the necessary commands for running the mapper, and it
# will check if the input files have qualities or not as well as deciding the
# type of qualities at some stage.
# After the first round of mapping it will check the resulting file and remap
# those reads that are unmapped.
# The mapper will first map to the transcriptome and the junctions in order to
# select those reads that are paired uniquely and using this information will
# after map the rest of the reads
# To DO
# Fix so that his mapping is not repeated if the files already exist


use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file);
use RNAseq_GEM3 ('check_index','determine_quality_type','get_mapper_routines',
		 'check_split_input','trim_ambiguous','get_unmapped',
		 'split_map','trim_reads','parse_gem_line');

my $infile;
my $index;
my $outdir;

my $mapper='GEM';
my $threads=8;

my $qualities;
my $ignore_quals;

my $tempdir;
my $filetype;
my $mismatches;
my $readlength;
my $maxmismatch;
my $minmapthreshold=10;
my $trimthreshold=10;
my $debug=0;

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'infile|i=s' => \$infile,
	   'outdir|o=s' => \$outdir,
	   'noqual' => \$ignore_quals);

my %options=%{read_config_file()};
$index=$options{'GENOMEINDEX'};
$mapper=$options{'MAPPER'};
$threads=$options{'THREADS'} unless ($threads);
$tempdir=$options{'LOCALDIR'};
$mismatches=$options{'MISMATCHES'};
$readlength=$options{'READLENGTH'};
$maxmismatch=int($readlength / 25);

# Get a log_fh
my $log_fh=get_log_fh('run_recursive_mapper.log',
		      $debug);

# Check if we have the indices
if ($index && check_index($index)) {
    print $log_fh "Index supplied and valid\n";
} else {
    die "I'm missing a valid index\n";
}

# Check the input file
my $input_ok=check_split_input($infile,
			       \$filetype);
unless ($input_ok) {
    die $infile,"\tIs not readable\n";
}

# check the outdir
unless ($outdir) {
    die "No output directory supplied\n";
}

# Make the outfile
my $outfile=$infile;
$outfile=~s/.*\///;
$outfile=$tempdir.'/'.$outfile;
$outfile=~s/(.fa(stq)*)$/.gem/;

if (-r "$outfile.map" ||
    -r "$outfile.map.gz") {
    print STDERR "$outfile Exists.. Skipping\n";
    exit;
}

if (($filetype eq 'fastq') ||
    $qualities) {
    # The readfile has qualities
    unless ($qualities) {
	print $log_fh 'Guessing quality format...';
	$qualities=determine_quality_type($infile);
	print $log_fh "qualities set to $qualities\n";
    }
}

if ($ignore_quals) {
    print STDERR "Ignoring qualities\n";
    $qualities='ignore';
}

# get the mapper subroutines
my %mapper=%{get_mapper_routines()};

# Run the mapper
$mapper=uc($mapper);
unless (exists $mapper{$mapper}) {
    die "I am not familiar with $mapper. The mapper is not implemented\n";
}

# Before mapping trim the ambiguous nucleotides from the reads, this should
# save some mapping time. All these reads are mapped without quality as they
# are done after the split mapping
my $basename=$infile;
$basename=~s/.*\///;
my $trimmed=$tempdir.'/'.$basename.".noambiguous.$$";
my $left=trim_ambiguous($infile,
			$trimmed,
			$log_fh);

unless($left) {
    die "No unmapped reads found\n";
}

# This should sidestep the multithread problem that may make gem get stuck
# when the number of reads in the mapping file is small
if ($threads > int(($left / 10000) + 0.5)) {
    $threads=int(($left / 10000) + 0.5);
    print $log_fh "Setting threads to $threads\n";
}

my @mapping_files;

# Map with increased mismatches
print $log_fh "Mapping with increased number of mismatches\n";
my $unmapped=$trimmed;

while ($mismatches < $maxmismatch) {    
    $mismatches++;
    print $log_fh "Mapping for $mismatches mismatches...\n";
    my $command;
    my $mapped=$tempdir.'/'.$basename.".mapped.$mismatches.$$";
    print $log_fh "Mapping $unmapped with $mismatches mismatches\n";
    # map the unmapped reads
    $mapper{$mapper}->($index,
		       $unmapped,
		       $mapped,
		       $threads,
		       $qualities,
		       $tempdir,
		       $mismatches);

    $command="rm $unmapped";
    print $log_fh "Executing:\t$command\n";
    system($command);

    # collect the mapped file
    push @mapping_files, $mapped.".map";

    # get the remaining unmapped reads
    $unmapped=$tempdir.'/'.$basename.".unmapped.$mismatches.$$";
    my $left1=get_unmapped($mapped.".map",
			   $unmapped,
			   $log_fh);

    $mapped=$tempdir.'/'.$basename.".split.$mismatches.$$";
    # split map the still unmapped reads
    split_map($index,
	      $unmapped,
	      $mapped,
	      $threads,
	      $qualities,
	      $tempdir,
	      $mismatches);

    $command="rm $unmapped";
    print $log_fh "Executing:\t$command\n";
    system($command);

    # collect the mapped file
    push @mapping_files, $mapped.".split-map";

    # get the remaining unmapped reads
    $unmapped=$tempdir.'/'.$basename.".unmapped.split.$mismatches.$$";
    my $left2=get_unmapped($mapped.".split-map",
			   $unmapped,
			   $log_fh);

    # This should sidestep the multithread problem that may make gem get stuck
    # when the number of reads in the mapping file is small
    if ($threads > int(($left2 / 10000) + 0.5)) {
	$threads=int(($left2 / 10000) + 0.5);
	print $log_fh "Setting threads to $threads\n";
    }
}

# Map with the initial number of mismatches but trimming the reads
$mismatches=$options{'MISMATCHES'};
print $log_fh "Resetting mismatches to $mismatches and performing trimmed mapping\n";
my $round=0;
while (1) {
    my $command;
    print $log_fh "Triming round $round\n";
    my $trimmed=$tempdir.'/'.$basename.".trimmed.$mismatches.$round.$$";
    my $initial=trim_reads($unmapped,
			   $trimmed,
			   $trimthreshold);

    $command="rm $unmapped";
    print $log_fh "Executing:\t$command\n";
    system($command);

    # map the unmapped reads
    my $mapped=$tempdir.'/'.$basename.".mapped.$mismatches.$round.$$";
    $mapper{$mapper}->($index,
		       $trimmed,
		       $mapped,
		       $threads,
		       $qualities,
		       $tempdir,
		       $mismatches);
    $command="rm $trimmed";
    print $log_fh "Executing:\t$command\n";
    system($command);

    # collect the mapped file
    push @mapping_files, $mapped.".map";


    $unmapped=$tempdir.'/'.$basename.".unmapped.$mismatches.$round.$$";
    my $left1=get_unmapped($mapped.".map",
			  $unmapped,
			  $log_fh);

    if ($left1 == 0) {
	last;
    }

    # Split the unmapped reads
    my $split=$tempdir.'/'.$basename.".split.$mismatches.$round.$$";
    split_map($index,
	      $unmapped,
	      $split,
	      $threads,
	      $qualities,
	      $tempdir,
	      $mismatches);

    $command="rm $unmapped";
    print $log_fh "Executing:\t$command\n";
    system($command);

    # collect the mapped file
    push @mapping_files, $split.".split-map";
    
    # get the remaining unmapped reads
    $unmapped=$tempdir.'/'.$basename.".unmapped.split.$mismatches.$round.$$";
    my $left2=get_unmapped($split.".split-map",
			   $unmapped,
			   $log_fh);

    if ($initial == $left2) {
	print join("\t",
		   $initial,
		   $left2),"\n";
	last;
    } elsif ($left2 == 0) {
	last;
    }
    $round++;

    # This should sidestep the multithread problem that may make gem get stuck
    # when the number of reads in the mapping file is small
    if ($threads > int(($left2 / 10000) + 0.5)) {
	$threads=int(($left2 / 10000) + 0.5);
	print $log_fh "Setting threads to $threads\n";
    }
}

# Remove the final file
my $command="rm $unmapped";
print $log_fh "Executing:\t$command\n";
system($command);

my $final_mapping=$outdir.'/'.$basename.".recursive.map";
combine_mapping_files(\@mapping_files,
		      $final_mapping,
		      $tempdir,
		      $log_fh);
close($log_fh);

exit;

sub combine_mapping_files{
    my $files=shift;
    my $outfile=shift;
    my $tmpdir=shift;
    my $log_fh=shift;

    my $final_file="$$.recursive.map";

    print $log_fh "combining the mapped files...\n";
    $ENC{'LC_ALL'}='C';
    my $command='cat ';
    $command.=join(' ',@{$files});
    $command.=" |sort -T $tmpdir| uniq";
    $command.=" > $final_file";
    
    print $log_fh "Executing:\t$command\n";
    system($command);

    # Go through the mapped file and select for each of the reads the longest
    # hit. There should only be one hit in any case, so if there are more we
    # will complain

    my $mappedfh=get_fh($final_file);
    my $outfh=get_fh($outfile,1);

    my %oldread;
    my $oldline='';
    while (my $line=<$mappedfh>) {
	chomp($line);
	my %line=%{parse_gem_line($line)};
	my $id=$line{'id'};
	if (%oldread &&
	    ($oldread{'id'} eq $id)) {
	    # Check if we have matches in the old one
	    if ($oldread{'matches'}=~/^0(:0)*$/) {
		if ($line{'matches'}!~/^0(:0)*$/) {
		    $oldline=$line;
		    %oldread=%line;
		} elsif ($oldread{'length'} < $line{'length'}) {
		    $oldline=$line;
		    %oldread=%line;
		} else {
		    next;
		}
	    } elsif ($line{'matches'}!~/^0(:0)*$/) {
		warn "ERROR multimaps for:\n";
		warn $line,"\n";
		warn $oldline,"\n";
	    }

	} elsif ($oldline) {
	    print $outfh $oldline,"\n";
	    $oldline=$line;
	    %oldread=%line;
	    next;
	} else {
	    $oldline=$line;
	    %oldread=%line;
	    next;
	}
    }
    # Print the last read
    print $outfh $oldline,"\n";

    close($mappedfh);
    close($outfh);

    # clean up
    $command ='rm ';
    $command.=join(' ',@{$files});

    print $log_fh "Executing:\t$command\n";
    system($command);

    $command = "rm $final_file";
    print $log_fh "Executing:\t$command\n";
    system($command);

}
