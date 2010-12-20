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
# This script should take the output of the GEM mapper
# It will take the read.list file and from there decide if the mappings we are
# looking at are paired or not.

# After this for each genome mapped file it will take the corresponding junction
# mapped file and filter out those possible pseudogene hits (those cases where
# a read hits both genome and junctions)

# It will combine the outputs of both in the following way:
# It will read the files corresponding to the mappings of each pair and from
# those it will select the best hit. The way the best hit is selected is the
# following: If the read maps unsplit with no mismathes it is taken, if the
# read maps unsplit with 1 or 2 mismatchesand maps to the junctions in an
# overlapping possition split, but with no mismatches it is probably an over
# hanging read and the junctions map will be selected
# If only one or the other is present that will be the one selected.

# Once each end of the pair is selected, the best compatible read on the other
# end will be taken. We have to discuss the criteria to this, but in the
# meanwhile the files will be processed as single reads and merged at the end
# with the gem2sam tool (keeping only the unique hits

### TO DO
# The script has to be rewritten to take into account the paired end reads.

# Input is a set of mapping files
# output is a single mapping file
use RNAseq_pipeline3 ('get_fh','get_files_from_table_sub','print_gff');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_lane_id','get_mapping_fh','is_paired');
use RNAseq_GEM3 qw(process_gem_hit parse_gem_line line_reader get_closest_hit);

# Declare some variables
my $paired;
my $samdir;

my %config=%{read_config_file()};
$paired=$config{'PAIRED'};
$samdir=$config{'SAMDIR'};

# Get the subroutines
my %mergers=%{get_mergers()};

# Get the groups of files we will try to merge.
my %files=%{read_file_list()};
my %pairids=%{get_lane_id(\%files)};

# process each group of files
foreach my $pair_id (keys %pairids) {
    # Decide if the reads are paired or not, as that will change how we merge
    # them. This should be done from the read.list.files
    my $paired=is_paired(\%pairids,
			 $pair_id);

    print STDERR "Reads are $paired\n";
    $paired='single';
    print STDERR "Treating reads as $paired\n";

    if ($paired eq 'paired') {
	foreach my $single (@{$pairids{$pair_id}}) {
	    my $mergedfn=$samdir.'/'.$pair_id.'.merged.gem.map.gz';
	    if (-r $mergedfn) {
		print STDERR $mergedfn,"\tIs present. Skipping...\n";
		next;
	    }
	    
	    print STDERR "Processing $paired data for $single\n";
	    print STDERR "Outputting to $mergedfn\n";

	    my $outfh=get_fh($mergedfn,1);
	    $mergers{$paired}->($pairids{$single},
				$pair_id,
				$outfh);
	    close($outfh);
	}
    } elsif ($paired eq 'single') {
	my $mergedfn=$samdir.'/'.$pair_id.'.merged.gem.map.gz';
	if (-r $mergedfn) {
	    print STDERR $mergedfn,"\tIs present. Skipping...\n";
	    next;
	}

	print STDERR "Processing $paired data for $pair_id\n";
	print STDERR "Outputting to $mergedfn\n";

	# Open the output file
	my $outfh=get_fh($mergedfn,1);
	$mergers{$paired}->($pairids{$pair_id},
			    $pair_id,
			    $outfh);
	close($outfh);
    } else {
	die "Unknown read type $paired\n";
    }
}

exit;

### TODO
#  write a merger that will take into account the paired end information
sub get_mergers {
    my %mergers;

    $mergers{'paired'}= sub {
	my $pair=shift;
	my $lane=shift;
	# Here we are using the same entry of lane for both split files
	die "This has to be checked\n";
	my $outfh=shift;

	# Check the reads are actually paired
	if (@{$pair} != 2) {
	    die "We have a problem\n";
	}

	# Get the files forming the first half of the pair
	my $genomemapping1fh=get_mapping_fh($pair->[0],
					    'genome');
	my $junctionsmapping1fh=get_mapping_fh($pair->[0],
					       'junctionsgenomic');
	my $splitmapping1fh=get_mapping_fh($lane,
					   'split');
	
	# get the files forming the second half of the pair
	my $genomemapping2fh=get_mapping_fh($pair->[1],
					    'genome');
	my $junctionsmapping2fh=get_mapping_fh($pair->[1],
					       'junctionsgenomic');
	my $splitmapping2fh=get_mapping_fh($lane,
					   'split');

	# Get the best hit in the following order: genome, junctions,split,
	# recursive
	my $lastsplit1;
	my $lastsplit2;
	my $lastrecurs1;
	my $lastrecurs2;
	my ($genome1,$genome2);
	while (defined($genome1=<$genomemapping1fh>) &&
	       defined($genome2=<$genomemapping2fh>)) {
	    my %linegenome1=%{parse_gem_line($genome1)};
	    my %linegenome2=%{parse_gem_line($genome2)};
	    my $junc1=<$junctionsmapping1fh>;
	    my $junc2=<$junctionsmapping2fh>;
	    my %linejunc1=%{parse_gem_line($junc1)};
	    my %linejunc2=%{parse_gem_line($junc2)};
	    
	    my $read_id1=$linegenome1{'id'};
	    print STDERR $read_id1,"\n";
	}
	### TO DO close all files
    };

    $mergers{'single'}= sub {
	my $pair=shift;
	my $lane=shift;
	my $outfh=shift;

	# Check the reads are actually single
	if (@{$pair} != 1) {
	    die "We have a problem\n";
	}

	# Get the files to merge
	my $genomemappingfh=get_mapping_fh($pair->[0],
					    'genome');
	my $junctionsmappingfh=get_mapping_fh($pair->[0],
					      'junctionsgenomic');
	my $splitmappingfh=get_mapping_fh($lane,
					   'split');

	# Get the best hit in the following order: genome, junctions,split,
	# recursive
	my $lastsplit=<$splitmappingfh>;
	my %lastsplit=%{parse_gem_line($lastsplit)};

	my ($genome);
	my $genomelineno=0;
	my $junctionlineno=0;

	while ($genome=<$genomemappingfh>) {
	    my %linegenome=%{parse_gem_line($genome)};
	    my $junc=<$junctionsmappingfh>;
	    my %linejunc=%{parse_gem_line($junc)};
	    
	    $genomelineno++;
	    $junctionlineno++;

	    # get the read ids
	    my $read_id=$linegenome{'id'};
	    my $read_id_split=$lastsplit{'id'};

	    unless ($linegenome{'id'} eq $linejunc{'id'}) {
		print STDERR join("\t",
				  $genomelineno,
				  $read_id,
				  $junctionlineno,
				  $linejunc{'id'}),"\n";
		die "Incorrect read, are you sure files are sorted\n";
	    }

	    # Get the matches
	    my $genomematches=$linegenome{'matches'};
	    my $juncmatches=$linejunc{'matches'};
	    my $splitmatches=$lastsplit{'matches'};

	    # Choose the best hit:
	    # This criteria is a matter of discussion
	    if ($genomematches!~/^0(:0)*:0$/) {
		# If there is a match to the genome
		if ($juncmatches=~/^0(:0)*:0$/) {
		    # If there is a math to the junctions decide wich is the
		    # best
		    print $outfh $genome;
		} else {
		    # we have both junctions and genome, so we have to choose
		    my @genomematches=split(':',$genomematches);
		    my @juncmatches=split(':',$juncmatches);
		    
		    # Set the maximim number of allowed mismatches to the ones
		    # found in the genome
		    for (my $i=0;$i<@genomematches;$i++) {
			my $genmm=$genomematches[$i];
			my $juncmm=$juncmatches[$i];

			if (($genmm == 0) &&
			    ($juncmm == 0)) {
			    next;
			} elsif ($genmm == 0) {
			    # The best junc hit is better than the best genome
			    # hit
			    print $outfh $junc;
			    last;
			} else {
			    # In this case there are genome hits at least as
			    # good as the junction hits
			    print $outfh $genome;
			    last;
			}
		    }
		}
	    } elsif ($juncmatches!~/^0(:0)*:0$/) {
		# If there is a match to the junctions
		print $outfh $junc;
	    } elsif (($splitmatches!~/^0(:0)*:0$/) &&
		     ($read_id eq $read_id_split)) {
		# If there is a match to the split mapping add qualities and
		# print
		my @lastsplit=split("\t",$lastsplit);
		my $length=length($lastsplit[1]);
		my $quals=substr($linegenome{'qual'},0,$length);
		splice(@lastsplit,2,0,$quals);
		$lastsplit=join("\t",@lastsplit);
		print $outfh $lastsplit;
	    } else {
		# There is no match, so we take the longest
		print $outfh $genome;
	    }

	    if ($read_id eq $read_id_split) {
		my $remaining;
		if (defined($remaining=<$splitmappingfh>)) {
		    $lastsplit=$remaining;
		    %lastsplit=%{parse_gem_line($lastsplit)};
		}
	    }
	}
	close($genomemappingfh);
	close($junctionsmappingfh);
	close($splitmappingfh);
    };

    return(\%mergers);
}
