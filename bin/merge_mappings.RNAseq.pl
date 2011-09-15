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

# The script will separate the reads into unique amd multimapping reads

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
my %pairids=%{get_lane_id()};

# process each group of files
foreach my $pair_id (keys %pairids) {
    # Decide if the reads are paired or not, as that will change how we merge
    # them. This should be done from the read.list.files
    my $paired=is_paired($pair_id);

    print STDERR "Reads are $paired\n";
    print STDERR "Treating reads as $paired\n";

    if ($paired eq 'paired') {
	my %outfh1;
	my %outfh2;
	foreach my $single (keys %{$pairids{$pair_id}}) {
	    my $mergedfn1=$samdir.'/'.$single.'.unique.merged.gem.map.gz';
	    my $mergedfn2=$samdir.'/'.$single.'.multi.merged.gem.map.gz';
	    
	    if (-r $mergedfn1) {
		print STDERR $mergedfn1,"\tIs present. Skipping...\n";
		next;
	    }
	    
	    print STDERR "Processing $paired data for $single\n";
	    print STDERR "Outputting to $mergedfn1 and mergedfn2\n";

	    $outfh1{$single}=get_fh($mergedfn1,1);
	    $outfh2{$single}=get_fh($mergedfn2,1);
	}

	$mergers{$paired}->([keys %{$pairids{$pair_id}}],
			    \%outfh1,
			    \%outfh2,
			    $pairids{$pair_id});

	# Close the filehandles
	foreach my $single (keys %outfh1) {
	    close($outfh1{$single});
	    close($outfh2{$single});
	}
    } elsif ($paired eq 'single') {
	my $mergedfn1=$samdir.'/'.$pair_id.'.unique.merged.gem.map.gz';
	my $mergedfn2=$samdir.'/'.$pair_id.'.multi.merged.gem.map.gz';
	if (-r $mergedfn1) {
	    print STDERR $mergedfn1,"\tIs present. Skipping...\n";
#	    next;
	}

	print STDERR "Processing $paired data for $pair_id\n";
	print STDERR "Outputting to $mergedfn1\n";

	# Open the output file
	my $outfh1=get_fh($mergedfn1,1);
	my $outfh2=get_fh($mergedfn2,1);
	$mergers{$paired}->([keys %{$pairids{$pair_id}}],
			    $outfh1,
			    $outfh2,
			    $pairids{$pair_id});
	close($outfh1);
	close($outfh2);
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
	my $uniqfh=shift,
	my $multifh=shift;
	my $pairs=shift;

	# Here we are using the same entry of lane for both split files

	# Check the reads are actually paired
	if (@{$pair} != 2) {
	    die "We have a problem\n";
	}

	# Get the files forming the first half of the pair
	my $genomemapping1fh=get_mapping_fh($pairs->{$pair->[0]},
					    'genome');
	my $junctionsmapping1fh=get_mapping_fh($pairs->{$pair->[0]},
					       'junctionsgenomic');
	my $splitmapping1fh=get_mapping_fh($pair->[0],
					   'split');
	my $recmapping1fh=get_mapping_fh($pair->[0],
					'recursive');
	my $uniqfh1=$uniqfh->{$pair->[0]};
	my $multifh1=$multifh->{$pair->[0]};

	*get_best_map1=get_best_map_sub($genomemapping1fh,
					$junctionsmapping1fh,
					$splitmapping1fh,
					$recmapping1fh);
	
	# get the files forming the second half of the pair
	my $genomemapping2fh=get_mapping_fh($pairs->{$pair->[1]},
					    'genome');
	my $junctionsmapping2fh=get_mapping_fh($pairs->{$pair->[1]},
					       'junctionsgenomic');
	my $splitmapping2fh=get_mapping_fh($pair->[1],
					   'split');
	my $recmapping2fh=get_mapping_fh($pair->[1],
					'recursive');
	my $uniqfh2=$uniqfh->{$pair->[1]};
	my $multifh2=$multifh->{$pair->[1]};

	*get_best_map2=get_best_map_sub($genomemapping2fh,
					$junctionsmapping2fh,
					$splitmapping2fh,
					$recmapping2fh);

	while (my $record1=get_best_map1()) {
	    my $record2=get_best_map2() || die "Missing entry\n";
	    my $type1=classify_pair($record1);
	    my $type2=classify_pair($record2);
	    if (($type1 eq $type2) &&
		($type1 eq 'unique')) {
		print $uniqfh1 $record1;
		print $uniqfh2 $record2;
	    } else {
		print $multifh1 $record1;
		print $multifh2 $record2;
	    }
	}

	close($genomemapping1fh);
	close($junctionsmapping1fh);
	close($splitmapping1fh);
	close($recmapping1fh);
	close($genomemapping2fh);
	close($junctionsmapping2fh);
	close($splitmapping2fh);
	close($recmapping2fh);
    };

    $mergers{'single'}= sub {
	my $pair=shift;
	my $uniquefh=shift;
	my $multifh=shift;
	my $pairs=shift;

	my $lane;
	# Check the reads are actually single
	if (@{$pair} != 1) {
	    die "We have a problem\n";
	} else {
	    ($lane)=@{$pair};
	}
	print STDERR $lane,"\n";

	# Get the files to merge
	my $genomemappingfh=get_mapping_fh($pairs->{$lane},
					    'genome');
	my $junctionsmappingfh=get_mapping_fh($pairs->{$lane},
					      'junctionsgenomic');
	my $splitmappingfh=get_mapping_fh($lane,
					   'split');
	my $recmappingfh=get_mapping_fh($lane,
					'recursive');

	*get_best_map=get_best_map_sub($genomemappingfh,
				       $junctionsmappingfh,
				       $splitmappingfh,
				       $recmappingfh);
	while (my $record=get_best_map()) {
	    my $type=classify_pair($record);
	    if ($type eq 'unique') {
		print $uniquefh $record;
	    } else {
		print $multifh $record;
	    }
	}

	close($genomemappingfh);
	close($junctionsmappingfh);
	close($splitmappingfh);
	close($recmappingfh);
    };

    return(\%mergers);
}

sub get_best_map_sub {
    my $genomefh=shift;
    my $junctionsfh=shift;
    my $splitfh=shift;
    my $recfh=shift;

    my $lastsplit=<$splitfh>  || print STDERR "Uninitialized split\n";
    my %lastsplit=%{parse_gem_line($lastsplit)};

    my $lastrec=<$recfh> || print STDERR "Uninitialized rec\n";
    my %lastrec=%{parse_gem_line($lastsplit)};

    my $genomelineno=0;
    my $junctionlineno=0;

    my $bestmap_sub=sub {
	# Get the best hit in the following order: genome, junctions,split,
	# recursive
	
	my $genome=<$genomefh> || return();
	my %linegenome=%{parse_gem_line($genome)};
	my $junc=<$junctionsfh>;
	my %linejunc=%{parse_gem_line($junc)};
	    
	$genomelineno++;
	$junctionlineno++;
	    
	# get the read ids
	my $read_id=$linegenome{'id'};
	my $read_id_split=$lastsplit{'id'};
	my $read_id_rec=$lastrec{'id'};

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
	my $recmatches=$lastrec{'matches'};
	    
	# Choose the best hit:
	my $record;
	# This criteria is a matter of discussion
	if ($genomematches!~/^0(:0)*:0$/) {
	    # If there is a match to the genome
	    if ($juncmatches=~/^0(:0)*:0$/) {
		# If there is a match to the junctions decide wich is the
		# best
		$record=$genome;
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
			$record=$junc;
			last;
		    } else {
			# In this case there are genome hits at least as
			# good as the junction hits
			$record=$genome;
			last;
		    }
		}
	    }
	} elsif ($juncmatches!~/^0(:0)*:0$/) {
	    # If there is a match to the junctions
	    $record=$junc;
	} elsif (($splitmatches!~/^0(:0)*:0$/) &&
		 ($read_id eq $read_id_split)) {
	    # If there is a match to the split mapping add qualities and
	    # print
	    my @lastsplit=split("\t",$lastsplit);
	    my $length=length($lastsplit[1]);
	    
	    # Include the qualities only if present
	    if ($linegenome{'qual'}) {
		my $quals=substr($linegenome{'qual'},0,$length);
		splice(@lastsplit,2,0,$quals);
	    }
	    
	    $lastsplit=join("\t",@lastsplit);
	    $record=$lastsplit;
	} elsif (($recmatches!~/^0(:0)*:0$/) &&
		 ($read_id eq $read_id_rec)) {
	    # If there is a match to the recursive mapping add qualities and
	    # print
	    my @lastrec=split("\t",$lastrec);
	    my $length=length($lastrec[1]);
	    
	    # Include the qualities only if present
	    if ($linegenome{'qual'}) {
		my $quals=substr($linegenome{'qual'},0,$length);
		splice(@lastrec,2,0,$quals);
	    }
	    
	    $lastrec=join("\t",@lastrec);
	    $record=$lastrec;
	} elsif ($genome) {
	    # There is no match, so we take the longest
	    $record=$genome;
	}

	# Update splitmappings and recursive mappings
	if ($read_id eq $read_id_split) {
	    my $remaining;
	    if ($splitfh && ($remaining=<$splitfh>)) {
		$lastsplit=$remaining;
		%lastsplit=%{parse_gem_line($lastsplit)};
	    }
	}

	if ($read_id eq $read_id_rec) {
	    my $remaining;
	    if ($recfh && ($remaining=<$recfh>)) {
		$lastrec=$remaining;
		%lastrec=%{parse_gem_line($lastrec)};
	    }
	}
	return($record);
    };
    return($bestmap_sub);
}


sub select_pair {
    my $record1=shift;
    my $record2=shift;
   
    my $type1=classify_pair($record1);
    my $type2=classify_pair($record2);

    my $type='unique';
    if ($type1 ne $type) {
	$type='uncertain';
    } elsif (($type1 eq $type2) &&
	     ($type1 eq 'multi')) {
	$type='multi';
    }
    return($type);
}

sub classify_pair {
    my $record=shift;
    my %record=%{parse_gem_line($record)};

    my $type='unique';
    unless ($record{'matches'}=~/^(0:)*1[^\d].*$/) {
	$type='multi';
    }
    return($type);
}
