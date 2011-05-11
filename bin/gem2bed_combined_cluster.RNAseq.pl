#!/soft/bin/perl
# DGK 2009

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# This script should parse the output of the GEM splitmapper at the moment it
# returns only the unique maps
### TO DO:
# Add the recursive mapping to the script 

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh');
use RNAseq_GEM3 ('parse_gem_line','gem_2_coords','gem_2_coords_split',
		 'coords2gtf','get_split_coord_parser',
		 'coords2bedSimple','coords2bedJunctions',
		 'get_junction_coords','line_reader');

# Get command line options
#my $read_length;
my $file_list;
my $localdir;

my %options=%{read_config_file()};
$file_list=$options{'FILELIST'};
$localdir=$options{'LOCALDIR'};

# Get some subroutines
my $dbh=get_dbh();
my $gen_map_table=$options{'PREFIX'}.'_genome_mapping';
my $junc_map_table=$options{'PREFIX'}.'_junctions_mapping';
my $split_map_table=$options{'PREFIX'}.'_split_mapping';
my $rec_map_table=$options{'PREFIX'}.'_recursive_mapping';

# Get some subroutines
*get_gen_file=get_files_from_table_sub($dbh,
				       $gen_map_table);
*get_junc_file=get_files_from_table_sub($dbh,
					$junc_map_table);
*get_split_file=get_files_from_table_sub($dbh,
					 $split_map_table);
*get_rec_file=get_files_from_table_sub($dbh,
				       $rec_map_table);
my $coord_parser=get_split_coord_parser();

my %files=%{read_file_list($file_list)};

# Get for each lane the files corresponding to the junctions mapping
my %lane_files=%{get_lane_files(\%files)};

# Process the files
print STDERR "Processing files\n";
my @map_files;
foreach my $pair (keys %lane_files) {
    my $bedfn=$localdir.'/'.$pair.'.combined.bed';
    if (-r $bedfn) {
	print STDERR $bedfn,"\tIs already present. Skipping...\n";
	next;
    }
    print STDERR "Building $bedfn...\n";

    my $bedfh=get_fh($bedfn,1);
    foreach my $lane (keys %{$lane_files{$pair}}) {
	# parse the genome mapping
	my $genome_map=$lane_files{$pair}{$lane}->[0];
	print STDERR "Parsing genome reads from $genome_map..";
	my $genfh=get_fh($genome_map);
	while (my $line=<$genfh>) {
	    my %line=%{parse_gem_line($line)};
	    my @coords=@{gem_2_coords(\%line)};
	    
	    foreach my $coords (@coords) {
		my $bed=coords2bedSimple($coords,
					 $genome_map);
		print $bedfh $bed,"\n";
	    }
	}
	close($genfh);
	
	#  parse the junctions mapping
	my $junction_map=$lane_files{$pair}{$lane}->[1];
	print STDERR "Parsing junction maps from $junction_map...";
	my $juncfh=get_fh($junction_map);
	my $interchrom=0;
	while (my $line=<$juncfh>) {
	    my %line=%{parse_gem_line($line)};
	    my @coords=@{gem_2_coords(\%line)};
	    
	    # Skip cases with no hits
	    if ($line{'hits'} eq '-') {
		next;
	    }
	    
	    my $oldbed='';
	    foreach my $coords (@coords) {
		my $bed=coords2bedJunctions($coords,
					    $junction_map);
		
		if ($bed) {
		    if ($bed eq $oldbed) {
			next;
		    } elsif ($oldbed &&
			     ($bed ne $oldbed)) {
			print $bedfh $oldbed,"\n";
			$oldbed=$bed;
		    } else {
			$oldbed=$bed;
		    }
		}
		if ($oldbed) {
		    print $bedfh $oldbed,"\n";
		}
	    }
	}
	close($juncfh);

	# parse the split mapping This will take only those cses with a single
	# split map possibility (no where one half is unique and the other multi
	my $split_map=$lane_files{$pair}{$lane}->[2];
	print STDERR "Parsing split maps from $split_map...";
	my $splitfh=get_fh($split_map);
	while (my $line=<$splitfh>) {
	    my %line=%{parse_gem_line($line)};
	    my @coords=@{gem_2_coords_split(\%line,
					    $coord_parser,
					    $line{'length'})};
	    
	    # Skip cases with no hits
	    if ($line{'hits'} eq '-') {
		next;
	    }
	    
	    my $oldbed='';
	    foreach my $coords (@coords) {
		my $bed=coords2splitJunctions($coords);
		
		if ($bed) {
		    if ($bed eq $oldbed) {
			next;
		    } elsif ($oldbed &&
			     ($bed ne $oldbed)) {
			print $bedfh $oldbed,"\n";
			$oldbed=$bed;
			next;
		    } else {
			$oldbed=$bed;
			next;
		    }
		}
		if ($oldbed) {
		    print $bedfh $oldbed,"\n";
		    $oldbed='';
		}
	    }
	}
	close($splitfh);

	# parse the recursive mapping This will take only those cases with a
	# single
	# split map possibility (no where one half is unique and the other multi
	my $rec_map=$lane_files{$pair}{$lane}->[3];
	print STDERR "Parsing rec maps from $rec_map...";
	my $recfh=get_fh($rec_map);
	while (my $line=<$recfh>) {
	    my %line=%{parse_gem_line($line)};
	    my $oldbed='';
	    if ($line{'type'} eq 'S') {
		my @coords=@{gem_2_coords_split(\%line,
						$coord_parser,
						$line{'length'})};
	    
		# Skip cases with no hits
		if ($line{'hits'} eq '-') {
		    next;
		}
	    
		foreach my $coords (@coords) {
		    my $bed=coords2splitJunctions($coords);
		    
		    if ($bed) {
			if ($bed eq $oldbed) {
			    next;
			} elsif ($oldbed &&
				 ($bed ne $oldbed)) {
			    print $bedfh $oldbed,"\n";
			    $oldbed=$bed;
			    next;
			} else {
			    $oldbed=$bed;
			    next;
			}
		    }
		    if ($oldbed) {
			print $bedfh $oldbed,"\n";
			$oldbed='';
		    }
		}
	    } else {
		my @coords=@{gem_2_coords(\%line)};
		
		foreach my $coords (@coords) {
		    my $bed=coords2bedSimple($coords,
					     $genome_map);
		    print $bedfh $bed,"\n";
		}
	    }
	}
	close($recfh);

	print STDERR "done\n";
    }
    close($bedfh);
}

exit;

sub coords2splitJunctions {
    my $coords=shift;
    
    my $bed=get_split_bed_coords2($coords);
    if ($bed) {
	return($bed);
    } else {
	return();
    }
}

sub get_split_bed_coords2 {
    my $coords=shift;

    my ($rstart1,$rend1,$rstart2,$rend2,$chr1,$chr2);
    $rstart1=$coords->{'up_start'}->[0]->[0];
    $rend1=$coords->{'up_start'}->[0]->[1];
    $rstart2=$coords->{'down_end'}->[0]->[0];
    $rend2=$coords->{'down_end'}->[0]->[1];
    $chr1=$coords->{'up_chr'};
    $chr2=$coords->{'down_chr'};

    my $strand1=$coords->{'up_strand'};
    my $strand2=$coords->{'down_strand'};

    # Skip cases where the strands are not the same as they cannot be encoded
    # in bed format
    if ($chr1 ne $chr2) {
#	print STDERR "Skipping interchromosomal hit from ",$coords->{'id'},"\n";
	return();
    }

    if ($strand1 ne $strand2) {
#	print STDERR "Skipping interchromosomal hit from ",$coords->{'id'},"\n";
	return();
    }

    # These are cases where the order is incorrect for the splits, so they
    # cannot be encoded in the bed format
    if ($rstart1 > $rstart2) {
	return();
    }

    # Rebuild the mismatch list
    my $mismatchmap='';
    my $mismatchno=0;
    if ($coords->{'mismatch_list'} &&
	@{$coords->{'mismatch_list'}}) {
	$mismatchno=@{$coords->{'mismatch_list'}};
    }
    foreach my $mismatch (@{$coords->{'mismatch_list'}}) {
	$mismatchmap.=$mismatch->[0];
	$mismatchmap.=$mismatch->[1];
    }

    my ($size1,$size2)=(0,0);
    if ($rstart1) {
	$size1=$rend1 - $rstart1 + 1;
    }
    if ($rstart2) {
	$size2=$rend2 - $rstart2 + 1;
    }

    my $size=join(',',$size1,$size2);
    my $start=join(',',0,$rstart2 - $rstart1);

    if ($size1 + $size2 != $coords->{'length'}) {
	warn "coordinate_problem: added half lengths differ from total\n";
	print STDERR join("\t",
			  $coords->{'length'},
			  $size1,
			  $size2,
			  $chr1,
			  $rstart1,
			  $rend1,
			  $rstart2,
			  $rend2,
			  $coords->{'id'},
			  $strand1,
			  $strand2),"\n";
	print STDERR join("\t",
			  $chr1,
			  $rstart1 - 1,
			  $rend2,
			  $coords->{'id'},
			  0,
			  $strand1,
			  $rstart1 - 1,
			  $rend2,
			  '255,0,0',
			  2,
			  $size,
			  $start),"\n";
	print STDERR "Please press return to continue...";
	<STDIN>;
    }

    if (($size1 > 0) &&
	($size2 > 0)) {
	my $bed=join("\t",
		     $chr1,
		     $rstart1 - 1,
		     $rend2,
		     $coords->{'id'},
		     0,
		     $strand1,
		     $rstart1 - 1,
		     $rend2,
		     '255,0,0',
		     2,
		     $size,
		     $start);

	# Run a check to make sure the bed format is correct
	my @sizes=split(',',$size);
	my @starts=split(',',$start);
	if ($rend2 != ($rstart1 + $sizes[-1] + $starts[-1] - 1)) {
	    die "Incorrect bed format\n$bed\n";
	}

	return($bed);
    } else {
	print STDERR "This should never happen:\t";
	die join("\t",
		 $chr1,$rstart1,$rend1,$strand1,$size1,
		 $chr2,$rstart2,$rend2,$strand2,$size2),"\n";
	return();
    }
}

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $junctionfile=get_junc_file($lane);
	my $genomefile=get_gen_file($lane);
	my $splitfile=get_split_file($lane);
	my $recfile=get_rec_file($lane);

	my $genomefilepath=$options{'GENOMEDIR'}.'/'.$genomefile;
	my $junctionfilepath=$options{'JUNCTIONSDIR'}.'/'.$junctionfile;
	my $splitfilepath=$options{'SPLITMAPDIR'}.'/'.$splitfile;
	my $recfilepath=$options{'RECMAPDIR'}.'/'.$recfile;

	# Get the genome mapping
	if (-r $genomefilepath) {
	    $lane_files{$pair}{$lane}->[0]=$genomefilepath;
	} elsif (-r $genomefilepath.'.gz') {
	    print STDERR "$genomefile is gzipped\n";
	    $lane_files{$pair}{$lane}->[0]=$genomefilepath.'.gz';
	} else {
	    die "Can't find $genomefile or $genomefile.gz\n";
	}

	# Get the junctions mapping
	if (-r $junctionfilepath) {
	    $lane_files{$pair}{$lane}->[1]=$junctionfilepath;
	} elsif (-r $junctionfilepath.'.gz') {
	    print STDERR "$junctionfile is gzipped\n";
	    $lane_files{$pair}{$lane}->[1]=$junctionfilepath.'.gz';
	} else {
	    die "Can't find $junctionfile or $junctionfile.gz\n";
	}

	# Get the split mapping
	if (-r $splitfilepath) {
	    $lane_files{$pair}{$lane}->[2]=$splitfilepath;
	} elsif (-r $splitfilepath.'.gz') {
	    print STDERR "$splitfile is gzipped\n";
	    $lane_files{$pair}{$lane}->[2]=$splitfilepath.'.gz';
	} else {
	    die "Can't find $splitfile or $splitfile.gz\n";
	}

	# Get the rec mapping
	if (-r $recfilepath) {
	    $lane_files{$pair}{$lane}->[3]=$recfilepath;
	} elsif (-r $recfilepath.'.gz') {
	    print STDERR "$recfile is gzipped\n";
	    $lane_files{$pair}{$lane}->[3]=$recfilepath.'.gz';
	} else {
	    die "Can't find $splitfile or $recfile.gz\n";
	}

    }

    return(\%lane_files);
}

sub get_closest_hit {
    my $read1=shift;
    my $read2=shift;

    my $ref_coord=int(($read1->[1] + $read1->[2]) / 2);

    my $closest;
    my $dist= 1e10;
    foreach my $read (@{$read2}) {
	if ($read1->[0] ne $read->[0]) {
	    # Exclude interchromosomal
	    next;
	} elsif ($read1->[3] eq $read->[3]) {
	    # Both reads have the same strand, which should not happed in
	    # paired unless there is something strange going on
	    next;
	} elsif (($read->[3] eq '+' && $read->[2] > $read1->[1]) ||
		 ($read->[3] eq '-' && $read->[1] < $read1->[2])) {
	    # Check if the orientation is correct
	    next;
	}
	my $prob_coord=int(($read->[1] + $read->[2]) / 2);
	my $prob_dist=abs($ref_coord - $prob_coord);
	if ($prob_dist < $dist) {
	    $dist=$prob_dist;
	    $closest=$read;
	}
    }
    return($closest);
}


