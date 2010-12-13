#!/soft/bin/perl

use strict;
use warnings;

BEGIN {
    unshift @INC, '/users/rg/dgonzalez/lib/Perl';
}

# This script will extrac the unique mappings to the juctions
# This script will report as uniquely mapping those reads that may map to
# several junctions, that all share the same junction position

use Getopt::Long;
use RNAseq_GEM3 qw(parse_gem_line process_gem_hit);
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub print_gff);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh');

# Declare some variables
my $mismatches;
my $file_list;
my $closest_hit=0;
my $mismatch_strict;

# Read command line options
GetOptions(
	   'closest=i' => \$closest_hit,
	   'mismatch_strict' => \$mismatch_strict,
    );

my %options=%{read_config_file()};
$mismatches=$options{'MISMATCHES'};
$file_list=$options{'FILELIST'};

# Get a database connection
my $dbh=get_dbh();

my %files=%{read_file_list($file_list)};

# Get some subroutines
my $junc_map_table=$options{'PREFIX'}.'_junctions_mapping';
*get_junctions_file=get_files_from_table_sub($dbh,
					     $junc_map_table);
*unique_parser=get_unique_parser($mismatches);

# Get for each lane the files corresponding to the junctions mapping
my %lane_files=%{get_lane_files(\%files)};

# Process the files
print STDERR "Processing files as single reads\n";
foreach my $pair (keys %lane_files) {
    foreach my $lane (keys %{$lane_files{$pair}}) {
	print STDERR "Processing $pair\n";
	unique_parser($lane_files{$pair}{$lane},
		      $lane);
    }
}

exit;

sub get_unique_parser {
    my $mismatches=shift;
    my $multi;
    my $parser;
    my $qualities;
    my %stats;

    $parser=sub {
	my $infile=shift;
	my $lane=shift;

	my $type='junction';

	unless($infile) {
	    die "A junctions mapping file is required\n";
	}

	my $outfileunique=$lane.'.single.unique.gtf.gz';
	my $outfilemulti=$lane.'.single.multi.gtf.gz';

	print STDERR "Printing results to $lane.single.[unique|multi].gtf.gz\n";

	my %stats;

	# Open the required files
	# Output files
	my $outuniquefh=get_fh($outfileunique,1);
	my $outmultifh=get_fh($outfilemulti,1);

	# Inputfiles
	my $readsfh=get_fh($infile);
	
	my $reader=line_reader($readsfh);

	print STDERR 'Filtering...';
	while (my $read=&$reader) {
	    # Filter the reads for maps
	    if (@{$read->{'matches'}}) {
		# We have matches for both ends
		if (@{$read->{'matches'}} == 1) {
		    $stats{'unique'}++;
		    print_gff($outuniquefh,
			      $infile,
			      $read->{'id'},
			      $read->{'matches'}->[0],
			      $type);
		} else {
		    my %multi;
		    foreach my $match (@{$read->{'matches'}}) {
			my @junction=split('_',$match->[0]);
			my $junctionpoint=join('_',@junction[2,3,5,6]);
			push @{$multi{$junctionpoint}},[$infile,
							$read->{'id'},
							$match];
		    }
		    my @positions=keys %multi;
		    if (@positions > 1) {
			$stats{'multi'}++;
			foreach my $key (@positions) {
			    foreach my $hit (@{$multi{$key}}) {
				print_gff($outmultifh,
					  @{$hit},
					  'multi');
			    }
			}
		    } else {
			$stats{'shared'}++;
			foreach my $hit (@{$multi{$positions[0]}}) {
			    print_gff($outuniquefh,
				      @{$hit},
				      'shared');
			}
		    }
		}
	    } else {
		# We have no match
		$stats{'no_hit'}++;
	    }
	}
	print STDERR "done\n";
    };

    return($parser);
}

# This is the main body of the parsers
sub line_reader {
    my $fh=shift;

    my $reader=sub {
	if (my $line=<$fh>) {
	    my %read=%{parse_gem_line($line)};
	    
	    my @map=split(':',$read{'matches'});

	    # Get best match position
	    my $best_match=0;
	    for (my $i=0;$i<=$mismatches;$i++) {
		if ($map[$i] eq '!') {
		    $best_match=-1;
		    last;
		} elsif ($map[$i] != 0) {
		    $best_match=$i;		    
		    last;
		} 
	    }

	    # Skip hits with too many matches
	    if ($best_match < 0) {
		next;
	    }

	    # From the @map array remove those cases that have more mismatches
	    # than the mismatch threshold
	    if ($mismatch_strict) {
		splice(@map,$mismatches + 1);
	    }

	    # From the @map array get the number of matches we have to take
	    # First remove the first part of the @map array up to the first
	    # match

	    splice(@map,0,$best_match);
	    # After remove the second part from the position after the first
	    # match plus the closest hit number we want
	    my $offset=$closest_hit + 1;
	    unless ($offset > @map) {
		splice(@map,$closest_hit + 1);
	    }
	    my $matches=0;
	    foreach my $hits (@map) {
		$matches+=$hits;
	    }

	    if ($matches < 1) {
		# If there are no matches skip the rest
		$read{'matches'}=[];
		return(\%read);
	    }

	    my @matches=split(',',$read{'hits'});
	    if ($matches > @matches) {
		if (@matches) {
		    $matches=@matches;
		} else {
		    $read{'matches'}=[];
		    return(\%read);
		}
	    }
	    splice(@matches,$matches);
	    my @processed;
	    foreach my $match (@matches) {
		my $processed=process_gem_hit($match,
					      $read{'length'},
					      $read{'matches'});
		if ($processed) {
		    push @processed, $processed;
		}
	    }
	    $read{'matches'}=[@processed];
	    return(\%read);
	} else {
	    return();
	}
    };
    return($reader);
}

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $juncfile=get_junctions_file($lane);
	my $filepath=$options{'JUNCTIONSDIR'}.'/'.$juncfile;
	if (-r $filepath) {
	    $lane_files{$pair}{$lane}=$filepath;
	} elsif (-r $filepath.'.gz') {
	    $lane_files{$pair}{$lane}=$filepath.'.gz';
	} else {
	    die "Can't find $filepath or $filepath.gz\n";
	}
    }

    return(\%lane_files);
}

