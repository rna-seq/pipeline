#!/soft/bin/perl

use strict;
use warnings;

BEGIN {
    unshift @INC, '/users/rg/dgonzalez/lib/Perl';
}

# This script is based on the awk script of the same name but with increased
# functionality as it combines gem2gff_unique and gem2gff_veryunique

use Getopt::Long;
use RNAseq_GEM3 qw(parse_gem_line process_gem_hit line_reader get_closest_hit);
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
my $trans_map_table=$options{'PREFIX'}.'_transcriptome_mapping';
*get_transcript_file=get_files_from_table_sub($dbh,
					      $trans_map_table);
#*unique_parser=get_unique_parser($mismatches);

# Get for each lane the files corresponding to the genome mapping and the
# junctions mapping
my %lane_files=%{get_lane_files(\%files)};

# Define the parsers
my %parsers=%{get_transcriptome_parser_subs($mismatches,
					    $closest_hit,
					    $mismatch_strict)};

# Process the files
foreach my $pair (keys %lane_files) {
    print STDERR "Processing $pair\n";
    process_files($lane_files{$pair},
		  $pair);
}

exit;

sub process_files {
    my $set=shift;
    my $pair=shift;

    my @paired=keys %{$set};

    if (@paired == 1) {
	print STDERR "Set identified as single reads\n";
	$parsers{'single'}->($set,
			     $paired[0]);
    } elsif (@paired == 2) {
	print STDERR "Set identified as paired reads\n";
	$parsers{'paired'}->($set,
			     $pair);
    } else {
	# This should never happen
	die "Apparently there are more than two files grouped\n";
    }
}

sub get_unique_parser {
    my $mismatches=shift;
    my $multi;
    my $parser;
    my $qualities;
    my %stats;

    $parser=sub {
	my $infile=shift;
	my $lane=shift;

	my $type='transcript';

	unless($infile) {
	    die "A transcriptome mapping file is required\n";
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
	
	my $reader=line_reader($readsfh,
			       $mismatches,
			       $closest_hit);

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
		    $stats{'multi'}++;
		    foreach my $match (@{$read->{'matches'}}) {
			print_gff($outmultifh,
				  $infile,
				  $read->{'id'},
				  $match,
				  'multi');
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

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $transfile=get_transcript_file($lane);
	my $filepath=$options{'TRANSDIR'}.'/'.$transfile;
	if (-r $filepath) {
	    $lane_files{$pair}{$lane}=$filepath;
	} elsif (-r $filepath.'.gz') {
	    print STDERR "$transfile is gzipped\n";
	    $lane_files{$pair}{$lane}=$filepath.'.gz';
	} else {
	    die "Can't find $transfile or $transfile.gz\n";
	}
    }

    return(\%lane_files);
}

# Get the actual gem parsers
sub get_transcriptome_parser_subs {
    my $mismatches=shift;
    my $closest=shift;
    my $mismatch_strict=shift;
    my $multi;
    my %parser;
    my $qualities;
    my %stats;

    $parser{'single'}=sub {
	my $set=shift;
	my $pair=shift;

	my @paired=keys %{$set};			     
	my $infile=$set->{$paired[0]};

	my $type='single';

	my $outfileunique=$pair.'.single.unique.gtf.gz';
	my $outfilemulti=$pair.'.single.multi.gtf.gz';

	print STDERR "Printing results to $pair.single.[unique|multi].gtf.gz\n";

	my %stats;

	# Open the required files
	# Output files
	my $outuniquefh=get_fh($outfileunique,1);
	my $outmultifh=get_fh($outfilemulti,1);
	
	# Inputfiles
	my $readsfh=get_fh($infile);
	
	my $reader=line_reader($readsfh,
			       $mismatches,
			       $closest,
			       $mismatch_strict);
	
	print STDERR 'Filtering...';
	while (my $read=&$reader) {
	    $stats{'total'}++;
	    
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
		    $stats{'multi'}++;
		    foreach my $match (@{$read->{'matches'}}) {
			print_gff($outmultifh,
				  $infile,
				  $read->{'id'},
				  $match,
				  'multi');
		    }
		}
	    } else {
		# We have no match
		$stats{'no_hit'}++;
	    }
	}
	close($readsfh);
	print STDERR "done\n";
    };

    $parser{'paired'}=sub {
	my $set=shift;
	my $pair=shift;

	my @paired=keys %{$set};			     
	my $infile1=$set->{$paired[0]};
	my $infile2=$set->{$paired[1]};

	my $type='paired';

	unless($infile1 && $infile2) {
	    die "Two files required for paired ends\n";
	}

	my $outfileunique=$pair.'.paired.unique.gtf.gz';
	my $outfilemulti=$pair.'.paired.multi.gtf.gz';

	print STDERR "Printing results to $pair.paired.[unique|multi].gtf.gz\n";
	my %stats;
	my %distances;
	
	# Open the required files
	# Output files
	my $outuniquefh=get_fh($outfileunique,1);
	my $outmultifh=get_fh($outfilemulti,1);

	# Inputfiles
	my $reads1fh=get_fh($infile1);
	my $reads2fh=get_fh($infile2);

	my $reader1=line_reader($reads1fh,
				$mismatches,
				$closest,
				$mismatch_strict);
	my $reader2=line_reader($reads2fh,
				$mismatches,
				$closest,
				$mismatch_strict);
	
	print STDERR 'Filtering...';
	while ((my $read1=&$reader1) &&
	       (my $read2=&$reader2)){
	    $stats{'total'}++;

	    # Filter the reads for maps
	    if ((@{$read1->{'matches'}}) &&
		(@{$read2->{'matches'}})) {
		# We have matches for both ends
		if ((@{$read1->{'matches'}} == 1) && 
		    (@{$read2->{'matches'}} == 1)) {
		    $stats{'double_unique'}++;
		    print_gff($outuniquefh,
			      $infile1,
			      $read1->{'id'},
			      $read1->{'matches'}->[0],
			      $type);
		    print_gff($outuniquefh,
			      $infile2,
			      $read2->{'id'},
			      $read2->{'matches'}->[0],
			      $type);
		} elsif (@{$read1->{'matches'}} == 1) {
		    $stats{'read1_unique'}++;
		    # select the closest hit
		    my $closest=get_closest_hit($read1->{'matches'}->[0],
						$read2->{'matches'});
		    if ($closest) {
			$stats{'recovered_read2'}++;
			print_gff($outuniquefh,
				  $infile1,
				  $read1->{'id'},
				  $read1->{'matches'}->[0],
				  $type);
			print_gff($outuniquefh,
				  $infile2,
				  $read2->{'id'},
				  $closest,
				  $type);
		    } else {
			print_gff($outuniquefh,
				  $infile1,
				  $read1->{'id'},
				  $read1->{'matches'}->[0],
				  'unmatched');
			foreach my $match (@{$read2->{'matches'}}) {
			    print_gff($outmultifh,
				      $infile2,
				      $read2->{'id'},
				      $match,
				      'multi');
			}
		    }
		} elsif (@{$read2->{'matches'}} == 1) {
		    $stats{'read2_unique'}++;
		    # select the closest hit
		    my $closest=get_closest_hit($read2->{'matches'}->[0],
						$read1->{'matches'});
		    if ($closest) {
			$stats{'recovered_read1'}++;
			print_gff($outuniquefh,
				  $infile1,
				  $read1->{'id'},
				  $closest,
				  $type);
			print_gff($outuniquefh,
				  $infile2,
				  $read2->{'id'},
				  $read2->{'matches'}->[0],
				  $type);
		    } else {
			print_gff($outuniquefh,
				  $infile2,
				  $read2->{'id'},
				  $read2->{'matches'}->[0],
				  'unmatched');
			foreach my $match (@{$read1->{'matches'}}) {
			    print_gff($outmultifh,
				      $infile1,
				      $read1->{'id'},
				      $match,
				      'multi');
			}
		    }
		} else {
		    $stats{'both_multi'}++;
		    foreach my $match (@{$read1->{'matches'}}) {
			print_gff($outmultifh,
				  $infile1,
				  $read1->{'id'},
				  $match,
				  'multi');
		    }
		    foreach my $match (@{$read2->{'matches'}}) {
			print_gff($outmultifh,
				  $infile2,
				  $read2->{'id'},
				  $match,
				  'multi');
		    }
		}
	    } elsif (@{$read1->{'matches'}}) {
		# We have matches only for the first end
		if (@{$read1->{'matches'}} == 1) {
		    $stats{'read1_unique'}++;
		    print_gff($outuniquefh,
			      $infile1,
			      $read1->{'id'},
			      $read1->{'matches'}->[0],
			      'unmatched');
		} else {
		    $stats{'read1_multi'}++;
		    foreach my $match (@{$read1->{'matches'}}) {
			print_gff($outmultifh,
				  $infile1,
				  $read1->{'id'},
				  $match,
				  'multi');
		    }
		}
		
	    } elsif (@{$read2->{'matches'}}) {
		# We have matches only for the second
		if (@{$read2->{'matches'}} == 1) {
		    $stats{'read2_unique'}++;
		    print_gff($outuniquefh,
			      $infile2,
			      $read2->{'id'},
			      $read2->{'matches'}->[0],
			      'unmatched');
		} else {
		    $stats{'read2_multi'}++;
		    foreach my $match (@{$read2->{'matches'}}) {
			print_gff($outmultifh,
				  $infile2,
				  $read2->{'id'},
				  $match,
				  'multi');
		    }
		}
	    } else {
		# We have no match
		$stats{'no_hit'}++;
	    }
	}
	print STDERR "done\n";
	
	# Close the files
	close($reads1fh);
	close($reads2fh);
	close($outuniquefh);
	close($outmultifh);
	
	foreach my $key (keys %stats) {
	    print STDERR join("\t",
			      $key,
			      $stats{$key}),"\n";
	}
    };
    
    return(\%parser);
}

