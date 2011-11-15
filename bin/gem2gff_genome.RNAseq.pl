#!/soft/bin/perl
# DGK

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

#    Author : David Gonzalez, david.gonzalez@crg.eu

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
# This script should parse the output of the GEM mapper
# It will take the read.list file and from there decide if the mappings we are
# looking at are paired or not.

# After this for each genome mapped file it will take the corresponding junction
# mapped file and filter out those possible pseudogene hits (those cases where
# a read hits both genome and junctions)

# It will print out the hits in gtf format separating them into two files, the
# unique mappings and the multimappings. The unique mappings will include those
# that have been resolved thanks to PE information

use Getopt::Long;
use RNAseq_pipeline3 ('get_fh','get_files_from_table_sub','print_gff');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh');
use RNAseq_GEM3 qw(process_gem_hit parse_gem_line line_reader get_closest_hit);

# Declare variables & Get command line options
my $mismatches;
my $closest_hit=0;
my $mismatch_strict;
my $unique_strict;
my $file_list;
my $genomedir;
my $stranded;

# Read command line options
GetOptions(
    'closest=i' => \$closest_hit,
    'mismatch_strict' => \$mismatch_strict,
    'unique_strict' => \$unique_strict
    );

if ($unique_strict) {
    print STDERR "WARNING: Unsing strict unique.\n";
    print STDERR "\tNo recovery of multi hits will be allowed\n";
}

my %options=%{read_config_file()};
$mismatches=$options{'MISMATCHES'};
$file_list=$options{'FILELIST'};
$genomedir=$options{'GENOMEDIR'};
$stranded=$options{'STRANDED'};

# Get a database connection
my $dbh=get_dbh();

my $junctions_map_table=$options{'PREFIX'}.'_junctions_mapping';
my $genome_map_table=$options{'PREFIX'}.'_genome_mapping';

# First ge the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get some subroutines
*get_junction_file=get_files_from_table_sub($dbh,
					    $junctions_map_table);
*get_genome_file=get_files_from_table_sub($dbh,
					  $genome_map_table);

# Get for each lane the files corresponding to the genome mapping and the
# junctions mapping
my %lane_files=%{get_lane_files(\%files)};

# Define the parsers
my %parsers=%{get_genome_parser_subs($mismatches,
				     $closest_hit,
				     $mismatch_strict,
				     $unique_strict,
				     $genomedir)};

# Process the files
foreach my $pair (keys %lane_files) {
    print STDERR "Processing $pair\n";
    if ($stranded) {
	print STDERR "Set identified as stranded\n";
    }
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

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $junctionfile=get_junction_file($lane);
	$junctionfile=~s/.map(.gz)?$/.map.gen.coords/;

	my $genomefile=get_genome_file($lane);

	my $genomefilepath=$options{'GENOMEDIR'}.'/'.$genomefile;
	my $junctionfilepath=$options{'JUNCTIONSDIR'}.'/'.$junctionfile;
	if (-r $genomefilepath) {
	    $lane_files{$pair}{$lane}->[0]=$genomefilepath;
	} elsif (-r $genomefilepath.'.gz') {
	    print STDERR "$genomefile is gzipped\n";
	    $lane_files{$pair}{$lane}->[0]=$genomefilepath.'.gz';
	} else {
	    die "Can't find $genomefile or $genomefile.gz\n";
	}

	if (-r $junctionfilepath) {
	    $lane_files{$pair}{$lane}->[1]=$junctionfilepath;
	} elsif (-r $junctionfilepath.'.gz') {
	    print STDERR "$junctionfile is gzipped\n";
	    $lane_files{$pair}{$lane}->[1]=$junctionfilepath.'.gz';
	} else {
	    die "Can't find $junctionfile or $junctionfile.gz\n";
	}

    }

    return(\%lane_files);
}

sub get_genome_parser_subs {
    my $mismatches=shift;
    my $closest=shift;
    my $mismatch_strict=shift;
    my $unique_strict=shift;
    my $genomedir=shift;
    my $multi;
    my %parser;
    my $qualities;
    my %stats;

    $parser{'single'}=sub {
	my $set=shift;
	my $pair=shift;

	my @paired=keys %{$set};			     
	my $infile=$set->{$paired[0]}->[0];
	my $juncfile=$set->{$paired[0]}->[1];

	my $type='single';

	unless($infile && $juncfile) {
	    die "A genome and a  junction file are required for pseudogene filtering\n";
	}

	my $outfileunique=$genomedir.'/'.$pair.'.single.unique.gtf.gz';
	my $outfilemulti=$genomedir.'/'.$pair.'.single.multi.gtf.gz';
	my $outfilepseudo=$genomedir.'/'.$pair.'.single.pseudo.gtf.gz';

	print STDERR "Printing results to $pair.single.[unique|multi|pseudo].gtf.gz\n";

	my %stats;

	# Open the required files
	# Output files
	my $outuniquefh=get_fh($outfileunique,1);
	my $outmultifh=get_fh($outfilemulti,1);
	my $outpseudofh=get_fh($outfilepseudo,1);

	# Inputfiles
	# check if the infile is there or if it is zipped
	if (-r $infile) {
	    print STDERR $infile,"\tPresent\n";
	} elsif (-r $infile.'.gz') {
	    $infile.='.gz';
	    print STDERR $infile,"\tPresent\n";
	} else {
	    die "I can't find $infile\n";
	}

	# check if the junctionfile is there or if it is zipped
	if (-r $juncfile) {
	    print STDERR $juncfile,"\tPresent\n";
	} elsif (-r $juncfile.'.gz') {
	    $juncfile.='.gz';
	    print STDERR $juncfile,"\tPresent\n";
	} else {
	    die "I can't find $juncfile\n";
	}

	my $readsfh=get_fh($infile);
	my $juncfh=get_fh($juncfile);
	
	my $reader=line_reader($readsfh,
			       $mismatches,
			       $closest,
			       $mismatch_strict,
			       $unique_strict);
	my $juncreader=line_reader($juncfh,
				   $mismatches,
				   $closest,
				   $mismatch_strict,
				   $unique_strict);

	print STDERR 'Filtering...';
	while (my $read=&$reader) {
	    my $junc=&$juncreader;

	    unless ($junc) {
		die "ERROR:Insuficient lines in the junctions mapping\n";
	    }
	    # Check if the files are sorted correctly
	    if ($read->{'id'} ne $junc->{'id'}) {
		print STDERR $read->{'id'},' paired with ', $junc->{'id'},"\n";
		print STDERR "Please check sorting of $infile and juncfile\n";
		die "Sorting of files is inconsistent\n";
	    }


	    # Extract junction mappings that also map to the genome uniquely
	    # in a separate file
	    $stats{'total'}++;
	    if ((@{$read->{'matches'}} == 1) &&
		(@{$junc->{'matches'}})) {
		$stats{'read1_pseudo'}++;
		    print_gff($outpseudofh,
			      $infile,
			      $read->{'id'},
			      $read->{'matches'}->[0],
			      'pseudo');
		# Do not include these in the final mappings
		next;
	    }
	    # Filter the reads for maps
	    if (@{$read->{'matches'}}) {
		if (@{$read->{'matches'}} == 1) {
		    my $match=$read->{'matches'}->[0];
		    if ($stranded && $read->{'id'}=~/1$/) {
			if ($match->[3] eq '-') {
			    $match->[3]='+';
			} elsif ($match->[3] eq '+') {
			    $match->[3]='-';
			} else {
			    warn "Unknown strand\n";
			    $match->[3]='.';
			}
		    }
		    $stats{'unique'}++;
		    print_gff($outuniquefh,
			      $infile,
			      $read->{'id'},
			      $match,
			      $type);
		} else {
		    $stats{'multi'}++;
		    foreach my $match (@{$read->{'matches'}}) {
			# If the reads are stranded change the strand of the
			# read 1
			if ($stranded && $read->{'id'}=~/1$/) {
			    if ($match->[3] eq '-') {
				$match->[3]='+';
			    } elsif ($match->[3] eq '+') {
				$match->[3]='-';
			    } else {
				warn "Unknown strand\n";
				$match->[3]='.';
			    }
			}
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

    $parser{'paired'}=sub {
	my $set=shift;
	my $pair=shift;

	my @paired=keys %{$set};			     
	my $infile1=$set->{$paired[0]}->[0];
	my $infile2=$set->{$paired[1]}->[0];
	my $juncfile1=$set->{$paired[0]}->[1];
	my $juncfile2=$set->{$paired[1]}->[1];

	my $type='paired';

	unless($infile1 && $infile2) {
	    die "Two files required for paired ends\n";
	}

	unless($juncfile1 && $juncfile2) {
	    die "Two junction files required for pseudogene filtering\n";
	}

	my $outfileunique=$genomedir.'/'.$pair.'.paired.unique.gtf.gz';
	my $outfilemulti=$genomedir.'/'.$pair.'.paired.multi.gtf.gz';
	my $outfilepseudo=$genomedir.'/'.$pair.'.paired.pseudo.gtf.gz';

	print STDERR "Printing results to $pair.paired.[unique|multi|pseudo].gtf.gz\n";
	my %stats;
	my %distances;
	
	# Open the required files
	# Output files
	my $outuniquefh=get_fh($outfileunique,1);
	my $outmultifh=get_fh($outfilemulti,1);
	my $outpseudofh=get_fh($outfilepseudo,1);

	# Inputfiles
	# check if the infiles are there or if it is zipped
	if (-r $infile1) {
	    print STDERR $infile1,"\tPresent\n";
	} elsif (-r $infile1.'.gz') {
	    $infile1.='.gz';
	    print STDERR $infile1,"\tPresent\n";
	} else {
	    die "I can't find $infile1\n";
	}

	if (-r $infile2) {
	    print STDERR $infile2,"\tPresent\n";
	} elsif (-r $infile2.'.gz') {
	    $infile2.='.gz';
	    print STDERR $infile2,"\tPresent\n";
	} else {
	    die "I can't find $infile2\n";
	}

	# check if the junctionfile is there or if it is zipped
	if (-r $juncfile1) {
	    print STDERR $juncfile1,"\tPresent\n";
	} elsif (-r $juncfile1.'.gz') {
	    $juncfile1.='.gz';
	    print STDERR $juncfile1,"\tPresent\n";
	} else {
	    die "I can't find $juncfile1\n";
	}

	# check if the junctionfile is there or if it is zipped
	if (-r $juncfile2) {
	    print STDERR $juncfile2,"\tPresent\n";
	} elsif (-r $juncfile2.'.gz') {
	    $juncfile2.='.gz';
	    print STDERR $juncfile2,"\tPresent\n";
	} else {
	    die "I can't find $juncfile2\n";
	}

	my $reads1fh=get_fh($infile1);
	my $reads2fh=get_fh($infile2);

	my $junc1fh=get_fh($juncfile1);
	my $junc2fh=get_fh($juncfile2);
	
	my $reader1=line_reader($reads1fh,
				$mismatches,
				$closest,
				$mismatch_strict,
				$unique_strict);
	my $reader2=line_reader($reads2fh,
				$mismatches,
				$closest,
				$mismatch_strict,
				$unique_strict);
	my $juncreader1=line_reader($junc1fh,
				    $mismatches,
				    $closest,
				    $mismatch_strict,
				    $unique_strict);
	my $juncreader2=line_reader($junc2fh,
				    $mismatches,
				    $closest,
				    $mismatch_strict,
				    $unique_strict);

	print STDERR 'Filtering, mixing & matching...';
	while ((my $read1=&$reader1) &&
	       (my $read2=&$reader2)){
	    my $junc1=&$juncreader1;
	    my $junc2=&$juncreader2;
	    unless ($junc1 && $junc2) {
		die "ERROR:Insuficient lines in the junctions mapping\n";
	    }

	    # Check if the files are sorted correctly
	    if ($read1->{'id'} ne $junc1->{'id'}) {
		print STDERR $read1->{'id'},' paired with ',$junc1->{'id'},"\n";
		die "Sorting of files is inconsistent\n";
	    }
	    if ($read2->{'id'} ne $junc2->{'id'}) {
		print STDERR $read2->{'id'},' paired with ',$junc2->{'id'},"\n";
		die "Sorting of files is inconsistent\n";
	    }

	    # Extract junction mappings that also map to the genome uniquely
	    # in a separate file
	    $stats{'total'}++;
	    if (((@{$read1->{'matches'}} == 1) &&
		(@{$junc1->{'matches'}})) &&
		((@{$read2->{'matches'}} == 1) &&
		 (@{$junc2->{'matches'}}))) {
		$stats{'both_pseudo'}++;
		my $match1=$read1->{'matches'}->[0];
		my $match2=$read2->{'matches'}->[0];
		if ($stranded) {
		    if ($read1->{'id'}=~/1$/) {
			if ($match1->[3] eq '-') {
			    $match1->[3]='+';
			} elsif ($match1->[3] eq '+') {
			    $match1->[3]='-';
			} else {
			    warn "Unknown strand\n";
			    $match1->[3]='.';
			}
		    } elsif ($read2->{'id'}=~/1$/) {
			if ($match2->[3] eq '-') {
			    $match2->[3]='+';
			} elsif ($match2->[3] eq '+') {
			    $match2->[3]='-';
			} else {
			    warn "Unknown strand\n";
			    $match2->[3]='.';
			}
		    }
		}
		print_gff($outpseudofh,
			  $infile1,
			  $read1->{'id'},
			  $match1,
			  'pseudo');
		print_gff($outpseudofh,
			  $infile2,
			  $read2->{'id'},
			  $match2,
			  'pseudo');
		next;
	    } else {
		if ((@{$read1->{'matches'}} == 1) &&
		    (@{$junc1->{'matches'}})) {
		    $stats{'read1_pseudo'}++;
		    my $match1=$read1->{'matches'}->[0];
		    if ($stranded) {
			if ($read1->{'id'}=~/1$/) {
			    if ($match1->[3] eq '-') {
				$match1->[3]='+';
			    } elsif ($match1->[3] eq '+') {
				$match1->[3]='-';
			    } else {
				warn "Unknown strand\n";
				$match1->[3]='.';
			    }
			}			
		    }
		    print_gff($outpseudofh,
			      $infile1,
			      $read1->{'id'},
			      $match1,
			      'pseudo');
		    next;
		}
		if ((@{$read2->{'matches'}} == 1) &&
		    (@{$junc2->{'matches'}})) {
		    $stats{'read2_pseudo'}++;
		    my $match2=$read2->{'matches'}->[0];
		    if ($stranded) {
			if ($read2->{'id'}=~/1$/) {
			    if ($match2->[3] eq '-') {
				$match2->[3]='+';
			    } elsif ($match2->[3] eq '+') {
				$match2->[3]='-';
			    } else {
				warn "Unknown strand\n";
				$match2->[3]='.';
			    }
			}			
		    }
		    print_gff($outpseudofh,
			      $infile2,
			      $read2->{'id'},
			      $match2,
			      'pseudo');
		    next;
		}
	    }
	    # Filter the reads for maps
	    if ((@{$read1->{'matches'}}) &&
		(@{$read2->{'matches'}})) {
		# We have matches for both ends
		if ((@{$read1->{'matches'}} == 1) && 
		    (@{$read2->{'matches'}} == 1)) {
		    $stats{'double_unique'}++;
		    my $match1=$read1->{'matches'}->[0];
		    my $match2=$read2->{'matches'}->[0];
		    if ($stranded) {
			if ($read1->{'id'}=~/1$/) {
			    if ($match1->[3] eq '-') {
				$match1->[3]='+';
			    } elsif ($match1->[3] eq '+') {
				$match1->[3]='-';
			    } else {
				warn "Unknown strand\n";
				$match1->[3]='.';
			    }
			} elsif ($read2->{'id'}=~/1$/) {
			    if ($match2->[3] eq '-') {
				$match2->[3]='+';
			    } elsif ($match2->[3] eq '+') {
				$match2->[3]='-';
			    } else {
				warn "Unknown strand\n";
			    $match2->[3]='.';
			    }
			}
		    }
		    print_gff($outuniquefh,
			      $infile1,
			      $read1->{'id'},
			      $match1,
			      $type);
		    print_gff($outuniquefh,
			      $infile2,
			      $read2->{'id'},
			      $match2,
			      $type);
		} elsif (@{$read1->{'matches'}} == 1) {
		    $stats{'read1_unique'}++;
		    # select the closest hit
		    my $closest=get_closest_hit($read1->{'matches'}->[0],
						$read2->{'matches'},
						$unique_strict,
						$stranded);
		    if ($closest) {
			$stats{'recovered_read2'}++;
			my $match1=$read1->{'matches'}->[0];
			my $match2=$closest;
			if ($stranded) {
			    if ($read1->{'id'}=~/1$/) {
				if ($match1->[3] eq '-') {
				    $match1->[3]='+';
				} elsif ($match1->[3] eq '+') {
				    $match1->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match1->[3]='.';
				}
			    } elsif ($read2->{'id'}=~/1$/) {
				if ($match2->[3] eq '-') {
				    $match2->[3]='+';
				} elsif ($match2->[3] eq '+') {
				    $match2->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match2->[3]='.';
				}
			    }
			}
			print_gff($outuniquefh,
				  $infile1,
				  $read1->{'id'},
				  $match1,
				  $type);
			print_gff($outuniquefh,
				  $infile2,
				  $read2->{'id'},
				  $match2,
				  $type);
		    } else {
			my $match1=$read1->{'matches'}->[0];
			if ($stranded) {
			    if ($read1->{'id'}=~/1$/) {
				if ($match1->[3] eq '-') {
				    $match1->[3]='+';
				} elsif ($match1->[3] eq '+') {
				    $match1->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match1->[3]='.';
				}
			    }
			}
			print_gff($outuniquefh,
				  $infile1,
				  $read1->{'id'},
				  $match1,
				  'unmatched');
			foreach my $match (@{$read2->{'matches'}}) {
			    if ($stranded) {
				if ($read2->{'id'}=~/1$/) {
				    if ($match->[3] eq '-') {
					$match->[3]='+';
				    } elsif ($match->[3] eq '+') {
					$match->[3]='-';
				    } else {
					warn "Unknown strand\n";
					$match->[3]='.';
				    }
				}
			    }
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
						$read1->{'matches'},
						$unique_strict,
						$stranded);
		    if ($closest) {
			my $match1=$closest;
			my $match2=$read2->{'matches'}->[0];
			if ($stranded) {
			    if ($read1->{'id'}=~/1$/) {
				if ($match1->[3] eq '-') {
				    $match1->[3]='+';
				} elsif ($match1->[3] eq '+') {
				    $match1->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match1->[3]='.';
				}
			    } elsif ($read2->{'id'}=~/1$/) {
				if ($match2->[3] eq '-') {
				    $match2->[3]='+';
				} elsif ($match2->[3] eq '+') {
				    $match2->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match2->[3]='.';
				}
			    }
			}
			$stats{'recovered_read1'}++;
			print_gff($outuniquefh,
				  $infile1,
				  $read1->{'id'},
				  $match1,
				  $type);
			print_gff($outuniquefh,
				  $infile2,
				  $read2->{'id'},
				  $match2,
				  $type);
		    } else {
			my $match2=$read2->{'matches'}->[0];
			if ($stranded) {
			    if ($read2->{'id'}=~/1$/) {
				if ($match2->[3] eq '-') {
				    $match2->[3]='+';
				} elsif ($match2->[3] eq '+') {
				    $match2->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match2->[3]='.';
				}
			    }
			}
			print_gff($outuniquefh,
				  $infile2,
				  $read2->{'id'},
				  $match2,
				  'unmatched');
			foreach my $match (@{$read1->{'matches'}}) {
			    if ($stranded) {
				if ($read2->{'id'}=~/1$/) {
				    if ($match->[3] eq '-') {
					$match->[3]='+';
				    } elsif ($match->[3] eq '+') {
					$match->[3]='-';
				    } else {
					warn "Unknown strand\n";
					$match->[3]='.';
				    }
				}
			    }
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
			if ($stranded) {
			    if ($read1->{'id'}=~/1$/) {
				if ($match->[3] eq '-') {
				    $match->[3]='+';
				} elsif ($match->[3] eq '+') {
				    $match->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match->[3]='.';
				}
			    }
			}
			print_gff($outmultifh,
				  $infile1,
				  $read1->{'id'},
				  $match,
				  'multi');
		    }
		    foreach my $match (@{$read2->{'matches'}}) {
			if ($stranded) {
			    if ($read2->{'id'}=~/1$/) {
				if ($match->[3] eq '-') {
				    $match->[3]='+';
				} elsif ($match->[3] eq '+') {
				    $match->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match->[3]='.';
				}
			    }
			}
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
		    my $match1=$read1->{'matches'}->[0];
		    if ($stranded) {
			if ($read1->{'id'}=~/1$/) {
			    if ($match1->[3] eq '-') {
				$match1->[3]='+';
			    } elsif ($match1->[3] eq '+') {
				$match1->[3]='-';
			    } else {
				warn "Unknown strand\n";
				$match1->[3]='.';
			    }
			}
		    }
		    print_gff($outuniquefh,
			      $infile1,
			      $read1->{'id'},
			      $match1,
			      'unmatched');
		} else {
		    $stats{'read1_multi'}++;
		    foreach my $match (@{$read1->{'matches'}}) {
			if ($stranded) {
			    if ($read1->{'id'}=~/1$/) {
				if ($match->[3] eq '-') {
				    $match->[3]='+';
				} elsif ($match->[3] eq '+') {
				    $match->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match->[3]='.';
				}
			    }
			}
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
		    my $match2=$read2->{'matches'}->[0];
		    if ($stranded) {
			if ($read2->{'id'}=~/1$/) {
			    if ($match2->[3] eq '-') {
				$match2->[3]='+';
			    } elsif ($match2->[3] eq '+') {
				$match2->[3]='-';
			    } else {
				warn "Unknown strand\n";
				$match2->[3]='.';
			    }
			}
		    }
		    print_gff($outuniquefh,
			      $infile2,
			      $read2->{'id'},
			      $match2,
			      'unmatched');
		} else {
		    $stats{'read2_multi'}++;
		    foreach my $match (@{$read2->{'matches'}}) {
			if ($stranded) {
			    if ($read2->{'id'}=~/1$/) {
				if ($match->[3] eq '-') {
				    $match->[3]='+';
				} elsif ($match->[3] eq '+') {
				    $match->[3]='-';
				} else {
				    warn "Unknown strand\n";
				    $match->[3]='.';
				}
			    }
			}
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
	close($junc1fh);
	close($junc2fh);
	close($outuniquefh);
	close($outmultifh);
	close($outpseudofh);
	
	foreach my $key (keys %stats) {
	    print STDERR join("\t",
			      $key,
			      $stats{$key}),"\n";
	}
    };
    
    return(\%parser);
}
