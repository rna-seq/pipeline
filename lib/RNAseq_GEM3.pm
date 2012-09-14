package RNAseq_GEM3;

#  GRAPE
#  Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# This package contains subroutines useful to parse GEM as well as output
# Different formats

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('parse_gem_line','gem_2_coords','gem_2_coords_split',
	    'coords2gtf','coords2bedSimple',
	    'coords2bedJunctions','check_index','check_input','get_fh',
	    'determine_quality_type','get_mapper_routines','get_best_hits',
	    'process_gem_hit','get_unmapped','check_split_input',
	    'get_junction_coords','trim_ambiguous','trim_reads','split_map',
	    'parse_gem_hit','get_split_coord_parser','combine_mapping_files',
	    'coords2splitcoords','get_split_mapper_routines');
push @EXPORT_OK,('line_reader','get_closest_hit');

use strict;
use warnings;
use Bio::Range;
use RNAseq_pipeline3 qw(get_fh print_gff run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file);


# This set of subroutines will parse the GEM output and return a hash containing
# the information for each hit
### OK
sub parse_gem_line {
    my $line=shift;
    chomp($line);
    my %line;

    my @line=split(/\t/,$line);
    if (@line == 5) {
	# We have qualities
	$line{'qual'}=splice(@line,2,1);
    } elsif (@line == 4) {
	# We have no qualities
    } else {
	die "Error in GEM format: $line\n";
    }
    my $type='F';
    if ($line[3]=~/=/o) {
	# Set to split-map if the read is split
	$type='S';
    }

    $line{'id'}=$line[0];
    # Remove spaces in the identifier and substitute them with underscores. This
    # Is done to avoid problems later with programs that do not admit the spaces
    $line{'id'}=~s/ +/_/g;

    # Fix the identifiers if they do not contain the standard /1 /2 ending
    $line{'id'}=~s/\|p1$/\/1/o;
    $line{'id'}=~s/\|p2$/\/2/o;

    # Fix Hi Seq  reads problem
    if ($line{'id'}=~s/[ _]1:/_X:/o) {
	$line{'id'}=~s/$/\/1/o;
    } elsif ($line{'id'}=~s/[ _]2:/_X:/o) {
	$line{'id'}=~s/$/\/2/o;
    }

    $line{'seq'}=$line[1];
    $line{'matches'}=$line[2];
    $line{'hits'}=$line[3];
    $line{'length'}=length($line{'seq'});
    $line{'type'}=$type;

    return(\%line);
}

sub get_best_hits {
    my $coords=shift;
    my $mismatches=shift;
    my @coords=@{$coords};
    my $type=0;

    my $hitstring=$coords->[0]->{'matches'};
    my @hitstring=split(':',$hitstring);

    for (my $i=0;$i<$mismatches;$i++) {
	my $hits=$hitstring[$i];
	if ($hits > 0) {
	    # This takes the hits with one more mismatch also
	    $hits+=$hitstring[$i + 1];
	    @coords=splice(@coords,0,$hits);
	    last;
	}
    }

    return(\@coords);
}

# This sub should take a gem standard (not split) hit and parse it
### TO DO we need a single sub that will extract all the hit info
sub process_gem_hit {
    my $match=shift;
    my $length=shift;
    my $matchsum=shift;

    if ($match =~/^-$/o) {
	return();
    }

    # Get chromosome
    my ($chr,$start)=split(':',$match);
    my @chr=split(/\s+/,$chr);
    if (@chr > 1) {
	$chr=$chr[0];
    }

    # Get strand
    my $strand=substr($start,0,1,'');
    if ($strand eq 'F') {
	$strand='+';
    } elsif ($strand eq 'R') {
	$strand='-';
    } else {
	warn "Strange strand for $match\n";
    }

    # Get qualities
    my ($coord,$qualities);
    if ($start=~/@/) {
	($coord,$qualities)= split(/\@/,$start);
    } else {
	$coord=$start;
    }
    $coord=~/(\d+)([^\d].+)*/o;
    $start=$1;
    my $mismatch_string = $2;
    my $end=$start + $length - 1;
    my $mismatch_no=0;
    if ($mismatch_string) {
	my $string=$mismatch_string;
	# Deal with insertions
	$string=~s/<.\d+>/I/g;
	$mismatch_no=$string=~s/[^\d]\d+/-/og;
    }

    unless ($qualities) {
	$qualities='None';
    }
    return([$chr,$start,$end,$strand,$mismatch_no,$qualities,$matchsum]);
}

# TO DO this sub should substitute the other one for split mapping
sub get_split_mapper_routines {
    my %mappers;

    $mappers{'GEM'}= sub {
	my $index=shift;
	my $infile=shift;
	my $outfile=shift;
	my $tmpdir=shift;
	my $mismatches=shift;
	my $zip=shift;
	
	# For the moment do not split map with more than 2 mismatches, as the
	# time is too long
	if ($mismatches > 2) {
	    print STDERR "WARNING: Setting mismatches > 2 during split mapping may take very long\n";
	}
	
	my $tmpfile=$outfile;
	if ($tmpdir) {
	    $tmpfile=$outfile;
	    $tmpfile=~s/.*\///;
	    $tmpfile=$tmpdir.'/'.$tmpfile;
	}
	
	# Build and execute the mapping command line
	my $command ='gem-split-mapper ';
	$command.="-I $index ";
	$command.="-i $infile ";
	$command.="-o $tmpfile ";
	$command.="-t 250 ";
	$command.="--min-split-size 12 ";
	$command.='-s \"GT\"+\"AG\",\"CT\"+\"AC\",\"GC\"+\"AG\",\"CT\"+\"GC\",\"ATATC\"+\"A.\",\".T\"+\"GATAT\",\"GTATC\"+\"AT\",\"AT\"+\"GATAC\" ';
	if ($mismatches) {
	    $command.="-m $mismatches ";
	}
	$command.=" > $tmpfile.log";
	print STDERR "Executing $command\n";
	system($command);

	# Wait for some seconds in case the files are not complete
	sleep(10);

	# Set LC_ALL to make sure it is "a la C"
	$ENV{'LC_ALL'}='C';

	# Collect the output into one file
	$command ="cat $tmpfile.?.split-map ";
	if ($zip) {
	    $command.="| sort -k1,1 -T $tmpdir |gzip -9 -c ";
	    $command.="> $outfile.split-map.gz";
	} else {
	    $command.="> $outfile.split-map";
	}

	sleep(10);
	
	print STDERR "Executing $command\n";
	system($command);
	$command ="rm $tmpfile.?.split-map";
	system($command);
    };
    
    return(\%mappers);
}

sub get_mapper_routines {
    my %mappers;

    $mappers{'GEM'}= sub {
	my $index=shift;
	my $infile=shift;
	my $outfile=shift;
	my $threads=shift;
	my $qualities=shift;
	my $tmpdir=shift;
	my $mismatches=shift;
	my $zip=shift;

	my $tmpfile=$outfile;
	if ($tmpdir) {
	    $tmpfile=$outfile;
	    $tmpfile=~s/.*\///;
	    $tmpfile=$tmpdir.'/'.$tmpfile;
	}

	# Build and execute the mapping command line
	my $command ='gem-mapper ';
	$command.="-I $index ";
	$command.="-i $infile ";
	$command.="-o $tmpfile ";
	if ($mismatches) {
	    $command.="-m $mismatches ";
	}
	$command.="-t $threads";
	if ($qualities) {
	    $command.=" -q $qualities";
	}
	$command.=" > $tmpfile.log";
	print STDERR "Executing: $command\n";
	system($command);

	# Set LC_ALL to make sure it is "a la C"
	$ENV{'LC_ALL'}='C';

	sleep(10);
	# Collect the output into one zipped sorted file
	# maybe use better -k2,2 -k1,1
	$command ="cat $tmpfile.?.map | sort -k1,1 -T $tmpdir ";
	if ($zip) {
	    $command.='|gzip -7 -c ';
	    $command.="> $outfile.map.gz";
	} else {
	    $command.="> $outfile.map";
	}

	sleep(10);

	print STDERR "Executing $command\n";
	system($command);
	$command ="rm $tmpfile.?.map";
	system($command);
    };

    return(\%mappers);
}

### TO DO
# The best thresholds should be selected for the split mapping
sub split_map {
    my $index=shift;
    my $infile=shift;
    my $outfile=shift;
    my $threads=shift;
    my $qualities=shift;
    my $tmpdir=shift;
    my $mismatches=shift;
    my $zip=shift;

    # For the moment do not split map with more than 2 mismatches, as the time
    # is too long
    if ($mismatches > 2) {
	print STDERR "WARNING: Setting mismatches > 2 during split mapping may take very long\n";
    }

    my $tmpfile=$outfile;
    if ($tmpdir) {
	$tmpfile=$outfile;
	$tmpfile=~s/.*\///;
	$tmpfile=$tmpdir.'/'.$tmpfile;
    }

    # Build and execute the mapping command line
    my $command ='gem-split-mapper ';
    $command.="-I $index ";
    $command.="-i $infile ";
    $command.="-o $tmpfile ";
    $command.="-t 250 ";
    $command.="--min-split-size 12 ";
    $command.='-s \"GT\"+\"AG\",\"CT\"+\"AC\",\"GC\"+\"AG\",\"CT\"+\"GC\",\"ATATC\"+\"A.\",\".T\"+\"GATAT\",\"GTATC\"+\"AT\",\"AT\"+\"GATAC\" ';
    if ($mismatches) {
	$command.="-m $mismatches ";
    }
    print STDERR "Executing $command\n";
    system($command);

    sleep(10);

    # Set LC_ALL to make sure it is "a la C"
    $ENV{'LC_ALL'}='C';

    # Collect the output into one file
    $command ="cat $tmpfile.?.split-map |sort -k1,1 -T $tmpdir ";
    if ($zip) {
	$command.='|gzip -7 -c ';
	$command.="> $outfile.split-map";
    } else {
	$command.="> $outfile.split-map";
    }
    
    sleep(10);

    print STDERR "Executing $command\n";
    system($command);
    $command ="rm $tmpfile.?.split-map";
    system($command);
}

sub get_unmapped {
    my $mapped=shift;
    my $unmapped=shift;
    my $log_fh=shift;

    my $left=0;
    my $total=0;

    print $log_fh "Getting unmapped reads from $mapped\n";

    my $mappedfh=get_fh("$mapped");
    my $unmappedfh=get_fh($unmapped,1);
    while (my $line=<$mappedfh>) {
	my %line=%{parse_gem_line($line)};
	$total++;
	if ($line{'matches'}=~/^0(:0)*$/o) {
	    if ($line{'qual'}) {
		print $unmappedfh '@',$line{'id'},"\n";
		print $unmappedfh $line{'seq'},"\n";
		print $unmappedfh '+',"\n";
		print $unmappedfh $line{'qual'},"\n";
	    } else {
		print $unmappedfh '>',$line{'id'},"\n";
		print $unmappedfh $line{'seq'},"\n";
	    }
	    $left++;
	}
    }
    close($mappedfh);
    close($unmappedfh);
    print $log_fh $total,"\tReads unmapped\n";
    print $log_fh $left,"\tReads still unmapped printed to $unmapped\n";

    return($left);
}

sub trim_reads {
    my $mapped=shift;
    my $unmapped=shift;
    my $trim=shift;
    my $threshold=40;

    unless($trim) {
	$trim=15;
    }

    my $left=0;

    print STDERR "Trimming 3' of $mapped reads by $trim nt\n";

    my $mappedfh=get_fh("$mapped");
    my $unmappedfh=get_fh($unmapped,1);
    while (my $line=<$mappedfh>) {
	chomp($line);
	if ($line=~/^@/o) {
	    # Entry is fastq
	    my $id=$line;
	    my $seq=<$mappedfh>;
	    chomp($seq);
	    <$mappedfh>; # get rid of the id for the fastq line
	    my $qual=<$mappedfh>;
	    chomp($qual);

	    # Trim the 3' end
	    my $length=length($seq) - $trim;
	    $seq=substr($seq,0,$length);
	    # Remove trailing Ns
	    $seq=~s/[N.]+$//o;
	    $length=length($seq);
	    $qual=substr($qual,0,$length);

	    if ($length > $threshold) {
		print $unmappedfh $id,"\n";
		print $unmappedfh $seq,"\n";
		print $unmappedfh '+',"\n";
		print $unmappedfh $qual,"\n";
		$left++;
	    } 
	} elsif ($line=~/^>/o) {
	    # Entry is fasta
	    my $id=$line;
	    my $seq=<$mappedfh>;

	    # Trim the 3' end
	    my $length=length($seq) - $trim;
	    $seq=substr($seq,0,$length);
	    # Remove trailing Ns
	    $seq=~s/[N.]+$//o;
	    $length=length($seq);

	    if ($length > $threshold) {
		print $unmappedfh $id,"\n";
		print $unmappedfh $seq,"\n";
		$left++;
	    }

	} else {
	    die "Unknown sequence format\n";
	}
    }
    close($mappedfh);
    close($unmappedfh);
    print STDERR $left,"\tReads unmapped\n";

    return($left);
}

sub trim_ambiguous {
    my $infile=shift;
    my $trimmed=shift;
    my $log_fh=shift;
    my $threshold=40;

    my $left=0;
    my $total=0;

    print $log_fh "Trimming ambiguous ends from unmapped reads in $infile\n";
    my $infh=get_fh("$infile");
    my $trimmedfh=get_fh($trimmed,1);
    while (my $line=<$infh>) {
	my %line=%{parse_gem_line($line)};

	if ($line{'matches'}=~/^0(:0)*$/o) {
	    $total++;

	    # Trim the 3' end
	    $line{'seq'}=~s/[N.]+$//o;
	    if ($line{'qual'}) {
		my $length3=length($line{'seq'});
		$line{'qual'}=substr($line{'qual'},0,$length3);
	    }

	    # Trim the 5' end
	    $line{'seq'}=~s/^[N.]+//o;
	    if ($line{'qual'}) {
		my $length5=length($line{'seq'});
		$line{'qual'}=substr($line{'qual'},-$length5);
	    }
	    
	    if (length($line{'seq'}) > $threshold) {
		if ($line{'qual'}) {
		    print $trimmedfh '@',$line{'id'},"\n";
		    print $trimmedfh $line{'seq'},"\n";
		    print $trimmedfh '+',"\n";
		    print $trimmedfh $line{'qual'},"\n";
		} else {
		    print $trimmedfh '>',$line{'id'},"\n";
		    print $trimmedfh $line{'seq'},"\n";
		}
		$left++;
	    }
	}
    }
    close($infh);
    close($trimmedfh);
    print $log_fh $total,"\tReads unmapped\n";
    print $log_fh $left,"\tReads in $trimmed after trimming\n";

    return($left);
}

### OK
### TO DO
# Allow for the quality to be specified on the command line.
sub determine_quality_type {
    my $infile=shift;

    print STDERR "This step is obsolete\n";
    # The Sanger format encodes a Phred quality score from 0 to 93 using ASCII
    # 33 to 126. Solexa/Illumina 1.0 format encodes a Solexa/Illumina quality
    # score from -5 to 40 using ASCII 59 to 104. Illumina 1.3 format encodes a
    # Phred quality score from 0 to 40 using ASCII 64 to 104.
    # This means that if we find scores lower than 59 the encoding is phred, and
    # if the scores are higher than 93 thay are likely to be solexa because it
    # is highly unlikely to have an unassembled base directly from a read with
    # a quality score greater than 60 (solexa reports only up to 40 and they did
    # the sequencer, so are probablly estimating upward)

    my $infh=get_fh($infile);
    my $line_type=0;
    my $quality_type='solexa';
    my $certainty=0;
    my $total=0;

    while (my $line=<$infh>) {
	chomp($line);
	if ($line=~/^\@/o) {
	    $line_type=1;
	    next;
	} elsif ($line_type == 1) {
	    if ($line=~/^\+/o) {
		$line_type = 2;
		next;
	    } else {
		next;
	    }
	} elsif ($line_type == 2) {
	    my @scores=split('',$line);
	    foreach my $score (@scores) {
		$total++;
		if (ord($score) < 59) {
		    $quality_type='phred';
		    return($quality_type);
		} elsif (ord($score) > 104) {
		    $quality_type='phred';
		    return($quality_type);
		} elsif (ord($score) > 95) {
		    # We will se this to 95 to be more strict with the higher 
		    # qualities from newer machines
		    print STDERR "Probably solexa\n";
		    $quality_type='solexa';
		    return($quality_type);
		}
	    }
	} else {
	    die "This does not look like fastq...\n";
	}
    }
    close($infh);

    warn "Unable to guess quality type. Assuming $quality_type\n";

   return($quality_type);
}

sub check_split_input {
    my $input=shift;
    my $filetype=shift;
    my $valid_input;
    if (($input &&
	 -r $input) &&
	($input=~/.split-map(.gz)?$/)) {
	$valid_input=1;
	$$filetype='gem.split';
    } else {
	warn "Problem with input\n";
	if ($input) {
	    warn "$input is not valid\n";
	} else {
	    die "No infile supplied\n";
	}
    }
    return($valid_input);
}

### OK
sub check_input {
    my $input=shift;
    my $filetype=shift;
    my $valid_input;
    if (($input &&
	 -r $input) &&
	($input=~/.(fa(stq)*)$/)) {
	$valid_input=1;
	$$filetype=$1;
    } else {
	warn "Problem with input\n";
	if ($input) {
	    warn "$input is not valid\n";
	} else {
	    die "No infile supplied\n";
	}
    }
    return($valid_input);
}

### OK
### TO DO
# Determine a more complete check to do in order to find not only if the files
# are there, but also if they are complete
sub check_index {
    my $index=shift;
    my $valid_index;

    print STDERR "Checking index $index...";
    if (-r "$index.blf") {
	$valid_index=1;
	print STDERR "Present\n";
    } else {
	print STDERR "Missing\n";
    }
    return($valid_index);
}

# This should transform a gem entry into a coordinates entry
sub gem_2_coords {
    my $entry=shift;
    my $mismatches=shift;
    my @coords;

    if ($mismatches) {
	print STDERR "Im using a useless argument gem_2_coords\n";
    }

    if ($entry->{'hits'} ne '-') {
	my $length=$entry->{'length'};
	my @matches=split(',',$entry->{'hits'});
	foreach my $match (@matches) {
	    my $hit=parse_gem_hit($match,
				  $length,
				  $entry->{'matches'});
	    $hit->{'id'}=$entry->{'id'};
	    $hit->{'matches'}=$entry->{'matches'};
	    push @coords,$hit;
	}
    }
    return(\@coords);
}


# This sub should take a coords hash and return a gtf string
sub coords2gtf {
    my $coords=shift;
    my $source=shift;
    my $info=join('; ',
		  'ID "'.$coords->{'id'}.'"',
		  'Identity "'.$coords->{'matches'}.'"');
    my $gtf=join("\t",
		 $coords->{'chr'},
		 $source,
		 'map',
		 $coords->{'start'},
		 $coords->{'end'},
		 '.',
		 $coords->{'strand'},
		 '.',
		 $info);
    return($gtf);
}

# This sub should take a coords hash and return a bed string
sub coords2bedSimple {
    my $coords=shift;
    my $source=shift;
    my $size=$coords->{'end'} - $coords->{'start'} + 1;
    my $start=$coords->{'start'} - 1;
    my $end=$coords->{'end'};
    my $gtf=join("\t",
		 $coords->{'chr'},
		 $start,
		 $end,
		 $coords->{'id'},
		 0,
		 $coords->{'strand'},
		 $start,
		 $end,
		 '0,255,0',
		 1,
		 $size,
		 0);
    return($gtf);
}

# This sub should take a coords hash with junctions information and return a
# bed string
sub coords2bedJunctions {
    my $coords=shift;
    my $source=shift;
    my $start=$coords->{'start'} - 1;
    my $end=$coords->{'end'};
    my $size=$end - $start + 1;
    
    my $bed=get_split_bed_coords($coords);
    if ($bed) {
	return($bed);
    } else {
	return();
    }
}
# This subroutine should transforma junction mapping to a juntion that is named
# using chr_start_end_strand_splice_chr_start_end_strand into a mapping that
# will be either a split mapping or a full mapping, and it will always keep the
# full mapping over the split if more than one
sub coords2splitcoords {
    my $coords=shift;
    my $start=$coords->{'start'} - 1;
    my $end=$coords->{'end'};
    my $size=$end - $start + 1;
    
    my ($hits,$mismatchno)=get_junction_coords($coords);
    if ($hits) {
	return($hits,$mismatchno);
    }
}

sub get_junction_coords {
    my $coords=shift;

    # To do rewrite in a more compact way
    my ($exon1,$exon2)=split('_splice_',$coords->{'chr'});
    my @coords1=split('_',$exon1);
    my ($chr1,$start1,$end1,$strand1);
    my @coords2=split('_',$exon2);
    my ($chr2,$start2,$end2,$strand2);

    $strand1=pop(@coords1);
    $end1=pop(@coords1);
    $start1=pop(@coords1);
    $chr1=join('_',@coords1);

    $strand2=pop(@coords2);
    $end2=pop(@coords2);
    $start2=pop(@coords2);
    $chr2=join('_',@coords2);

    unless($strand1 eq $strand2) {
	die "Conflicting strands: ".$coords->{'chr'}."\n";
    }

    unless($chr1 eq $chr2) {
	die "Conflicting Chromosomes: ".$coords->{'chr'}."\n";
    }

    # Get the length of each junction half. This should be the end minus the
    # start plus 1, and also the total length
    my $half_length1=$end1 - $start1 + 1;
    my $half_length2=$end2 - $start2 + 1;
    my $total_length=$half_length1 + $half_length2;
    my $inslength=$coords->{'inslength'};

    # Using the length of each fragment calculate where the map starts and
    # ends
    # Adjust start1
    my ($rstart1,$rend1,$rstart2,$rend2,$roffset)=(0,0,0,0,0);
    if ($strand1 == 1) {
	if ($coords->{'end'} + $inslength <= $half_length1) {
	    # The read is completely included in the first fragment
	    $rstart1=$start1 + $coords->{'start'} - 1;
	    $rend1=$start1 + $coords->{'end'} - 1 + $inslength;
	} elsif ($coords->{'start'} > $half_length1) {
	    # The read is completely included in the second fragment
	    $rstart1=$start2 + $coords->{'start'} - $half_length1 - 1;
	    $rend1=$start2 + $coords->{'end'} - $half_length1 - 1 + $inslength;
	} elsif (($coords->{'start'} <= $half_length1) &&
		 ($coords->{'end'} + $inslength > $half_length1)) {
	    # The read spans the junction
	    $rstart1=$start1 + $coords->{'start'} - 1;
	    $rend1=$end1;
	    $rstart2=$start2;
	    $rend2=$start2 + $coords->{'end'} - $half_length1 - 1 + $inslength;
	    $roffset=$rend1 - $rstart1 + 1
	} else {
	    die "Read problem\n";
	}
    } elsif ($strand1 == -1) {
	# Here we basically want to mirror the coordinates, so we need the total
	# length of the junction and we will substract each coordinate from this
	# That way we can actually position the reads the same way as we would
	# in the plus strand
	my $alt_start=$total_length - $coords->{'end'} - $inslength + 1;
	my $alt_end=$total_length - $coords->{'start'} + 1;
	if ($alt_end <= $half_length1) {
	    # The read is completely included in the first fragment (which
	    # is the lower coordinate one)
	    $rstart1=$start1 + $alt_start - 1;
	    $rend1=$start1 + $alt_end - 1;
	} elsif ($alt_start > $half_length1) {
	    # The read is completely included in the second fragment
	    $rstart1=$start2 + $alt_start - $half_length1 - 1;
	    $rend1=$start2 + $alt_end - $half_length1 - 1;
	} elsif (($alt_start <= $half_length1) &&
		 ($alt_end > $half_length1)) {
	    # The read spans the junction
	    $rstart1=$start1 + $alt_start - 1;
	    $rend1=$end1;
	    $rstart2=$start2;
	    $rend2=$start2 + $alt_end - $half_length1 - 1;
	    $roffset=$rend1 - $rstart1 + 1
	} else {
	    die "Read problem\n";
	}
    } else {
	warn "Unknown strand\n";
    }

    # Calculate the strand to output, as this should be with regard to the
    # genome and not to the junction
    my $outstrand;
    if ($strand1 == 1) {
	if ($coords->{'strand'} eq '+') {
	    $outstrand='F';
	} elsif ($coords->{'strand'} eq '-') {
	    $outstrand='R';
	}
    } elsif ($strand1 == -1) {
	if ($coords->{'strand'} eq '+') {
	    $outstrand='R';
	} elsif ($coords->{'strand'} eq '-') {
	    $outstrand='F';
	}
    } else {
	die "Strand problem\n";
    }

    # Need this for those cases that are split in order to know where they
    # belong in the list
    # Rebuild the mismatch list
    my $mismatchmap='';
    my $mismatchno=0;
    if (@{$coords->{'mismatch_list'}}) {
	$mismatchno=@{$coords->{'mismatch_list'}};
    }
    foreach my $mismatch (@{$coords->{'mismatch_list'}}) {
	$mismatchmap.=$mismatch->[0];
	$mismatchmap.=$mismatch->[1];
    }

    my $hitstring='';
    if ($rstart1 && $rstart2) {
	# Read spans the junctions
	$hitstring ="[$roffset]=$chr1:$outstrand".$rstart1.'~';
	$hitstring.="$chr2:$outstrand".$rstart2;
    } elsif ($rstart1) {
	# Build the quality string
	my $qualstring='';
	$qualstring.='@'.$coords->{'qual'}.'/';
	$qualstring.=$coords->{'mismatches'};
	
	# Read included in the first half
	$hitstring=$chr1.':'.$outstrand.$rstart1.$mismatchmap.$qualstring;

	# Skip those reads that are not spanning a junction
	return()
    } else {
	die "Problem with read:\t\n";
    }

    my ($size1,$size2)=(0,0);
    if ($rstart1) {
	$size1=$rend1 - $rstart1 + 1;
    }
    if ($rstart2) {
	$size2=$rend2 - $rstart2 + 1;
    }

    my $size=join(',',$size1,$size2);
    my $start=join(',',0,$start2 - $start1);

    if ($size1 + $size2 != $coords->{'length'} + $inslength) {
	warn "coordinate_problem: added half lengths longer than total\n";
	print STDERR join("\t",
			  $coords->{'length'},
			  $size1,
			  $size2,
			  $chr1,
			  $start1 - 1,
			  $end2,
			  $coords->{'id'},
			  $strand1),"\n";
	print STDERR "Please press return to continue...";
	<STDIN>;
    }
    return($hitstring,$mismatchno);
}

# This subroutione is essentially the same as the previous, but it returns a bed
# file instead of a split-map entry
sub get_split_bed_coords {
    my $coords=shift;
    my ($chr1,$start1,$end1,$strand1,$splice,
	$chr2,$start2,$end2,$strand2)=split('_',$coords->{'chr'});

    unless($strand1 eq $strand2) {
	die "Conflicting strands: ".$coords->{'chr'}."\n";
    }

    unless($chr1 eq $chr2) {
	die "Conflicting Chromosomes: ".$coords->{'chr'}."\n";
    }

    # Get the length of each junction half. This should be the end minus the
    # start plus 1, and also the total length
    my $half_length1=$end1 - $start1 + 1;
    my $half_length2=$end2 - $start2 + 1;
    my $total_length=$half_length1 + $half_length2;
    my $inslength=$coords->{'inslength'};

    # Using the length of each fragment calculate where the map starts and
    # ends
    # Adjust start1
    my ($rstart1,$rend1,$rstart2,$rend2,$roffset)=(0,0,0,0,0);
    if ($strand1 == 1) {
	if ($coords->{'end'} + $inslength <= $half_length1) {
	    # The read is completely included in the first fragment
	    $rstart1=$start1 + $coords->{'start'} - 1;
	    $rend1=$start1 + $coords->{'end'} - 1 + $inslength;
	} elsif ($coords->{'start'} > $half_length1) {
	    # The read is completely included in the second fragment
	    $rstart1=$start2 + $coords->{'start'} - $half_length1 - 1;
	    $rend1=$start2 + $coords->{'end'} - $half_length1 - 1 + $inslength;
	} elsif (($coords->{'start'} <= $half_length1) &&
		 ($coords->{'end'} + $inslength > $half_length1)) {
	    # The read spans the junction
	    $rstart1=$start1 + $coords->{'start'} - 1;
	    $rend1=$end1;
	    $rstart2=$start2;
	    $rend2=$start2 + $coords->{'end'} - $half_length1 - 1 + $inslength;
	    $roffset=$rend1 - $rstart1 + 1
	} else {
	    die "Read problem\n";
	}
    } elsif ($strand1 == -1) {
	# Here we basically want to mirror the coordinates, so we need the total
	# length of the junction and we will substract each coordinate from this
	# That way we can actually position the reads the same way as we would
	# in the plus strand
	my $alt_start=$total_length - $coords->{'end'} - $inslength + 1;
	my $alt_end=$total_length - $coords->{'start'} + 1;
	if ($alt_end <= $half_length1) {
	    # The read is completely included in the first fragment (which
	    # is the lower coordinate one)
	    $rstart1=$start1 + $alt_start - 1;
	    $rend1=$start1 + $alt_end - 1;
	} elsif ($alt_start > $half_length1) {
	    # The read is completely included in the second fragment
	    $rstart1=$start2 + $alt_start - $half_length1 - 1;
	    $rend1=$start2 + $alt_end - $half_length1 - 1;
	} elsif (($alt_start <= $half_length1) &&
		 ($alt_end > $half_length1)) {
	    # The read spans the junction
	    $rstart1=$start1 + $alt_start - 1;
	    $rend1=$end1;
	    $rstart2=$start2;
	    $rend2=$start2 + $alt_end - $half_length1 - 1;
	    $roffset=$rend1 - $rstart1 + 1
	} else {
	    die "Read problem\n";
	}
    } else {
	warn "Unknown strand\n";
    }

    unless ($rstart1 && $rstart2) {
	# Skip cases that do not span the junctions
	return()
    }

    # Calculate the strand to output, as this should be with regard to the
    # genome and not to the junction
    my $outstrand;
    if ($strand1 == 1) {
	if ($coords->{'strand'} eq '+') {
	    $outstrand='+';
	} elsif ($coords->{'strand'} eq '-') {
	    $outstrand='-';
	}
    } elsif ($strand1 == -1) {
	if ($coords->{'strand'} eq '+') {
	    $outstrand='-';
	} elsif ($coords->{'strand'} eq '-') {
	    $outstrand='+';
	}
    } else {
	die "Strand problem\n";
    }

    # Need this for those cases that are split in order to know where they
    # belong in the list
    # Rebuild the mismatch list
    my $mismatchmap='';
    my $mismatchno=0;
    if (@{$coords->{'mismatch_list'}}) {
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
    my $start=join(',',0,$start2 - $rstart1);

    if ($size1 + $size2 != $coords->{'length'} + $inslength) {
	warn "coordinate_problem: added half lengths longer than total\n";
	print STDERR join("\t",
			  $coords->{'length'},
			  $size1,
			  $size2,
			  $chr1,
			  $start1 - 1,
			  $end2,
			  $coords->{'id'},
			  $strand1),"\n";
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
		     $outstrand,
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
		 $chr1,$start1,$end1,$strand1,$size1,
		 $chr2,$start2,$end2,$strand2,$size2),"\n";
	return();
    }
}

sub parse_gem_hit {
    my $hit=shift;
    my $length=shift;
    my $mismatch_map=shift;
    my %hit;

    my ($chr,$rest1)=split(':',$hit);
    # Keep only the first word in the Id
    my @chr=split(/\s+/,$chr);
    if (@chr > 1) {
	$chr=$chr[0];
    }
    my $strand=substr($rest1,0,1,'');
    my ($rest,$qual_string,$qual,$mismatches);
    if ($rest1=~/@/) {
	# File contains qualities
	($rest,$qual_string)=split('@',$rest1);
	($qual,$mismatches)=split('/',$qual_string);
    } else {
	$rest=$rest1;
	$qual=$mismatches=0;
    }
    $rest=~s/(\d+)//;
    my $start=$1;
    
    my @mismatches;
    my $inslength=0;
    while ($rest=~/([ATCG]|<[+-]\d+>)(\d+)/g) {
	my $type=$1;
	my $pos=$2;
	if ($type =~/<[+-](\d+)>/) {
	    $inslength+=$1;
	}
	
	push @mismatches,[$type,$pos];
    }

    # Set the right strand:
    if ($strand eq 'F') {
	$strand='+';
    } elsif ($strand eq 'R') {
	$strand='-';
    } else {
	warn "Unknown strand: $strand\n";
	$strand='.';
    }

    $hit{'chr'}=$chr;
    $hit{'strand'}=$strand;
    $hit{'qual'}=$qual;
    $hit{'mismatches'}=$mismatches;
    $hit{'start'}=$start;
    $hit{'end'}=$start + $length - 1;
    $hit{'length'}=$length;
    $hit{'mismatch_list'}=[@mismatches];
    $hit{'mismatch_map'}=$mismatch_map;
    $hit{'inslength'}=$inslength;

    return(\%hit);
}

sub parse_split_gem_hit {
    my $string=shift;
    my $pos=shift;
    my $coord_parser=shift;
    my $read_length=shift;
    my %hit;

    # First we need to parse everything
    # Get the strand of the upstream and of the downstream
    my ($up_string,$down_string)=split(/~/,$string);
    my ($up_coords,$down_coords);
    ($hit{'up_chr'},$up_coords)=split(':',$up_string);
    ($hit{'down_chr'},$down_coords)=split(':',$down_string);
    $hit{'up_strand'}=substr($up_coords,0,1,'');
    $hit{'down_strand'}=substr($down_coords,0,1,'');

    # Depending of the way the strands are combined we must parse in different
    # ways
    my $parser_key=$hit{'up_strand'}.$hit{'down_strand'};
    if (exists $coord_parser->{$parser_key}){
	($hit{'up_start'},
	 $hit{'down_end'})=$coord_parser->{$parser_key}->($up_coords,
							  $down_coords,
							  $pos,
							  $read_length);
    } else {
	warn "Unknown strand combination $parser_key\n";
    }

    return(\%hit);
}

### Gem parsers
# This is the main body of the parsers
sub line_reader {
    my $fh=shift;
    my $mismatches=shift;
    my $closest_hit=shift;
    my $mismatch_strict=shift;
    my $unique_strict=shift;

    # For debugging
#    print STDERR join("\t",
#		      $mismatches,
#		      $closest_hit),"\n";

    # Set default values if they are not avaliable
    unless($mismatches) {
	print STDERR "No mismatch number provided: setting to 2\n";
	$mismatches=2;
    }

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
		splice(@map,$offset);
	    }
	    my $matches=0;
	    foreach my $hits (@map) {
		$matches+=$hits;
	    }

	    if ($matches < 1) {
		# If there are no matches skip the rest
		$read{'matches'}=[];
		return(\%read);
	    } elsif ($unique_strict && ($matches > 1)) {
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

sub get_closest_hit {
    my $read1=shift;
    my $read2=shift;
    my $unique_strict=shift;
    my $stranded=shift;

    my $closest;

    if ($unique_strict) {
	return()
    }

    my $ref_coord=int(($read1->[1] + $read1->[2]) / 2);
    my $dist= 1e10;
    foreach my $read (@{$read2}) {
	unless ($read1->[0] eq $read->[0]) {
	    # Exclude interchromosomal
	    next;
	}
	# Check that the strands are correct (both members should be mapping to
	# different strands
	if ($stranded) {
	    if ($read1->[3] ne $read->[3]) {
		# Both reads have the same strand, which should happed in
		# paired if they are stranded
		next;
	    }
	} else {
	    if ($read1->[3] eq $read->[3]) {
		# Both reads have the same strand, which should not happed in
		# paired unless there is something strange going on
		next;
	    }
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

# Get some routines to parse the splitmappings
sub gem_2_coords_split {
    my $entry=shift;
    my $coord_parser=shift;
    my $read_length=shift;

    my @coords;

    if ($entry->{'hits'} ne '-') {
	my $length=$entry->{'length'};
	my @matches=split(',',$entry->{'hits'});
	
	# Parse the hits in this case we are taking only the best one
	foreach my $match (@matches) {
	    my ($pos,$coord_string)=split('=',$match);
	    my $hit=parse_split_gem_hit($coord_string,
					$pos,
					$coord_parser,
					$read_length);
	    $hit->{'id'}=$entry->{'id'};
	    $hit->{'matches'}=$entry->{'matches'};
	    $hit->{'length'}=$entry->{'length'};
	    
	    # Set the strands correctly
	    my $up_strand=$hit->{'up_strand'};
	    if ($up_strand eq 'F') {
		$hit->{'up_strand'}='+';
	    } elsif ($up_strand eq 'R') {
		$hit->{'up_strand'}='-';
	    } else {
		warn "Unknown strand ",$hit->{'id'},"\n";
	    }
	    my $down_strand=$hit->{'down_strand'};
	    if ($down_strand eq 'F') {
		$hit->{'down_strand'}='+';
	    } elsif ($down_strand eq 'R') {
		$hit->{'down_strand'}='-';
	    } else {
		warn "Unknown strand ",$hit->{'id'},"\n";
	    }

	    # Skip cases with more than one combination of sites
	    ### TO DO This should be modified to make the sub more general
	    # so it can be used by all the split parsers
	    unless (@{$hit->{'up_start'}} == 1 &&
		    @{$hit->{'down_end'}} == 1) {
		next;
	    }
	    push @coords,$hit;	
	}
    }
    return(\@coords);
}

sub get_split_coord_parser {
    my %parser;

    $parser{'FF'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_single_start($up_coord_string,
				       $read_length,
				       $pos);

	# Check the down string
	my $ranges2=check_multi_start($down_coord_string,
				      $read_length,
				      $pos,
				      1);

	return($ranges1,$ranges2);
    };
    
    $parser{'FR'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_single_start($up_coord_string,
				       $read_length,
				       $pos);
	
	# Check the down string
	my $ranges2=check_single_start($down_coord_string,
				       $read_length,
				       $pos,
				       1);

	return($ranges1,$ranges2);
    };

    $parser{'RR'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_multi_start($up_coord_string,
				      $read_length,
				      $pos);

	# Check the down string
	my $ranges2=check_single_start($down_coord_string,
				       $read_length,
				       $pos,
				       1);

	return($ranges1,$ranges2);
    };
    
    $parser{'RF'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_multi_start($up_coord_string,
				      $read_length,
				      $pos);

	# Check the down string
	my $ranges2=check_multi_start($down_coord_string,
				      $read_length,
				      $pos,
				      1);

	return($ranges1,$ranges2);
    };

    return(\%parser);
}

sub check_multi_start {
    my $coord_string=shift;
    my $read_length=shift;
    my $pos=shift;
    my $down=shift;
    my @ranges;

    my %starts;
    $pos=~s/(\[|])//g;
    $coord_string=~s/(\[|])//g;
    
    my @pos=sort {$a <=> $b} (split(/[-;]/,$pos));
    my ($pos_start,$pos_end)=($pos[0],$pos[-1]);
    my @splits=($pos_start .. $pos_end);
    if ($down) {
	foreach my $split (@splits) {
	    $split=$read_length - $split;
	}
	@splits=sort {$b <=> $a} @splits;
    } else {
	@splits=sort {$a <=> $b} @splits;
    }

    my @coords=split(/[-;]/,$coord_string);
    my @coords_sort=sort {$a <=> $b} @coords;
    my ($hit_start,$hit_end)=($coords_sort[0],$coords_sort[-1]);
    my @starts=($hit_start .. $hit_end);

    if ($coords[1] &&
	($coords[0] > $coords[1])) {
	# This means we are in the minus strand, so we invert the array
	@starts=reverse(@starts);
    }

    unless (@splits == @starts) {
	warn "Wrong number of starts $pos, $coord_string\n";
	print STDERR join("\t",
			  @splits),"\n";
	print STDERR join("\t",
			  @starts),"\n";
	return();
    }

    for (my $i=0;$i<@splits;$i++) {
	push @ranges, [$starts[$i],$starts[$i] + $splits[$i] - 1];
    }

    return(\@ranges);
}

sub check_single_start {
    my $coord_string=shift;
    my $read_length=shift;
    my $pos=shift;
    my $down=shift;
    my @ranges;

    my %starts;
    $pos=~s/(\[|])//g;
    $coord_string=~s/(\[|])//g;

    my @pos=sort {$a <=> $b} (split(/[-;]/,$pos));
    my ($pos_start,$pos_end)=($pos[0],$pos[-1]);
    my @splits=($pos_start .. $pos_end);
    if ($down) {
	foreach my $split (@splits) {
	    $split= $read_length - $split;
	}
	@splits=sort {$b <=> $a} @splits;
    } else {
	@splits=sort {$a <=> $b} @splits;
    }
    my @starts=split(/[-;]/,$coord_string);
    @starts=($starts[0] .. $starts[-1]);

    unless ((@splits == @starts) ||
	    (@starts == 1)){
	warn "Wrong number of starts $pos, $coord_string\n";
	return();
    }

    for (my $i=0;$i<@splits;$i++) {
	my $start;
	if (@starts > 1) {
	    $start=$starts[$i];
	} else {
	    $start=$starts[0]
	}
	push @ranges, [$start,$start + $splits[$i] - 1];
    }

    return(\@ranges);
}

# Combine a set of mapping files
sub combine_mapping_files{
    my $files=shift;
    my $outfile=shift;
    my $tmpdir=shift;
    my $log_fh=shift;

    my $final_file="$$.recursive.map";

    print $log_fh "combining the mapped files...\n";
    $ENV{'LC_ALL'}='C';
    my $command='cat ';
    $command.=join(' ',@{$files});
    $command.=" |sort -T $tmpdir| uniq";
    $command.=" > $final_file";
    
    run_system_command($command,
		       $log_fh);

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
	    if ($oldread{'matches'}=~/^0(:0)*$/o) {
		if ($line{'matches'}!~/^0(:0)*$/o) {
		    $oldline=$line;
		    %oldread=%line;
		} elsif ($oldread{'length'} < $line{'length'}) {
		    $oldline=$line;
		    %oldread=%line;
		} else {
		    next;
		}
	    } elsif ($line{'matches'}!~/^0(:0)*$/o) {
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
    run_system_command($command,
		       $log_fh);

    if (-r $final_file) {
	$command = "rm $final_file";
	run_system_command($command,
			   $log_fh);
    }

    if (-r "$tmpdir/*.unmapped.gem.split-map.*.$$.log") {
	$command = "rm $tmpdir/*.unmapped.gem.split-map.*.$$.log";
	run_system_command($command,
			   $log_fh);
    }
}

1;
