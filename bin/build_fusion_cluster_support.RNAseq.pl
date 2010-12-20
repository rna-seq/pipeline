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
# This script should take as input the output of the build_fusion_clusters
# script and it will select from the fusion reads and fusion clusters
# those clusters where at least $threshold reads support the p
# If will project the reads and select from the projection those cases where
# there are at least $threshold reads supporting the projection.
# For each splitmap it will determine which of these ranges it overlaps if any
# and it will select those pairs of ranges supported by at least $threshold
# split reads

use Getopt::Long;
use Bio::Range;
use RNAseq_pipeline3 qw(get_fh parse_gff_line);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_coords_from_junc_id_sub);

my $threshold=1;
my $species;
my $project;
my $clusterdir;
my $tmpdir;
my $workdir;
my $filetype;

GetOptions('-threshold=i' => \$threshold,
	   '-filetype=s' => \$filetype
    );

unless ($filetype) {
    die "No filetype supplied\n";
}

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$clusterdir=$options{'CLUSTERDIR'};
$tmpdir=$options{'LOCALDIR'};

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};

if ($filetype eq 'trans') {
    $workdir=$options{'TRANSDIR'};
} elsif ($filetype eq 'genome') {
    $workdir=$options{'GENOMEDIR'};
} else {
    die "Unknown filetype $filetype supplied\n";
}

# Get some subs
# We will use this sub, as the projectiosn have the same format as the junctions
*get_exon_coords=get_coords_from_junc_id_sub();

# Get the file with the projected split generated in the previous step
my $projected_reads_fn=$workdir.'/'.$project.'.unique.'.$filetype.'.cluster.gtf.gz';
unless (-r $projected_reads_fn) {
    die "Cannot find the unique clustering of the split maps at $projected_reads_fn\n";
}

# Select the splitmap projections supported by $threshold  or more read
my %projections=%{get_projections($projected_reads_fn,
				  $threshold)};

my @outfiles;
foreach my $lane (keys %lanes) {
    my $type;
    if ($lanes{$lane} == 1) {
	$type='single';
    } elsif ($lanes{$lane} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my $infilename=$workdir.'/'.$lane.'.'.$type.'.fusion.supporting.gtf.gz';
    print STDERR "Processing $infilename\n";

    unless (-r $infilename) {
	die "$infilename is not readable\n";
    }

    my %splits=%{get_split_maps($infilename)};

    # Find overlaps between the reads and the projections
    my %overlaps=%{get_overlaps(\%projections,
				$infilename)};

    # Find pairs of projections using the overlaps
    # Here we select those cases where a hit maps to two positions each of them
    # with at least $threshold hits
    my %splices;
    my %candidates=%{get_candidates(\%overlaps,
				    $threshold,
				    \%splices)};

    # Print results as a table with the different candidates
    my $outfile=$workdir.'/'.$lane.'.'.$type.'.fusion.clusters.gtf.gz';
    print_split_table($lane,
		      $outfile,
		      \%candidates,
		      \%splices,
		      $filetype,
	);

    print STDERR "done\n";
    push @outfiles, $outfile;
}

# Print resulting table to the STDOUT
foreach my $file (@outfiles) {
    my $fh=get_fh($file);
    while (my $line=<$fh>) {
	print $line;
    }
    close($fh);
}

exit;

sub print_split_table {
    my $lane=shift;
    my $outfn=shift;
    my $candidates=shift;
    my $splices=shift;
    my $filetype=shift;

    print STDERR "Printing table in $outfn\n";
    my $outfh=get_fh($outfn,1);

    foreach my $overlap (keys %{$candidates}) {
	foreach my $hit (keys %{$candidates->{$overlap}}) {
	    my $type=$overlap;
	    if ($overlap eq 'far') {
		my @junction_string=split('_',$hit);
		if ($junction_string[0] ne $junction_string[4]) {
		    $type='interchrom';
		}
	    }
	    my @coords=
	    print $outfh join("\t",
			      $type,
			      $hit,
			      scalar @{$candidates->{$overlap}{$hit}},
			      @{$splices->{$overlap}{$hit}},
			      $filetype,
			      $lane),"\n";
	}
    }
    close($outfh);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}++;
    }

    return(\%lanes);
}

sub get_candidates {
    my $overlaps=shift;
    my $threshold=shift;
    my $splices=shift;

    my %candidates;
    my %remove;

    print STDERR "Choosing candidates...";

    foreach my $hit (keys %{$overlaps}) {
	my @overlaps=sort keys %{$overlaps->{$hit}};

	if (@overlaps == 2) {
	    # If there are two keys the position is good
	    my $candidate_id=join('_splice_',
				  @overlaps);

	    my $upstream=get_exon_coords($overlaps[0]);
	    my $downstream=get_exon_coords($overlaps[1]);
	    $splices->{'far'}{$candidate_id}=[$upstream->[0],
				     $upstream->[2],
				     $downstream->[0],
				     $downstream->[1]];
	    push @{$candidates{'far'}{$candidate_id}}, $hit;
	} elsif (@overlaps == 1) {
	    # There is one key but both fragments of the read hit it
	    if ($overlaps->{$hit}->{$overlaps[0]} == 2) {
		$splices->{'close'}{$overlaps[0]}=['\N','\N','\N','\N'];
		push @{$candidates{'close'}{$overlaps[0]}}, $hit;
	    }
	} else {
	    ### TO DO check what these are
	    warn "Strange number of overlaps for $hit\n";
	    print STDERR join("\t",@overlaps),"\n";
	}
    }

    # remove those cases that do not have enough support
    foreach my $candidate (keys %{$candidates{'far'}}) {
	my $number=@{$candidates{'far'}{$candidate}};
	if ($number < $threshold) {
	    $remove{$candidate}= 1;
	}
    }

    foreach my $candidate (keys %{$candidates{'close'}}) {
	my $number=@{$candidates{'close'}{$candidate}};
	if ($number < $threshold) {
	    $remove{$candidate}=1;
	}
    }

    foreach my $distance ('far','close') {
	foreach my $candidate (keys %remove) {
	    delete $candidates{$distance}{$candidate};
	}
    }

    print STDERR "done\n";
    my $count=keys %remove;
    print STDERR $count,"\tCases removed due to low support (less than $threshold reads)\n";
    my $count_far=keys %{$candidates{'far'}};
    my $count_close=keys %{$candidates{'close'}};
    print STDERR $count_far,"\tClassified as far\n";
    print STDERR $count_close,"\tClassified as close\n";

    return(\%candidates);
}

sub get_overlaps {
    my $projections=shift;
    my $infile=shift;
    my %overlaps;

    # First print the set of projections with threshold or more overlaps
    my $project_file2=$tmpdir."/projected.2.$$.gff";
    my $outfh=get_fh($project_file2,1);

    foreach my $interval (keys %{$projections}) {
	print $outfh join("\t",
			  $projections->{$interval}->[0],
			  'Project',
			  'condensed',
			  $projections->{$interval}->[1],
			  $projections->{$interval}->[2],
			  '.',
			  '.',
			  '.'),"\n";
    }
    close($outfh);

    # Unzip the infile if it is zipped
    my $readfile2=$tmpdir."/reads.2.$$.gff";
    if ($infile=~/.gz/) {
	my $command="gunzip $infile -c > $readfile2";
	system($command);
    } else {
	my $command="ln -s $infile $readfile2";
	system($command);
    }

    my $overlap_file=$tmpdir."/overlap.$$.gff";
    my $command="overlap $project_file2 $readfile2 -ucsc -m 10 -o $overlap_file";
    print STDERR "Executing: $command\n";
    system($command);
    $command="rm $project_file2 $readfile2";
    print STDERR "Executing: $command\n";
    system($command);

    # Build the hash with the overlaps
    my $overlapfh=get_fh($overlap_file);
    my $total_overlaps=0;
    my $no_overlaps=0;

    while (my $line=<$overlapfh>) {
	# We cant use the standard gff parser as the output of overlap is not
	# standard
	my %line=%{parse_gff_line($line)};
	my $overlap_id=join('_',
			    $line{'chr'},
			    $line{'start'},
			    $line{'end'});

	# Get the hits
	my $hits=$line{'feature'}{'list_feat2:'};

	# Check there are hits, this is because we are using the global clusters
	# so not necessarity all will have read support in every file
	if ($hits eq '.') {
	    $no_overlaps++;
	    next;
	} else {
	    $total_overlaps++;
	}

	$hits=~s/\|?p?[12]//g;
	my @hits=split(',',$hits);

	# Load the hits in the hash
	foreach my $hit (@hits){
	    $overlaps{$hit}{$overlap_id}++;
	}
    }
    print STDERR $total_overlaps,"\tOverlaps found\n";
    print STDERR $no_overlaps,"\tClusters not overlapped found\n";

    # Check that no hit has more than 2 projections
    foreach my $hit (keys %overlaps) {
	my $projections=keys %{$overlaps{$hit}};
	if ($projections > 2){
	    warn "Too many hits for $hit. Sure you're using unique maps?\n";
	}
    }

    close($overlapfh);
    $command="rm $overlap_file";
    system($command);

    return(\%overlaps);
}

sub get_split_maps {
    my $file=shift;
    my %splits;

    my $infh=get_fh($file);
    while (my $line=<$infh>) {
	my %line=%{parse_gff_line($line)};
	my $split_id=$line{'feature'}{'read_id'};
	$split_id=~s/\|?p?[12]$//;
	unless ($split_id) {
	    warn "No id for $line\n";
	}
	my ($chr,$start,$end,$strand)=($line{'chr'},
				       $line{'start'},
				       $line{'end'});
	my $length=$end - $start + 1;
	push @{$splits{$split_id}},[$chr,$start,$end,$strand,$length];
    }
    close($infh);

    # Check the input and remove those cases where the read is not present as a
    # pair
    my %remove;
    foreach my $split_id (keys %splits) {
	unless (@{$splits{$split_id}} == 2) {
	    warn "Problem with $split_id\n";
	    $remove{$split_id}=1;
	}
    }

    foreach my $read (keys %remove) {
	delete $splits{$read};
    }

    my $count=keys %remove;
    print STDERR $count,"\tUnpaired entries removed\n";

    return(\%splits);
}

# Select the splitmap projections supported by $threshold  or more reads
sub get_projections {
    my $infile=shift;
    my $threshold=shift;
    my %interval;
    my $total=0;

    print STDERR "Extracting intervals supported by $threshold or more reads...\n";
    my $proj_fh=get_fh($infile);
    while (my $line=<$proj_fh>) {
	$total++;
	my %line=%{parse_gff_line($line)};
	my $coverage=$line{'feature'}{'included'};
	$coverage=~s/"|'//g;
	unless ($coverage) {
	    warn "No coverage in $line\n";
	}
	
	if ($coverage >= $threshold) {
	    my ($chr,$start,$end)=($line{'chr'},
				   $line{'start'},
				   $line{'end'});
	    my $interval_id=join('_',
				 $chr,$start,$end);
	    $interval{$interval_id}=[$chr,$start,$end];
	}
    }
    print STDERR "done\n";

    my $count=keys %interval;
    print STDERR $count,"\tclusters selected from $total\n";

    return(\%interval);
}
