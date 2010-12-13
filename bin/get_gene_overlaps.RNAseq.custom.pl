#!/soft/bin/perl

use strict;
use warnings;

# Objective
# This script should run overlap on a set of read clusters, coming from the
# genome mapping, the junctions mapping and the split mapping filtered for
# selecting those junctions between clusters supported by a certain threshold
# number of reads and in which the junction is also supported by that number
# of reads.

# Load some modules
use RNAseq_pipeline2 ('get_fh','parse_gff_line');
use RNAseq_pipeline_settings ('read_config_file','read_file_list');
use Getopt::Long;

# Declare some variables
my $threshold=1;
my $tmpdir;
my $clusterdir;
my $genomedir;
my $exondir;
my $splitdir;
my $junctionsdir;
my $projectid;
my $prefix;
my $stranded;
my $file_list;
my $parallel='default';

GetOptions('threshold|t=i' => \$threshold,
    );

# Read the options file
my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$genomedir=$options{'PROJECT'}.'/'.$options{'GENOMEDIR'};
$exondir=$options{'EXONDIR'};
$splitdir=$options{'SPLITMAPDIR'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$projectid=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$stranded=$options{'STRANDED'};
$file_list=$options{'FILELIST'};
if ($options{'PARALLEL'}) {
    $parallel='parallel';
}

# Get the required sub
*get_feature_overlap=get_feature_overlap_sub($parallel);

# First ge the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get for each lane the files corresponding to the genome mapping and the
# junctions mapping
my %lane_files=%{get_lane_files(\%files,
				$genomedir)};

my $annotation=$genomedir.'/'.$prefix.'.gene.gtf';

my %coverage;

foreach my $pair (keys %lane_files) {
    # Find the coverage granted by genome mapping;
    my $outfile=$lane_files{$pair};
    $outfile=~s/.gz$//;
    get_feature_overlap($lane_files{$pair},
			$annotation,
			$stranded,
			$outfile,
			$tmpdir);
}

exit;

sub get_lane_files {
    my $files=shift;
    my $dir=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	$lane_files{$pair}{$lane}=1;
    }

    # Determine if we are looking at paired or single reads
    foreach my $pair (keys %lane_files) {
	my @lanes=keys %{$lane_files{$pair}};
	if (@lanes == 1) {
	    # single files
	    print STDERR "$pair set identified as single reads\n";
	    $lane_files{$pair}=$dir.'/'.$lanes[0].'.single.unique.gtf.gz';
	} elsif (@lanes == 2) {
	    # paired files
	    print STDERR "$pair set identified as paired reads\n";
	    $lane_files{$pair}=$dir.'/'.$pair.'.paired.unique.gtf.gz';
	} else {
	    # This should never happen
	    die "Apparently there are more than two files grouped\n";
	}
    }

    return(\%lane_files);
}

###
# Here are a series of subroutines used to run sarahs overlap program on files
# regardless of their size
sub get_feature_overlap_sub {
    my $parallel=shift;
    my $size=1000000;
    my %subs;

    my $flags='-v -m -10 -ucsc';
    
    $subs{'default'}= sub {
	my $features=shift;
	my $annotation=shift;
	my $stranded=shift;
	my $outfile=shift;
	my $tmpdir=shift;

	if ($stranded) {
#	print STDERR "WARNING:Overlapping as Unstranded due to -m -st problem in overlap\n"
	    $flags.=' -st 1';
	}

	# Get the otput file name
	my $overlapfn=$outfile;
	$overlapfn=~s/.gz//;
	$overlapfn.='.overlap.gz';
	
	# Check if the overlap file exists
	if (-r $overlapfn) {
	    print STDERR $overlapfn,"\tAlready exists. Skipping...\n";
	    return();
	}

	# run overlap
	print STDERR "Running overlap on $annotation and $features\n";
	
	# Split the file
	if (-r $features) {
	    print STDERR "Splitting $features into $size lines fragments\n";
	} else {
	    die "Can't read $features\n";
	}
	my @files=split_file($features,
			     $size,
			     $tmpdir);
	my @overlap_files;
	
	# run overlap for each of them
	foreach my $file (@files) {
	    my $outfn=$file.'.overlap';
	    run_overlap($annotation,
			$file,
			$flags,
			$outfn);
	    push @overlap_files, $outfn;
	}
	
	# Combine the overlap files
	my @over_files=combine_overlap(\@overlap_files,
				       $overlapfn);
	
	# clean up
	clean_up(@files,
		 @over_files);
    };


    $subs{'parallel'}=sub {
	my $features=shift;
	my $annotation=shift;
	my $stranded=shift;
	my $outfile=shift;
	my $tmpdir=shift;

	# Check which responding nodes nodes have a threshold lower than 4
	my @available_clients=get_client_loads();
	my $client_number=0;
	$client_number=@available_clients;

	if ($stranded) {
#	print STDERR "WARNING:Overlapping as Unstranded due to -m -st problem in overlap\n"
	    $flags.=' -st 1';
	}

	# Get the otput file name
	my $overlapfn=$outfile;
	$overlapfn=~s/.gz//;
	$overlapfn.='.overlap.gz';
	
	# Check if the overlap file exists
	if (-r $overlapfn) {
	    print STDERR $overlapfn,"\tAlready exists. Skipping...\n";
	    return();
	}
	
	# run overlap
	print STDERR "Running overlap on $annotation and $features\n";
	
	# Split the file
	if (-r $features) {
	    print STDERR "Splitting $features into $size lines fragments\n";
	} else {
	    die "Can't read $features\n";
	}
	my @files=split_file($features,
			     $size,
			     $tmpdir);
	my @overlap_files;
	
	# run overlap for each of them
	my @left=@files;
	my @pids;
	my %nodes;
	while(@left ||
	      (@available_clients < $client_number)) {
	    # While will not exit until all jobs are done and all clients are
	    # free
	    if (@available_clients && @left) {
		my $file=shift(@left);
		my $node=shift(@available_clients);
		my $child_pid=fork();
		$nodes{$child_pid}=$node;
		my $outfn=$file.'.overlap';
		# Temporary files seem named by overlap using the time, so 
		# if we don't sleep we can have problems with this
		sleep(1);
		if ($child_pid) {
		    push @pids, $child_pid;
		    push @overlap_files, $outfn;
		} else {  # child exec date		
		    print STDERR "Procesing $file in $node\n";
		    run_overlap($annotation,
				$file,
				$flags,
				$outfn,
				$node);
		}
	    } else {
		# Put in a while loop in case a child dies while another corpse
		# is being collected, to avoid leaving zombies
		my $ended=wait();
		sleep(1);
		if ($ended > 0) {
		    push @available_clients, $nodes{$ended};
		    print STDERR "Job\t",$ended,"\tFinished in $nodes{$ended}\n";
		}
	    }
	    sleep(1);
	}
	
	sleep(1);
	print STDERR "All partial jobs finished for $overlapfn\n";
	
	# Combine the overlap files
	my @over_files=combine_overlap(\@overlap_files,
				       $overlapfn);
	
	# clean up
	clean_up(@files,
		 @over_files);
    };

    if ($subs{$parallel}) {
	return($subs{$parallel});
    } else{
	print STDERR "Unknown type $parallel. Setting to default (serial)\n";
	return($subs{'default'});
    }
}

sub run_overlap {
    my $annotation=shift;
    my $features=shift;
    my $flags=shift;
    my $outfile=shift;
    my $host=shift;

    my $command;
    if ($host) {
	$command="ssh $host nice overlap $annotation $features $flags > $outfile";
	print STDERR "Executing: $command\n";
	exec($command);
    } else {
	$command="overlap $annotation $features $flags > $outfile";
	print STDERR "Executing: $command\n";
	system($command);
    }
}

sub clean_up {
    my @files=@_;
    print STDERR "Cleaning up...\n";
    my $initial=@files;
    my $result;
    print STDERR $initial,"\tFiles to be removed\n";
    $result=unlink @files;
    print STDERR 'unlink ',@files,"\n";
    print STDERR $result,"\tFiles removed\n";
    if ($initial == $result) {
	print STDERR "successfully removed temporary files\n";
    } else {
	warn "Problem removing files\n";
    }
}

sub split_file {
    my $file=shift;
    my $size=shift;
    my $tmpdir=shift;
    my @splitfiles;


    my $infh=get_fh($file);
    my $lineno=0;
    my $fileno=0;
    my $outfn=$file.'.split.'.$fileno;
    if ($tmpdir) {
	$outfn=~s/.*\///;
	$tmpdir=~s/\/$//;
	$outfn=$tmpdir.'/'.$outfn;
    }
    my $outfh=get_fh($outfn,1);
    push @splitfiles,$outfn;
    while (my $line=<$infh>) {
	$lineno++;
	unless ($lineno % $size) {
	    close($outfh);
	    $fileno++;
	    $outfn=$file.'.split.'.$fileno;
	    if ($tmpdir) {
		$outfn=~s/.*\///;
		$tmpdir=~s/\/$//;
		$outfn=$tmpdir.'/'.$outfn;
	    }
	    push @splitfiles,$outfn;
	    $outfh=get_fh($outfn,1);
	}
	print $outfh $line;
    }
	    
    @splitfiles=sort @splitfiles;

    return(@splitfiles);
}

sub combine_overlap {
    my $files=shift;
    my $overlapfn=shift;
    
    # Because we are running overlap on the same file there will always be
    # the same number of lines
    my @lines;
    my %feats;
    my @overlap_files;
    print STDERR "Combining overlap files...\n";
    foreach my $file (@{$files}) {
	my $overlapfile=$file;
	if (-r $overlapfile) {
	    print STDERR "Processing $overlapfile\n";
	} else {
	    warn "Unable to find readable file $overlapfile\n";
	}
	push @overlap_files,$overlapfile;

	my $overfh=get_fh($overlapfile);
	my $line_no=0;
	while (my $line=<$overfh>) {
	    chomp($line);
	    my %line=%{parse_gff_line($line)};

	    # Add the number of overlaps
	    my $hits;

	    # Add the locations if they are present
	    $hits=$line{'feature'}{'list_feat2:'};
	    if ($hits ne '.') {
		$feats{$line_no}[1].=$hits;
	    }
	    $feats{$line_no}[0]+=$line{'feature'}{'nb_ov_feat2:'};
	    unless ($lines[$line_no]) {
		$lines[$line_no]=$line;
	    }

	    $line_no++;
	}
	close($overfh);
    }

    print STDERR "Merging...\n";
    my $overlapfh=get_fh($overlapfn,1);
    for (my $i=0;$i<@lines;$i++) {
	my $line=$lines[$i];
	my $feats=$feats{$i}[0];
	my $hits=$feats{$i}[1];

	# Add a '.' for the cases with no hits in any file
	unless ($hits) {
	    $hits='.';
	}

	if ($hits) {
	    unless ($line=~s/(nb_ov_([^\s])+: )"\d+";( list_([^\s])+: )"([^\s])+";$/$1"$feats";$3"$hits";/) {
		warn "Problem with overlap line\n";
	    }
	} else {
	    warn "Problem no hits field in: $line\n";
	}

	print $overlapfh $line,"\n";
    }
    print STDERR "done\n";
    return (@overlap_files);
}

sub get_client_loads {
    my @out = ();
    my $TIMEOUT = 10;
    my $threshold = 7;
    my %thresholds=('corb'=> 6,
		    'icarus'=> 4,
		    'cel'=> 4,
		    'tro' => 1,
		    'cuc' => 7,
		    'vesc' => 1,
		    'tmatik' => 3
	);
    my @clients = ('cel',
		   'corb',
#		   'tro',
#		   'vesc',
#		   'icarus', 
#		   'cuc', 
		   'tmatik');
    # check load on all clients...
    print STDERR "Checking client load for the last 1/5/15 minute(s)\n";
    foreach my $client (@clients) {
	my $load = '';
	eval {
	    local $SIG{ALRM} = sub { die 'alarm time out' };
	    alarm $TIMEOUT;
	    $load = `ssh $client uptime`;
	    alarm 0;
	    1; # return value from eval on normalcy
	} or warn "load check on $client timed out after $TIMEOUT seconds.\n";
	unless ($load) {
	    next;
	}
	$load =~ /load average: (.*?), (.*?), (.*?)\s+/;
	my ($min_1,$min_5,$min_15) = ($1,$2,$3);
	my $max = 0;
	$max = $min_1 if ($min_1 > $max);
	$max = $min_5 if ($min_5 > $max);
	$max = $min_15 if ($min_15 > $max);
	while ($max < $thresholds{$client}) {
	    push @out, $client;
	    $thresholds{$client}--;
	}
    }

    # # Print out the available resources
    my $client_number=@out;
    print STDERR $client_number,"\tClients available\n";
    my %free;
    foreach my $client (@out){
	$free{$client}++;
    }
    foreach my $client (keys %free) {
	print STDERR join("\t",
			  $client,
			  $free{$client}),"\n";
    }

    return (@out)
}
