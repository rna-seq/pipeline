#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script will take a file with a set of transcript lengths, and a gtf file
# with a set of mappings, and it will generate a distribution of read mappings
# along the transcripts by dividing them into bins (fractions of their length
# and assigning each read to one of these bins (or more depending on the frac
# tion included

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use Bio::SeqIO;
use Getopt::Long;

my $transcripts;
my $reads;
my $length;
my $bin_no=10;
my $offset=2; # This gives the interval for which a read will be reported
my $species;
my $project;

GetOptions(
	   'bins|b=i' => \$bin_no,
	   'length=i' => \$length,
	   'offset|o=i' => \$offset,
	   'reads|r=s' => \$reads);

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};

# Get the lengths file and read it
$transcripts=$options{'TRANSDIR'}.'/transcript.lengths';
unless ($transcripts) {
    die "A transcript lengths file is required\n";
}
my %lengths=%{get_transcript_lengths($transcripts)};


# Set the outfile
my $outfile=$options{'PREFIX'}.'_unique_maps_transcripts';
my $disttable=$options{'PREFIX'}.'_read_dist_transcripts';
my $disttablefh=get_fh($disttable.'.txt',1);

my %distribution;
unless ($outfile) {
    die "No output file name\n";
}

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};

foreach my $lane (keys %lanes) {
    my $type;
    if ($lanes{$lane} == 1) {
	$type='single';
    } elsif ($lanes{$lane} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my $infilename=$lane.'.'.$type.'.unique.gtf.gz';
    print STDERR "Processing $infilename\n";

    unless (-r $infilename) {
	die "$infilename is not readable\n";
    }

    build_distribution($infilename,
		       \%distribution,
		       $lane);

    print STDERR "done\n";

    # Generate the data for R and build the plots
    data_4_R($infilename,
	     \%lengths,
	     $lane,
	     $disttablefh);
    
}
close($disttablefh);

build_table($outfile,
	    \%distribution);

exit;

sub build_table {
    my $graph=shift;
    my $dist=shift;

    my $stats_file=$graph;
    unless ($stats_file=~/\.txt$/) {
	$stats_file.='.txt';
    }
    my $tmpfh=get_fh($stats_file,1);
    foreach my $filename (keys %{$dist}) {
	foreach my $chr (keys %{$dist->{$filename}}) {
	    my $chr_id=$chr;
	    print $tmpfh join("\t",
			      $filename,
			      $chr_id,
			      $dist->{$filename}->{$chr}),"\n";
	}
    }
    # Close the file to make sure buffer is flushed
    close($tmpfh);
} 

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}++;
    }

    return(\%lanes);
}

sub build_distribution {
    my $file=shift;
    my $dist=shift;
    my $lane=shift;

    my $infh=get_fh($file);

    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);
	my $filename=$file;
	$filename=~s/.*\///;
	$dist->{$lane}->{$line[0]}++;
    }

    close($infh);
}

# This subroutine will get 5 filehandles one which will contain the stats
# for all the reads and another four that will contain subsets of the data
sub get_file_handles {
    my %files=("read_dist.1_99.$$.stats" => '',
	       "read_dist.100_999.$$.stats" => '',
	       "read_dist.1000_2999.$$.stats" => '',
	       "read_dist.3000_5999.$$.stats" => '',
	       "read_dist.6000_8999.$$.stats" => '',
	       "read_dist.9000_n.$$.stats" => '',
	       "read_dist.all.$$.stats" => '');
    foreach my $file (keys %files) {
	my $fh=get_fh($file,1);
	$files{$file}=$fh;
    }
    return(\%files);
}

# This will just print each value in one of the files depending on the
# size of the transcript (length)
sub breakdown_by_size {
    my $files=shift;
    my $value=shift;
    my $length=shift;
    my $fh;

    if ($length < 100) {
	$fh=$files->{"read_dist.1_99.$$.stats"};
    } elsif ($length < 1000) {
	$fh=$files->{"read_dist.100_999.$$.stats"};
    } elsif ($length < 3000) {
	$fh=$files->{"read_dist.1000_2999.$$.stats"};
    } elsif ($length < 6000) {
	$fh=$files->{"read_dist.3000_5999.$$.stats"};
    } elsif ($length < 9000) {
	$fh=$files->{"read_dist.6000_8999.$$.stats"};
    } else {
	$fh=$files->{"read_dist.9000_n.$$.stats"};
    }

    return($fh);
}

# Get the size category to which the transcript belongs
sub get_size_cat {
    my $length=shift;
    my $cat;

    if ($length < 100) {
	$cat='1_99';
    } elsif ($length < 1000) {
	$cat='100_999';
    } elsif ($length < 3000) {
	$cat='1000_2999';
    } elsif ($length < 6000) {
	$cat='3000_5999';
    } elsif ($length < 9000) {
	$cat='6000_8999';
    } else {
	$cat='9000_n';
    }

    return($cat);
}

# This sub will get each nt of each of the reads an asign it a position along
# the transcript adding it to a plot
sub data_4_R {
    my $readsfile=shift;
    my $lengths=shift;
    my $lane=shift;
    my $disttablefh=shift;
    my $fh=get_fh($readsfile);

    my %files=%{get_file_handles()};
    my %length_dist;

    while (my $line=<$fh>) {
	chomp($line);
	my @line=split("\t",$line);
	my ($read_id,$start,$end)=@line[0,3,4];
	my $length=$end - $start + 1;
	my $trans_length=$lengths->{$read_id};
	for (my $i=0;$i<$length;$i+=$offset) {
	    unless ($trans_length) {
		die "No length for $line\n";
	    }
	    my $dist_value=sprintf "%.2f",($start / $trans_length) * 100;
	    my $fh_all=$files{"read_dist.all.$$.stats"};
	    print $fh_all join("\t",
			       $dist_value,
			       $trans_length),"\n";
	    $start+=$offset;
	    my $fh_group=breakdown_by_size(\%files,
					   $dist_value,
					   $trans_length);

	    my $size_cat=get_size_cat($trans_length);
	    my $round=int($dist_value + 0.5);
	    $length_dist{$size_cat}{$round}++;
	    $length_dist{'0_All'}{$round}++;

	    print $fh_group join("\t",
				 $dist_value,
				 $trans_length),"\n";	    
	}
    }
    foreach my $file (keys %files) {
	close($file);
    }

    # Print out the breakdown table
    foreach my $cat (keys %length_dist) {
	foreach my $value (sort {$a <=> $b} keys %{$length_dist{$cat}}) {
	    my ($start)=split('_',$cat);
	    my $key=$cat;
	    if ($cat eq '0_All') {
		$key='All';
	    }
	    print $disttablefh join("\t",
				    $key,
				    $start,
				    $value,
				    $length_dist{$cat}{$value},
				    $lane),"\n";
	}
    }

    # build the plots using R
    # We will build one figure with the distribution of all the reads and one
    # with the breakdown for each of the files
    
    $readsfile=~s/.*\///;
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;

    # Print the Breakdown of the read distribution
    $r_string.="postscript(\"$readsfile.breakdown.stats.ps\")\n";
    $r_string.="par(oma=c(2,0,2,0))\n";
    $r_string.="par(mfrow=c(2,3))\n";
    foreach my $statsfile ("read_dist.1_99.$$.stats",
			   "read_dist.100_999.$$.stats",
			   "read_dist.1000_2999.$$.stats",
			   "read_dist.3000_5999.$$.stats",
			   "read_dist.6000_8999.$$.stats",
			   "read_dist.9000_n.$$.stats") {
	# First check the number of points, if there are more than 100M we will
	# skip this step, as it will be problematic
	my $count=`wc -l $statsfile`;
	chomp($count);
	($count)=split(/\s+/,$count);
	$count=~s/\s*//g;

	if ($count > 100000000) {
	    warn "File $statsfile is too big. Skipping...\n";
	    next;
	}
	my $file=$statsfile;
	$file=~s/\.$$\.stats//;
	my $range=(split(/\./,$file))[1];
	$range=~s/_/-/;
	$r_string.="stats<-read.table(\"$statsfile\",sep=\"\t\")\n";
	$r_string.='hist(stats$V1,col="cyan",';
	$r_string.="main=\"Transcript length:$range\",";
	$r_string.="xlab=\"Transcript length\",ylab=\"Read coverage\")\n";
    }
    $r_string.="title(main=\"$species $project $lane Transcript read distribution\",sub=\"File:$readsfile\",outer=T)\n";
    $r_string.="dev.off()\n";

    # Print the total read distribution
    # First check the number of points, if there are more than 100M we will skip
    # this step, as it will be problematic
    my $statsfile="read_dist.all.$$.stats";
    my $count=`wc -l $statsfile`;
    chomp($count);
    ($count)=split(/\s+/,$count);
    $count=~s/\s*//g;

    if ($count > 100000000) {
	warn "$statsfile is too big to create all distribution. Skipping...\n";
    } else {
	my $range='All';
	$range=~s/read_dist.//;
	$r_string.="postscript(\"$readsfile.all.stats.ps\")\n";
	$r_string.="stats<-read.table(\"$statsfile\",sep=\"\t\")\n";
	$r_string.='hist(stats$V1,col="cyan",';
	$r_string.="main=\"$species $project $lane $range transcript read distribution\",";
	$r_string.="xlab=\"Transcript length\",ylab=\"Read coverage\")\n";
	$r_string.="dev.off()\n";
    }

    print $exec_fh $r_string;
    close($exec_fh);
	
    # execute the R file and remove it
    my $command="R --vanilla --quiet < $execution_file";
    system($command);
    $command="rm $execution_file";
    system($command);

    # remove the input files
    foreach my $file (keys %files) {
	$command="rm $file";
	system($command);
    }
}

# This script will extract the reads from the mapping file and assign them to
# bins depending on the trasncript they belong to
sub assign_reads {
    my $file=shift;
    my $bins=shift;
    my $fh=get_fh($file);
    my @dist;
    my $count;

    while (my $line=<$fh>) {
	chomp($line);
	my @line=split("\t",$line);
	my ($read_id,$start,$end)=@line[0,3,4];
	my $length=$end - $start + 1;
	# Skip short transcripts
	if ($length > $bins->{$read_id}->[$bin_no/2]) {
	    next;
	}
	for (my $i=0;$i<$bin_no;$i++) {
	    my $limit=$bins->{$read_id}->[$i];
	    if ($start < $limit) {
		while ($start < $end) {
		    my $frag=$limit - $start;
		    $start+=$frag;
		    $dist[$i]+= ($frag / $length);
		    $i++;
		    $limit=$bins->{$read_id}->[$i];
		}
		last;
	    }
	}
	$count++;
    }
    close ($fh);
    foreach my $frag (@dist) {
	printf "%.2f\n", ($frag / 1);
    }

}


sub get_transcript_lengths {
    my $file=shift;
    my $fh=get_fh($file);

    my %lengths;
    while (my $line=<$fh>) {
	chomp($line);
	my ($trans_id,$trans_length)=split("\t",$line);
	$lengths{$trans_id}=$trans_length;
    }
    close($fh);

    return(\%lengths)
}

