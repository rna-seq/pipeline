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
# This script should take the gff files resulting from the split-mappings and
# the junctions overlap total files, and build a list of junctions that combine
# those from the annotation with those from the splitmapping.

# Load some 
use RNAseq_pipeline3 ('get_fh','parse_gff_line');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh',
			       'get_junction_type_sub');
use Getopt::Long;

# Declare some variables
my $threshold=1;
my $splitdir;
my $junctionsdir;
my $projectid;
my $prefix;
my $stranded;
my $file_list;
my $parallel='default';
my $paralleltmp;
my $junctable;

GetOptions('threshold|t=i' => \$threshold,
    );

# Read the options file
my %options=%{read_config_file()};
$splitdir=$options{'SPLITMAPDIR'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$projectid=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$stranded=$options{'STRANDED'};
$file_list=$options{'FILELIST'};
$paralleltmp=$options{'PROJECT'}.'/'.$options{'LOCALPARALLEL'};
$junctable=$options{'JUNCTIONSTABLE'};
if ($options{'PARALLEL'}) {
    $parallel='parallel';
}

unless ($paralleltmp) {
    $paralleltmp=$options{'PROJECT'}.'/work';
}

# Get the requred subs
*junction_type=get_junction_type_sub($junctable);
*get_start_exons=get_exons_from_start_sub();
*get_end_exons=get_exons_from_end_sub();

# First get the files we are going to analyze
my %files=%{read_file_list($file_list)};

# Get for each lane the files corresponding to the genome mapping
my %lane_files=%{get_lane_files(\%files)};
my %groups=%{get_groups(\%files)};


my %coverage;

foreach my $group (keys %groups) {
    my %junctions;

    foreach my $pair (keys %lane_files) {
	# Find the coverage granted by genome mapping;
	foreach my $lane (keys %{$lane_files{$pair}}) {
	    # Process the junctions file
	    my $juncfile=$junctionsdir.'/'.$lane.'.single.unique.overlap.total';
	    process_junctions($juncfile,
			      $pair,
			      \%junctions);
	    
	    my $splitfile=$splitdir.'/'.$lane.'.single.unique.gtf.gz';
	    process_splits($splitfile,
			   $pair,
			   \%junctions);
	    
	    # process the split file

	}
    }

    # Print out the information
    my $outfile=$junctionsdir.'/'.$group.'.paired.all.junc.class.pooled';
    print STDERR "Processing $outfile\n";

    # Check if file exists already
    my $overfinalname=$outfile;
    $overfinalname=~s/.gz//;

    my $outfh=get_fh($outfile,1);
    foreach my $junc (keys %junctions) {
	print $outfh join("\t",
			  @{$junctions{$junc}}),"\n";
	
    }
    close($outfh);
    print STDERR "$group\t$outfile done\n";
    # Clear the memory
    %junctions=();
}
    
exit;

sub process_junctions {
    my $infile=shift;
    my $pair=shift;
    my $junctions=shift;

    print STDERR "Reading $infile\n";
    my $infh=get_fh($infile);
    while (my $line=<$infh>) {
	chomp($line);
	
	my ($junc,$hits)=split("\t",$line);
	# Get the junction type
	my $type=junction_type($junc);
	$junc=~s/_splice//;
	my @junc=split('_',$junc);
	my $end=pop(@junc);
	my $start=pop(@junc);
	my $chr=join('_',@junc);

	# get the start exons
	my $start_exons=get_start_exons($chr,
					$start);
	# get the end exons
	my $end_exons=get_end_exons($chr,
				    $end);

	if ($junctions->{$junc}) {
	    $junctions->{$junc}->[5]+=$hits;
	} else {
	    $junctions->{$junc}=[$chr,$start,$chr,$end,$type,$hits,
				 $start_exons,$end_exons,$pair];
	}

    }
}

sub process_splits {
    my $infile=shift;
    my $pair=shift;
    my $splits=shift;

    print STDERR "Reading $infile\n";
    my $infh=get_fh($infile);
    my $line_no=0;
    my ($start,$end,$chr1,$chr2);
    while (my $line=<$infh>) {
	$line_no++;
	chomp($line);

	my %line=%{parse_gff_line($line)};

	# Check if it is an end line or a start line
	if ($line_no % 2) {
	    # splice start line 
	    $start=$line{'end'};
	    $chr1=$line{'chr'};
	} else {
	    # splice end line
	    $end=$line{'start'};
	    $chr2=$line{'chr'};
	    # get the start exons
	    my $start_exons=get_start_exons($chr1,
					    $start);
	    # get the end exons
	    my $end_exons=get_end_exons($chr2,
					$end);

	    my ($type)='split';
	    if ($chr1 ne $chr2) {
		$type='split_interchro';
	    }
	    my $split_id=join("\t",$chr1,$chr2,$start,$end);

	    if ($splits->{$split_id}) {
		$splits->{$split_id}->[5]++;
	    } else {
		$splits->{$split_id}=[$chr1,$start,$chr2,$end,$type,0,
				      $start_exons,$end_exons,$pair];
	    }
	    ($chr1,$chr2,$start,$end)=('','','','');
	}
    }
    close($infh);
}

sub get_exons_from_end_sub {
    my %options=%{read_config_file()};
    my $dbh=get_dbh(1);
    my $table=$options{'JUNCTIONSTABLE'};

    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT exon2_id ';
    $query.="FROM $table ";
    $query.='WHERE chr = ? AND end = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $chr=shift;
	my $end=shift;

	my $juncid=join('_',$chr,$end);

	unless ($cache{$juncid}) {
	    $count=$sth->execute($chr,$end);

	    if ($count < 1) {
		$cache{$juncid}='-';
	    } else {
		my @exons;
		while (my ($exon)=$sth->fetchrow_array()) {
		    push @exons,$exon;
		}
		$cache{$juncid}=join(';',@exons);
	    }
	}
	return($cache{$juncid});
    };
    return($subroutine);
}

sub get_exons_from_start_sub {
    my %options=%{read_config_file()};
    my $dbh=get_dbh(1);
    my $table=$options{'JUNCTIONSTABLE'};

    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT exon1_id ';
    $query.="FROM $table ";
    $query.='WHERE chr = ? AND start = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $chr=shift;
	my $start=shift;

	my $juncid=join('_',$chr,$start);

	unless ($cache{$juncid}) {
	    $count=$sth->execute($chr,$start);

	    if ($count < 1) {
		$cache{$juncid}='-';
	    } else {
		my @exons;
		while (my ($exon)=$sth->fetchrow_array()) {
		    push @exons,$exon;
		}
		$cache{$juncid}=join(';',@exons);
	    }
	}
	return($cache{$juncid});
    };
    return($subroutine);
}

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	$lane_files{$pair}{$lane}=1;
    }

    return(\%lane_files);
}

sub get_groups {
    my $files=shift;
    my %groups;
    
    foreach my $file (keys %{$files}) {
	my $group=$files->{$file}->[2] || 'All';
	push @{$groups{$group}},$files->{$file}->[0];
    }

    return(\%groups);
}
