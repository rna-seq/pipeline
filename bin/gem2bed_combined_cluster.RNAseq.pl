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

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh');
use RNAseq_GEM3 ('parse_gem_line','gem_2_coords','coords2gtf',
		 'coords2bedSimple','coords2bedJunctions',
		 'get_junction_coords');

# Get command line options
my $read_length;
my $file_list;
my $localdir;

my %options=%{read_config_file()};
$file_list=$options{'FILELIST'};
$read_length=$options{'READLENGTH'};
$localdir=$options{'LOCALDIR'};

# Get some subroutines
my $dbh=get_dbh();
my $gen_map_table=$options{'PREFIX'}.'_genome_mapping';
my $junc_map_table=$options{'PREFIX'}.'_junctions_mapping';
*get_gen_file=get_files_from_table_sub($dbh,
				       $gen_map_table);
*get_junc_file=get_files_from_table_sub($dbh,
					$junc_map_table);

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
	print STDERR "done\n";
    }
    close($bedfh);
}

exit;

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $junctionfile=get_junc_file($lane);
	my $genomefile=get_gen_file($lane);

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
