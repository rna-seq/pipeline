#!/soft/bin/perl
# DGK

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
# This script should take three files one corresponding to the genome mapping,
# one to the junctions mapping and one to the split mapping, and it will
# combine them into a full bed file
# with one entry per read, this means the junctions mappings will be split, the
# same as the split reads.

use RNAseq_GEM3 ('parse_gem_line','gem_2_coords','coords2gtf',
		 'coords2bedSimple','coords2bedJunctions',
		 'get_junction_coords');
use RNAseq_pipeline3 qw(get_fh);
use Getopt::Long;

# Declare variables
my $genome_map;
my $junction_map;
my $split_map;

GetOptions('genome|g=s' => \$genome_map,
	   'junctions|j=s' => \$junction_map,
	   'split|s=s' => \$split_map);

unless ($genome_map || $junction_map || $split_map) {
    die "at least one mapping file should be provided\n";
}

if ($genome_map) {
    print STDERR "Parsing genome reads..";
    my $infh=get_fh($genome_map);
    while (my $line=<$infh>) {
	my %line=%{parse_gem_line($line)};
	my @coords=@{gem_2_coords(\%line)};

	foreach my $coords (@coords) {
	    my $bed=coords2bedSimple($coords,
				     $genome_map);
	    print $bed,"\n";
	}
    }
    close($infh);
    print STDERR "done\n";
}

if ($junction_map) {
    print STDERR "Parsing junction maps...";
    my $infh=get_fh($junction_map);
    my $interchrom=0;
    while (my $line=<$infh>) {
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
		    print $oldbed,"\n";
		    $oldbed=$bed;
		} else {
		    $oldbed=$bed;
		}
	    }
	    if ($oldbed) {
		print $oldbed,"\n";
	    }
	}
    }
    close($infh);
    print STDERR "done\n";
}

exit;
