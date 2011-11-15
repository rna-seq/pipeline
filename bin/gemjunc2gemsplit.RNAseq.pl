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

use Pod::Usage;
use Getopt::Long;
use RNAseq_pipeline3 ('get_fh');
use RNAseq_GEM3 ('parse_gem_line','parse_gem_hit','gem_2_coords',
		 'get_junction_coords','coords2splitcoords');

# Declare some variables
my $help;
my $man;
my $infile;
my $outfile;
my $mismatches=2;

# Get the command line options
GetOptions('help|h' => \$help,
	   'man|m' => \$man,
	   'in|i=s' => \$infile,
	   'out|o=s' => \$outfile);

# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Check if we have input and output files
pod2usage("I cant find $infile. In need an input file") unless ($infile && (-r $infile));
pod2usage("In need an output file") unless ($outfile);

if (-r $outfile ||
    -r $outfile.'.gz') {
    print STDERR $outfile,"\tis present. Skipping...\n";
    exit;
}

# Process the input junctions file
print STDERR "Parsing junction maps fron $infile...\n";
my $infh=get_fh($infile);
my $outfh=get_fh($outfile,1);


while (my $line=<$infh>) {
    my %line=%{parse_gem_line($line)};
    my @coords=@{gem_2_coords(\%line)};

    my $no_hits;
    if ($line{'hits'} eq '-') {
	$no_hits=$line{'matches'};
    }

    my %hits;
    my $max_mismatch;
    foreach my $coords (@coords) {
	my ($hitstring,$mismatchno)=coords2splitcoords($coords);
	unless($hitstring) {next;}
	# The mapping to the junctions will generate redundant mappings
	# when there are alternative exons with different starts or ends,
	# this should remove this redundancy
	# if we have a mix of split and unsplit hits we want to keep the best
	# type discarding all the hits corresponding to the other type. If they
	# are equally good we keep the unsplit hit
	$hits{$hitstring}->[0]++;
	if ($hits{$hitstring}->[1] &&
	    ($mismatchno < $hits{$hitstring}->[1])) {
	    $hits{$hitstring}->[1]=$mismatchno;
	} else {
	    $hits{$hitstring}->[1]=$mismatchno;
	}

	# Set the hit type
	my $hittype='full';
	if ($hitstring=~/~/) {
	    $hittype='split';
	}
	$hits{$hitstring}->[2]=$hittype;
    }

    # Collect all the hits as a string
    my @hits;
    my @mismatches;
    # Initialize the mismatches
    for (my $i=0;$i<=$mismatches;$i++) {
	$mismatches[$i]=0;
    }
    my $oldtype;
    my $oldmismatches;
    foreach my $hit (sort {$hits{$a}->[1] <=> $hits{$b}->[1]} keys %hits) {
	# make sure we do not mix types
	my $type=$hits{$hit}->[2];
	if ($oldtype) {
	    if ($oldtype eq $type) {
		push @hits,$hit;
		$mismatches[$hits{$hit}->[1]]++;
	    } elsif ($oldmismatches == $hits{$hit}->[1]) {
		if ($type eq 'full') {
		    die "Sorting problem\n";
		} else {
		    last;
		}
	    } else {
		warn "There is a problem\n";
	    }
	} else {
	    $oldtype=$type;
	    $oldmismatches=$hits{$hit}->[1];
	    push @hits,$hit;
	    $mismatches[$hits{$hit}->[1]]++;
	}
    }
    my $allhits=join(',',@hits);

    # Build the mismatch string
    foreach my $hits (@mismatches) {
	unless ($hits) {
	    $hits=0;
	}
    }
    my $mismatchstring=join(':',@mismatches);

    unless ($allhits) {
	# Here we need to distinguish between the cases with too many hits and
	# those with none.
	$allhits='-';
	if ($no_hits) {
	    $mismatchstring=$no_hits;
	}
    }

    # Print the complete read line
    my $quals=$line{'qual'};
    unless($quals) {
	$quals='B' x length($line{'seq'});
    }
    print $outfh join("\t",
		      $line{'id'},
		      $line{'seq'},
		      $quals,
		      $mismatchstring,
		      $allhits),"\n";
}
close($infh);
close($outfh);

exit;

__END__

=head1 NAME
    
    gemjunc2gemsplit.pl
    
=head1 SYNOPSIS
    
sample [options] 
    
  Options:
    -help            brief help message
    -man             full documentation
        
=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.
    
=back
    
=head1 DESCRIPTION
    
    gemjunc2gemsplit.pl
    The objective of this script is to take a junctions mapping file that
    contains the junction ids as chr start, end and strand and it will
    transform it into a split-mapping file with genome coordinates

=cut
