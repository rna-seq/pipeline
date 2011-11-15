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
# This script will take a gtf file as an argument and a list of types, and it
# will print each of the supplied types into their respective files.

use RNAseq_pipeline3 qw(get_fh gff_line_check);
use RNAseq_pipeline_settings3 qw(read_config_file);

# Declare some variables
my $prefix;
my $annotation;
my $projectdir;

my %types=('gene' => 'genome',
	   'exon' => 'exons');

# Collect options
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$annotation=$options{'ANNOTATION'};
$projectdir=$options{'PROJECT'};

unless(-r $annotation) {
    die "$annotation is not readable\n";
}

# Get the filehandles that are needed depending on the required types
my %files=%{get_fhs($prefix,
		    $projectdir,
		    \%types)};

my %counts;

my $infh=get_fh($annotation);
print STDERR "Reading Annotation from $annotation\n";
while (my $line=<$infh>) {
    if ($line=~/^#/) {
	next;
    }
    # This will process the line to make sure it complies with the gtf standard
    my ($parsed,$type)=gff_line_check($line);

    if ($files{$type}) {
	$counts{$type}++;
	my $outfh=$files{$type};
	print $outfh $parsed,"\n";
    }
}
close($infh);
print STDERR "done\n";

# close output files
foreach my $file (keys %files) {
    close($files{$file});
}

# Print the results
foreach my $type (keys %counts) {
    print join("\t",
	       $type,
	       $counts{$type}),"\n";
}
	      

exit;

sub get_fhs {
    my $prefix=shift;
    my $dir=shift;
    my $types=shift;
    my %files;

    foreach my $type (keys %{$types}) {
	my $filename;
	$filename.=$dir.'/'.$types{$type}.'/'.$prefix.'.'.$type.'.gtf';
	$files{$type}=get_fh($filename,1);
	print STDERR "Printing $type to $filename\n";
    }
    return(\%files);
}
