#!/soft/bin/perl

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

# gzip a list of files with the same ending if unzipped 

my $ending=shift;

my @files=`ls`;
my @unzipped;

# Taint check
if ($ending =~ /^([-\@\w.]+)$/) {
    $ending = $1;
} else {
    die "The provided ending containts potentially unsafe characters\n";
}

foreach my $file (@files) {
    chomp($file);
    if ($file=~/$ending$/) {
	push @unzipped, $file;
    }
}

if (@unzipped) {
    my $command="gzip -f -7 ";
    $command.=join(' ',@unzipped);
    print STDERR "Executing:\t$command\n";
    system($command);
} else {
    print STDERR "All files are gzipped\n";
}

exit;
