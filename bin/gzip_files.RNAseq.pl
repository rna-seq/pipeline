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
