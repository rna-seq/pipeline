#!/soft/bin/perl
# DGK

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
# This script should take a fasta or fastq file and it will remove some
# invalid characters from the id

use Pod::Usage;
use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh);

# Declare some variables
my $help;
my $man;
my $infile;

# Get the command line options
GetOptions('help|h' => \$help,
	   'man|m' => \$man);

# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

$infile=shift;
pod2usage("No input file") unless ($infile);

my $infh=get_fh($infile);
while (my $line=<$infh>) {
    if ($line=~/^[+>@]/o) {
	chomp($line);
	if ($line=~s/\|p?1$/\/1/o) {
	    print $line,"\n";
	} elsif ($line=~s/\|p?2$/\/2/o) {
	    print $line,"\n";
	}
    } else {
	print $line;
    }
}
close($infh);

exit;

__END__

=head1 NAME
    
    clean_read_ids.preprocess.pl
    
=head1 SYNOPSIS
    
sample [hm] input 
    
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
    
    This script should rempve the p characters from the endings of the fastq
    files (p1 and p2)

=cut
