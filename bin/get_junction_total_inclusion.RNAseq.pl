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
# This script will take a gff file with the mappings of reads on the junctions,
# and a read length, and it will count the number of reads completely included
# in the junction (this is where the length of the mapped read is equal to the
# read length)

# Added the stranded option
### TO DO must rewrite this

use RNAseq_pipeline_settings3 qw(read_config_file);
use RNAseq_pipeline3 qw(get_fh);

# Get the options
my %options=%{read_config_file()};
my $readlength=$options{'READLENGTH'};
my $stranded=$options{'STRANDED'};

# Get the file
my $infile=shift;
if ($infile) {
    print STDERR "Processing $infile\n";
} else {
    die "No input file\n";
}

# Count the junction occurrences. Here they are treated as chromosomes,
# as each junction is a separate sequence. (the chromosome field in the gtf
# file contains the junction ID)
my %junctions=%{get_junctions_count($infile,
				    $readlength,
				    $stranded)};

foreach my $junction (keys %junctions) {
    print join("\t",
	       $junction,
	       $junctions{$junction}),"\n";
}

exit;

sub get_junctions_count {
    my $infile=shift;
    my $readlength=shift;
    my $stranded=shift;
    
    my $infh=get_fh($infile);
    my %junctions;
    my %stranded;
    my $short=0;
    
    while (my $line=<$infh>) {
	# here we are processing pairs of lines, as each junction will
	# correspond to one pair, so we have to take the next line too
	chomp($line);
	my @line=split("\t",$line);
	my $chr=$line[0];
	my $length1=$line[4] - $line[3] + 1;

	my $line2=<$infh>;
	chomp($line2);
	my @line2=split("\t",$line2);
	my $length2=$line2[4] - $line2[3] + 1;
	my $length=$length1 + $length2;

	my $site=join('_',$chr,$line[4],'splice',$line2[3]);
	my $site_stranded=join('_',$chr,$line[4],'splice',$line2[3],$line[6]);
	
	if ($length == $readlength) {		
	    $junctions{$site}++;
	    $stranded{$site_stranded}++;
	} elsif ($length > $readlength) {
	    warn "Mapped read $length is longer then $readlength for $chr\n";
	    $junctions{$site}++;
	    $stranded{$site_stranded}++;
	} else {
	    $short++;
	}
    }
    close($infh);
    
    if ($short) {
	print STDERR $short,"\tRead maps were shorter than $readlength\n";
    }

    if ($stranded) {
	return(\%stranded);
    } else {
	return(\%junctions);
    }
}
