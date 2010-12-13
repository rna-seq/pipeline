#!/soft/bin/perl
# DGK

use strict;
use warnings;

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
