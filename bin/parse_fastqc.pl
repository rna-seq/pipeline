#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take the output of fastqc and parse it organizing it into
# 4 tables:
# Per base quality scores
# Per base nucleotide count
# Per file sequence counts
# General statistics for the sequences (tests by fastqc, etc...)

use Pod::Usage;
use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use Tools::FastQC qw(parse_fastqc_report);

# Declare some variables
my $help;
my $man;
my $localdir;
my $prefix;
my %report;

# Get the command line options
GetOptions('help|h' => \$help,
	   'man|m' => \$man);

# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Read configuration file
my %options=%{read_config_file()};
$localdir=$options{'LOCALDIR'};
$prefix=$options{'PREFIX'};

# Get a list of the files that should be present if fastqc is run
my %files=%{read_file_list()};
foreach my $file (keys %files) {
    print STDERR $file,"\n";
    print STDERR $files{$file}->[1],"\n";
    my $fastqcdir=$file;
    $fastqcdir=~s/(.fa.*$|.bam)/_fastqc/;
    print STDERR $fastqcdir,"\n";
    $fastqcdir=$localdir.'/'.$fastqcdir;
    if (-s $fastqcdir) {
	print "Processing data from $fastqcdir\n";
    } else {
	die "Unable to find fastqc data for $file\n";
    }
    my $report=parse_fastqc_report($fastqcdir);
    $report{$file}=$report;
}

foreach my $file (keys %report) {
    my $dataset=$files{$file}->[1];
    my %summary=%{$report{$file}};
    # print the basic statistics file
    my $statsfn=$prefix.'_basic_stats.txt';
    my $statsfh=get_fh($statsfn,1);
    foreach my $key (keys %summary) {
	print $statsfh join("\t",
			    $dataset,
			    @{$summary{$key}{'status'}}),"\n";
    }
    close($statsfh);
    
    # print the qualities statistics file
    my $qtfn=$prefix.'_qualitiespos.txt';
    my $qtfh=get_fh($qtfn,1);
    my $key='Perbasesequencequality';
    foreach my $line (@{$summary{$key}{'values'}}) {
	print $qtfh join("\t",
			    $dataset,
			    @{$line}[0,1]),"\n";
    }
    close($qtfh);

    # print the nucleotides statistics file
    my $ntfn=$prefix.'_ambiguous.txt';
    my $ntfh=get_fh($ntfn,1);
    my $key2='PerbaseNcontent';
    my $key3='Perbasesequencecontent';
    for (my $i=0; $i < @{$summary{$key2}{'values'}} ;$i++) {
	my $line1=$summary{$key2}{'values'}->[$i];
	my $line2=$summary{$key3}{'values'}->[$i];
	unless ($line1->[0] eq $line2->[0]) {
	    die "Problem with nt number\n";
	}
	print $ntfh join("\t",
			 $dataset,
			 @{$line1}[0,1],
			 @{$line2}[2,4,1,3]),"\n";
    }
    close($ntfh);
}

exit;

__END__

=head1 NAME
    
    ?
    
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
    
    This program is not documented yet, sorry

=cut
