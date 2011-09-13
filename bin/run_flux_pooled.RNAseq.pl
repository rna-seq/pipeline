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

# Objective
# This script will create the necessary commands for running the flux capacitor
# within the framework of the opipeline
# do flux.sh -s $file.bed -r fluxout/gencode_data.rel2b.exons_sorted.gtf --tmp tmp -n fluxoutpaired --pair; done


use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_groups);
use Tools::Flux qw(build_parameter_file run_flux);

my $localdir;
my $exonfile;
my $prefix;
my $bindir;
my $sorted=0; # In principle the files will not be sorted, so they will need
              # to be sorted by the flux
my $tmpdir;
my $debug=0;

my %options=%{read_config_file()};
$localdir=$options{'LOCALDIR'};
$prefix=$options{'PREFIX'};
$exonfile=$options{'EXONDIR'}.'/'.$prefix.'.exon.gtf';
$bindir=$options{'BIN'};
$tmpdir=$options{'LOCALDIR'};

# Get a log file
my $log_fh=get_log_fh('run_flux_pooled.RNAseq.log',
		      $debug);

# First get a list of the bed files we are going to process
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};
my %groups=%{get_groups()};

# get the necesary directories we are going to use in order to build the command
# line for the capacitor

# split the annotation
my %annot;
unless (-r $exonfile) {
    die "Can't read $exonfile\n";
}

# Run the capacitor for each of the files that we have and move the output to
# a file with the adequate name
foreach my $group (keys %groups) {
    my %reads;
    my $type;

    # Check if the file corresponding to this group exists already
    my $fluxfile=$group.'.flux.gtf';
    if (-r $fluxfile) {
	print STDERR $fluxfile,"\tPresent. Skipping...\n";
	next;
    }

    # Determine if paired
    foreach my $lane (keys %{$groups{$group}}) {
	if (@{$lanes{$lane}} == 1) {
	    $type='single';
	} elsif (@{$lanes{$lane}} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type\n";
	}
    }



    # Check if we have the combined bed file an if not create it
    my $infilename=$localdir.'/'.$group.'.combined.bed';
    if (-r $infilename) {
	print $log_fh "Processing $infilename\n";
    } else {
	print $log_fh "Building $infilename\n";
	my %files;
	foreach my $lane (@{$groups{$group}}) {
	    # Get the files corresponding to the both halves of the reads
	    my $infilename=$localdir.'/'.$lane.'.combined.bed';
	    if (-r $infilename) {
		print $log_fh "Processing $infilename\n";
		$files{$infilename}++;
	    } else {
		die "Can't read $infilename\n";
	    }
	}

	# Concatenate the different files
	my $command='cat '.join(' ',keys %files)." > $infilename";
	run_system_command($command);
	
    }
    
    my %results;
    my $outfile=run_flux($localdir,
			 $type,
			 $infilename,
			 $exonfile,
			 $bindir,
			 $tmpdir,
			 $sorted);
	
    if ($outfile) {
	my $command="mv $outfile $fluxfile";
	run_system_command($command);
    } else {
	die "Problem generating $outfile\n";
    }
}

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	push @{$lanes{$files->{$file}->[0]}},$files->{$file}->[1];
    }

    return(\%lanes);
}
