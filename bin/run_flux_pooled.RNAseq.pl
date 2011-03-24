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
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

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
my %groups=%{get_groups(\%files)};

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
    foreach my $lane (@{$groups{$group}}) {
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

sub build_parameter_file {
    my $annot=shift;
    my $maps=shift;
    my $type=shift;
    my $sorted_maps=shift;
    my $output=shift;
    my $tmpdir=shift;

    # Set the name for the parameter file
    my  $parameterfn=$output;
    $parameterfn=~s/.gtf$/.params/;

    # Check the type
    unless ($type=~/(Paired|Single)/i) {
	die "Incorrect read type $type\n";
    }
    # Make sure we have the correct casing (In case the flux is strict with
    # this)
    $type=lc($type);
    $type=ucfirst($type);
    my $descriptor=uc($type);

    if ($descriptor eq 'SINGLE') {
	$descriptor='SIMPLE';
    }

    # Print the entries of the parameter file
    my $paramfh=get_fh($parameterfn,1);
    print $paramfh join(" ",
			'ANNOTATION_FILE',
			$annot),"\n";
    print $paramfh join(" ",
			'MAPPING_FILE',
			$maps),"\n";
    print $paramfh join(" ",
			'SORTED_MAPPINGS_FILE',
			$sorted_maps),"\n";
    ### WARNING
    # if  this flag is not set the input annotation and read files will be
    # deleted (WORK_LOCAL)
    print $paramfh join(" = ",
			'WORK_LOCAL',
			'1'),"\n";
    print $paramfh join(" ",
			'STDOUT_FILE',
			$output),"\n";
    print $paramfh join(" ",
			'TMP_DIR',
			$tmpdir),"\n";
    print $paramfh join(" ",
			'READ_DESCRIPTOR',
			$descriptor),"\n";
    print $paramfh join(" ",
			'ANNOTATION_MAPPING',
			$type),"\n";

    close($paramfh);

    return($parameterfn)
}

sub run_flux {
    my $localdir=shift;
    my $type=shift;
    my $maps=shift;
    my $annot=shift;
    my $bindir=shift;
    my $tmpdir=shift;
    my $sorted=shift;

    # set the name of the outfile
    my $fileroot1=$annot;
    my $fileroot2=$maps;

    # We have to make sure we are only removing up to the first ending '.',
    # As we could have naming problems if not
    ### TO DO
    # Think of a more stable way of doing this that allows for '.' in the
    #  path name
    $fileroot1=~s/(\.\w+?)?$//;
    $fileroot2=~s/(\.\w+?)?$//;
    $fileroot2=~s/.*\///;

    my $outfile=join('_',
		     $fileroot1,
		     'flux',
		     $fileroot2).'.gtf';
    
    my $sorted_maps;
    if ($sorted) {
	$sorted_maps=$maps;
    } else {
	$sorted_maps=$maps.'.sorted';
    }

    # Check if the output exists already and if so skip the creation
    # /home/dgonzalez/bin/testPipeline3/work/chr7_sorted__chr7.gtf
    if (-r $outfile) {	
	print STDERR $outfile,"\tIs present. Skipping...\n";
	return($outfile);
    }

    # Add a subdirectory to the tmpdir
    $tmpdir.='/tmp';
    # create this second localdir if it doesn't exist
    my $command="mkdir -p $tmpdir";
    print STDERR "Executing: $command\n";
    system($command);

    my $parameterfile=build_parameter_file($annot,
					   $maps,
					   $type,
					   $sorted_maps,
					   $outfile,
					   $tmpdir);

    # build the flux command
    $command="$bindir/flux.sh $parameterfile -u";
    print STDERR $command,"\n";
    system($command);

    if (-r $outfile) {	
	print STDERR $outfile,"\tBuilt\n";
	# Clean up by removing the parameter file
	# The read file seems to be cleared by the flux
	$command="rm $parameterfile $maps";
	print STDERR "Executing: $command\n";
	system($command);
    } else {
	die "Problem building $outfile\n";
    }

    return($outfile);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	push @{$lanes{$files->{$file}->[0]}},$files->{$file}->[1];
    }

    return(\%lanes);
}

sub get_groups {
    my $files=shift;
    my %groups;
    
    foreach my $file (keys %{$files}) {
	my $group=$files->{$file}->[2] || 'All';
	push @{$groups{$group}},$files->{$file}->[0];
    }

    return(\%groups);
}
