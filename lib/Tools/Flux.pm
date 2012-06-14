package Tools::Flux;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('build_parameter_file','run_flux');

use strict;
use warnings;

use RNAseq_pipeline3 qw(get_fh run_system_command);

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
    print $paramfh join(" ",
			'SORT_IN_RAM',
			'false'),"\n";

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

    print STDERR $fileroot1,"\n";
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
    run_system_command($command);

    my $parameterfile=build_parameter_file($annot,
					   $maps,
					   $type,
					   $sorted_maps,
					   $outfile,
					   $tmpdir);

    # build the flux command
    $command="$bindir/flux -t capacitor -p $parameterfile -u";
    run_system_command($command);

    if (-r $outfile) {	
	print STDERR $outfile,"\tBuilt\n";
	# Clean up by removing the parameter file
	# The read file seems to be cleared by the flux
	$command="rm $parameterfile";
	run_system_command($command);
    } else {
	die "Problem building $outfile\n";
    }

    return($outfile);
}

1;
