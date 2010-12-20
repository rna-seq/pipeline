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

### TO DO
# Parallelize on the cluster
# Make the mapper take a list of files

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_log_fh);
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
my $log_fh=get_log_fh('run_flux.RNAseq.log',
		      $debug);

# First get a list of the bed files we are going to process
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# get the necesary directories we are going to use in order to build the command
# line for the capacitor

# split the annotation
my %annot;
if (-r $exonfile) {
    print STDERR "Processing $exonfile\n";
    my $exonfh=get_fh($exonfile);
    split_by_chrom($localdir,
		   'annot',
		   $exonfh,
		   \%annot,
		   $log_fh);
    close($exonfh);
    # close the file handles;
    foreach my $file (keys %annot) {
	close($annot{$file}->[0]);
    }
} else {
    die "Can't read $exonfile\n";
}

# Run the capacitor for each of the files that we have and move the output to
# a file with the adequate name
foreach my $lane (keys %lanes) {
    my $type;
    if (@{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (@{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }
    my %reads;

    # Check if the file corresponding to this lane exists already
    my $fluxfile=$lane.'.flux.gtf';
    if (-r $fluxfile) {
	print STDERR $fluxfile,"\tPresent. Skipping...\n";
	next;
    }

    # Get the files corresponding to the both halves of the reads and
    # split them by chromosome
    my $infilename=$localdir.'/'.$lane.'.combined.bed';
    if (-r $infilename) {
	print $log_fh "Processing $infilename\n";
    } else {
	die "Can't read $infilename\n";
    }
    my $infh=get_fh($infilename);
    split_by_chrom($localdir,
		   'reads',
		   $infh,
		   \%reads,
		   $log_fh);
    close($infh);
    # close the file handles;
    foreach my $file (keys %reads) {
	close($reads{$file}->[0]);
    }

    my %results;
    # Run the flux capacitor
    # This is the parallel step
    foreach my $chr (keys %annot) {
	unless ($reads{$chr}) {
	    print STDERR "No reads mapped to $chr\n";
	    next;
	}

	my $outfile=run_flux($localdir,
			     $type,
			     $reads{$chr}->[1],
			     $annot{$chr}->[1],
			     $bindir,
			     $tmpdir,
			     $sorted);
	
	if ($outfile) {
	    $results{$chr}=$outfile;
	} else {
	    die "Problem generating $outfile\n";
	}

	if ($debug) {
	    last;
	}
    }

    # merge the flux results
    my @files;
    foreach my $chr (keys %results) {
	if (-r $results{$chr}) {
	    push @files, $results{$chr};
	}
    }
    my $command='cat ';
    $command.=join(' ',@files);
    $command.=" > $fluxfile";
    print $log_fh "Executing:\t$command\n";
    system($command);

    # clean up files
    my @clean;
    foreach my $chr (keys %annot) {
	if ($results{$chr} && -r $results{$chr}) {
	    push @clean, $results{$chr};
	}
    }
    $command='rm ';
    $command.=join(' ',@clean);
    print $log_fh "Executing:\t$command\n";
    system($command);
}

# clean up the sorted annotation files
foreach my $chr (keys %annot) {
    # delete the files
    my $file=$annot{$chr}->[1];
    my $command="rm $file";
    print $log_fh "Executing $command\n";
    system($command);
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

### TO DO put into a module
sub split_by_chrom {
    my $localdir=shift;
    my $type=shift;
    my $fh=shift;
    my $hash=shift;
    my $log_fh=shift;

    unless($localdir) {
	$localdir='/tmp';
    }

    while (my $line=<$fh>) {
	chomp($line);
	my @line=split("\t",$line);
	my $chrom=$line[0];
	if ($chrom=~/^Un/) {
	    next;
	}
	unless (defined $hash->{$chrom}) {
	    my $tmpfile=$localdir.'/'.$chrom.'.'.$$.$type;
	    $hash->{$chrom}->[0]=get_fh($tmpfile,1);
	    $hash->{$chrom}->[1]=$tmpfile;
	    print $log_fh "Creating: ",$tmpfile,"\n";
	}

	my $outfh=$hash->{$chrom}->[0];
	# Fix the read Id if paired bed
	if ($line[3]=~s/\|p/\//) {
	    $line=join("\t",@line);
	}
	print $outfh $line,"\n"; 
    }
}
