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
# This script will take as an input the overlap.total result of an experiment
# and it will extract from it the features that have at least one (or the 
# threshold number) of completely included reads
# It will use this to build a saturation curve for the lanes

# Load some modules
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_unique_transcripts',
			       'get_saturation_curve','get_dbh');
use Getopt::Long;

# Declare some variables
my $threshold=1;
my $species;
my $project;
my $prefix;
my $transdir;
my $exon_number;
my $filetype='Transcript';

GetOptions('threshold|t=i' => \$threshold);

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$transdir=$options{'TRANSDIR'};
$exon_number=get_unique_transcripts();

print STDERR $exon_number,"\tUnique transcripts present in the annotation\n";

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};
my $lanes=keys %lanes;

print STDERR $lanes,"\tLanes present\n";

my %detected;

# Read and process the files
foreach my $lane (keys %lanes) {
    my $type;
    if ($lanes{$lane} == 1) {
	$type='single';
    } elsif ($lanes{$lane} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }
    my %genes;
#    my $filename=$lane.'.'.$type.'.unique.gtf.gz';
    my $filename=$transdir.'/'.$lane.'.'.$type.'.unique.gtf.gz';

    if (-r $filename) {
	print STDERR "Processing $filename\n";
    } else {
	die "Can't read $filename\n";
    }

    my $infh=get_fh($filename);
    while (my $line=<$infh>) {
	chomp($line);
	my @line=split(/\t/,$line);
	my ($trans_id)=$line[0];

	$detected{$filename}{$trans_id}=1;
    }
    close($infh);
}

# Build the curve
my @curve=@{get_saturation_curve(\%detected)};

# Plot the results
plot_saturation_point(\@curve,
		      $exon_number,
		      $filetype);

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}++;
    }

    return(\%lanes);
}

sub plot_saturation_point {
    my $curve=shift;
    my $total=shift;
    my $filetype=shift;

    my $filename=$filetype.'.saturation';

    # print a temporary file with the data
    my $tmpfile="$$.saturation.points.txt";
    my $outfh=get_fh($tmpfile,1);
    my $max=0;
    foreach my $point (@{$curve}) {
	my $fraction=sprintf "%.2f",$point->[1] / $total;
	$point->[0]=~s/.*\///;
	$point->[0]=~s/\.(paired|single)\.unique.gtf.overlap.total//;
	print $outfh join("\t",
			  @{$point},
			  $fraction,
			  $threshold),"\n";
	if ($point->[1] > $max) {
	    $max=$point->[1];
	}
    }
    close($outfh);

    # Set the ylimit
    $max+=1000;

    # Calculate the positions of the second axis
    my @second_axis=(0,
		     $total/4,
		     $total/2,
		     (3 * $total)/4,
		     $total);

    my @second_axis_labels=(0,0.25,0.5,0.75,1);
    my $points=0;
    foreach my $point (@second_axis) {
	if ($point <= $max) {
	    $points++;
	}
    }

    @second_axis=splice(@second_axis,0,$points);
    @second_axis_labels=splice(@second_axis_labels,0,$points);

    # build an x axis with the lanes
    my @lanes;
    foreach my $entry (@{$curve}) {
	push @lanes,$entry->[0];
    }

    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;
    my $lanes=@{$curve};
    my $lanenames='c("'.join('","',@lanes).'")';
    my $axis2='c('.join(',',@second_axis).')';
    my $axis2labels='c('.join(',',@second_axis_labels).')';
    $r_string ="cols=rainbow($lanes)\n";
    $r_string.="stats<-read.table(\"$tmpfile\",sep=\"\t\")\n";

    $r_string.="postscript(\"$filename.ps\")\n";
    $r_string.='plot(stats$V2,type="l",';
    $r_string.="main=\"Saturation\",";
    $r_string.="xlab=\"Lanes\",ylab=\"Features detected\",";
    $r_string.="ylim=c(0,$max),col=cols)\n";
    $r_string.="axis(4,at=$axis2,labels=$axis2labels)\n";
    $r_string.="axis(3,at=c(1:$lanes),labels=$lanenames)\n";
    $r_string.="dev.off()\n";

    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet < $execution_file";
    system($command);
    $command="rm $tmpfile";
    system($command);
    $command="rm $execution_file";
    system($command);
}
