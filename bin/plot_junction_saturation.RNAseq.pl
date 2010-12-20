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
			       'get_saturation_curve','get_junction_no');
use Getopt::Long;

# Declare some variables
my $threshold=1;
my $species;
my $project;
my $prefix;
my $junctiondir;
my $junction_number;
my $filetype='Junction';

GetOptions('threshold|t=i' => \$threshold);

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$junctiondir=$options{'JUNCTIONSDIR'};
$junction_number=get_junction_no();

print STDERR $junction_number,"\tJunctions present in the combinations\n";

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
    my $filename=$junctiondir.'/'.$lane.'.'.$type.'.unique.overlap.total';

    if (-r $filename) {
	print STDERR "Processing $filename\n";
    } else {
	die "Can't read $filename\n";
    }

    my $infh=get_fh($filename);
    while (my $line=<$infh>) {
	chomp($line);
	my @line=split(/\t/,$line);
	my ($gene_id,$overlaps)=@line[0,1];

	if ($overlaps >= $threshold) {
	    $detected{$filename}{$gene_id}=1;
	}
    }
    close($infh);
}

# Build the curve
my @curve=@{get_saturation_curve(\%detected)};

# Plot the results
plot_saturation_point(\@curve,
		      $junction_number,
		      $filetype);

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[1]}++;
    }

    return(\%lanes);
}

sub plot_saturation_point {
    my $curve=shift;
    my $total=shift;
    my $filetype=shift;

    my $filename=$filetype.'.saturation.ps';

    # print a temporary file with the data
    my $tmpfile="$$.saturation.points.txt";
    my $outfh=get_fh($tmpfile,1);
    foreach my $point (@{$curve}) {
	my $fraction=sprintf "%.2f",$point->[1] / $total;
	$point->[0]=~s/.*\///;
	$point->[0]=~s/\.(paired|single)\.unique.overlap.total//;
	print $outfh join("\t",
			  @{$point},
			  $fraction,
			  $threshold),"\n";
    }
    close($outfh);

    # Calculate the positions of the second axis
    my @second_axis=(0,
		     $total/4,
		     $total/2,
		     (3 * $total)/4,
		     $total);
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
    $r_string ="cols=rainbow($lanes)\n";
    $r_string.="stats<-read.table(\"$tmpfile\",sep=\"\t\")\n";
    $r_string.="postscript(\"$filename\")\n";
    $r_string.='plot(stats$V2,names.arg=c(stats$V1),type="l",';
    $r_string.="main=\"Saturation\",";
    $r_string.="xlab=\"Lanes\",ylab=\"Features detected\",legend=T,";
    $r_string.="ylim=c(1,$total),col=cols,log=\"y\")\n";
    $r_string.="axis(4,at=$axis2,labels=c(0,0.25,0.5,0.75,1))\n";
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
