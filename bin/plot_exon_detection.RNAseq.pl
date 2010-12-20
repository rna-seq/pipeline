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
# This script will take as an input a list of overlap files 
# From the database the script will determine the total number of
# of features, and it will determine in each of the other input files how many
# of these are included, and what fraction
# It will also print a graph with the expression level measured by number of
# reads for the different features as weel as the average expression across all
# lanes for each of the expressed features
# If an annotation is provided the fraction off features from the annotation
# that is detected in each file will be output.

use RNAseq_pipeline3 qw(get_fh get_Average);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use Getopt::Long;

# Declare some varianbles
my $species;
my $project;
my $paired;
my $prefix;
my $annotation;
my $exondir;

my %features;
my @files;
my %expressed;

my $threshold=1;
my $dataset='Exons';
my $max=0;

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$annotation=$exondir."/$prefix.exon.gtf";

# Read the command lin options
GetOptions('threshold=i' => \$threshold);

my %distribution;

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};
my $lanes=keys %lanes;

# Get the number of detected features in each lane
foreach my $lane (keys %lanes) {
    my $type;
    if ($lanes{$lane} == 1) {
	$type='single';
    } elsif ($lanes{$lane} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my $infilename=$exondir.'/'.$lane.'.'.$type.'.unique.gtf.overlap.total';
    print STDERR "Processing $infilename\n";

    unless (-r $infilename) {
	die "$infilename is not readable\n";
    }

    build_distribution($infilename,
		       \%distribution,
		       $lane,
		       \$max,
		       \%expressed,
		       \%features);

    print STDERR "done\n";
    
}

# Get the average number of reads for each of the features
my %averages=%{get_averages(\%features,
			    $lanes)};
    


# Extract from the annotationthe number of
# features in the file and use this as the limit to the y axis
my $no_features=0;
if ($annotation) {
    ($no_features)=split(/\s+/,`cut -f 1,4,5,7 $annotation|sort |uniq| wc -l`);
} else {
    ($no_features)=sort {$b <=> $a} values %expressed;
}

# Plot expressed
plot_expressed(\%expressed,
	       $dataset,
	       $no_features,
	       $prefix);

# Plot the expression of each feature in each lane
plot_by_feature(\%distribution,
		\%averages,
		$dataset,
		$no_features,
		$max,
		$prefix);

# clean up
my $command="rm $dataset.txt";
system($command);

exit;

sub build_distribution {
    my $file=shift;
    my $dist=shift;
    my $lane=shift;
    my $max=shift;
    my $expressed=shift;
    my $by_feature=shift;

    my $infh=get_fh($file);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$reads,$length,$type)=split("\t",$line);
	if ($$max < $reads) {
	    $$max=$reads;
	}
	if ($dist->{$lane}{$feature}) {
	    warn $dist->{$lane}{$feature}," duplicated\n";
	} else {
	    $dist->{$lane}{$feature}=$reads;
	    push @{$by_feature->{$feature}},$reads;
	}
	if ($reads >= $threshold) {
	    $expressed->{$lane}++;
	}
    }

    close($infh);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}++;
    }

    return(\%lanes);
}

sub plot_by_feature {
    my $dist=shift;
    my $averages=shift;
    my $graph=shift;
    my $no_features=shift;
    my $max=shift;
    my $prefix=shift;
    $max*=1.2;

    my $tmpfn="$graph.txt";
    # Print a temporary file with the expression levels
    my $tmpfh=get_fh($tmpfn,1);
    my @lanes=sort keys %{$dist};
    foreach my $feature (keys %{$averages}) {
	my @reads;
	foreach my $lane (@lanes) {
	    my $reads=$dist->{$lane}->{$feature} || 0;
	    push @reads,$reads;
	}

	print $tmpfh join("\t",
			  $averages->{$feature},
			  @reads),"\n";
    }   
    close($tmpfh);

    # This sorts the file according to the average value of the feature
    my $command="sort -o $graph.txt -n -r -k 1 $graph.txt";
    system($command);

    # Plot the graph
    my $names='"'.join('","',@lanes).'"';
    my $lanes=@lanes;
    print STDERR $names,"\n";
    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;
    my $filename=$prefix.'.'.$graph.'.split';

    # Initialize the colors & symbols
    $r_string.="cols=topo.colors($lanes)\n";

    # Read the data
    $r_string.="features<-read.table(\"$graph.txt\",sep=\"\t\")\n";

    $r_string.="postscript(\"$filename.ps\")\n";
    $r_string.="plot(c(1:length(features[,1])),features[[1]],log=\"x\",xlab=\"$graph expression Rank\",type=\"l\",col=\"red\",pch=0,ylim=c(0,$max),ylab=\"Number of reads (log10)\")\n";
    $r_string.="for (n in 2:length(features[1,])) {lines(features[[n]],col=cols[n - 1],pch=0)}\n";
    $r_string.="lines(features[[1]],col=\"red\")\n";
    $r_string.="legend(\"topright\",c($names),col=cols,lty=1,cex=0.6,ncol=2)\n";
    $r_string.="dev.off()\n";

    $r_string.="jpeg(\"$filename.jpeg\")\n";
    $r_string.="plot(c(1:length(features[,1])),features[[1]],log=\"x\",xlab=\"$graph expression Rank\",type=\"l\",col=\"red\",pch=0,ylim=c(0,$max),ylab=\"Number of reads (log10)\")\n";
    $r_string.="for (n in 2:length(features[1,])) {lines(features[[n]],col=cols[n - 1],pch=0)}\n";
    $r_string.="lines(features[[1]],col=\"red\")\n";
    $r_string.="legend(\"topright\",c($names),col=cols,lty=1,cex=0.6,ncol=2)\n";
    $r_string.="dev.off()\n";

    # Print execution file
    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    $command="R --vanilla < $execution_file";
    system($command);
    $command="rm $execution_file";
    system($command);
}


sub plot_expressed {
    my $dist=shift;
    my $graph=shift;
    my $no_features=shift;
    my $prefix=shift;

    my $tmpfn="$graph.txt";
    # Print a temporary file with the expression levels
    my $tmpfh=get_fh($tmpfn,1);
    my @lanes=sort keys %{$dist};
    foreach my $lane (@lanes) {
	print $tmpfh join("\t",
			  $lane,
			  $dist->{$lane}),"\n";
    }

    # Close the file to make sure buffer is flushed
    close($tmpfh);

    # Build the R command file
    my $names='"'.join('","',@lanes).'"';
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;
    my $lanes=@lanes;
    $r_string ="cols=topo.colors($lanes)\n";

    $r_string.="stats<-read.table(\"$graph.txt\",sep=\"\t\")\n";
    $r_string.="postscript(\"$prefix.$graph.detection.ps\")\n";
    $r_string.='barplot(stats$V2,beside=T,axes=T,';
    $r_string.="main=\"$graph detected from the $no_features in the annotation\",";
    $r_string.="xlab=\"Lanes\",ylab=\"$graph detected\",";
    $r_string.="col=cols,ylim=c(0,$no_features))\n";
    $r_string.="legend(\"topright\",c($names),col=cols,";
    $r_string.="cex=0.6,ncol=2,fill=cols)\n";
    $r_string.="dev.off()\n";

    $r_string.="jpeg(\"$prefix.$graph.detection.jpeg\")\n";
    $r_string.='barplot(stats$V2,beside=T,axes=T,';
    $r_string.="main=\"$graph detected from the $no_features in the annotation\",";
    $r_string.="xlab=\"Lanes\",ylab=\"$graph detected\",";
    $r_string.="col=cols,ylim=c(0,$no_features))\n";
    $r_string.="legend(\"topright\",c($names),col=cols,";
    $r_string.="cex=0.6,ncol=2,fill=cols)\n";
    $r_string.="dev.off()\n";

    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet < $execution_file";
    system($command);
    $command="rm $execution_file";
    system($command);
}

sub get_averages {
    my $features=shift;
    my $lanes=shift;

    my %averages;

    # First initialize any uninitialized value in the hash
    foreach my $feature (keys %{$features}) {
	for (my $i=0;$i<$lanes;$i++){
	    unless(defined $features->{$feature}->[$i]) {
		$features->{$feature}->[$i]=0;
	    }
	}
    }
    
    # After get the average for each
    foreach my $feature (keys %{$features}) {
	$averages{$feature}=get_Average($features{$feature});
    }
    
    return(\%averages);
}
