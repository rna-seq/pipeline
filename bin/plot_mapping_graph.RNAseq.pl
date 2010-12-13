#!/soft/bin/perl

use strict;
use warnings;

# Objective
# This script should build for a mapped file the following stats:
# Total number of reads, number of reads mapped, number of unique mappings
# within the mismatch threshold and number of 1:0:0 matches
# it will also use R to build a graphical summary with this info.

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file);
use Getopt::Long;

my $graph;
my $species;
my $project;
my $prefix;
my $type;

GetOptions('-graph=s' => \$graph,
	   '-type=s' => \$type);

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};

unless ($species && $project) {
    die "Species and project are unknown ???\n";
}

unless ($graph) {
    die "No input files\n";
}

build_graph($graph,
	    $species,
	    $project,
	    $type);

exit;

sub build_graph {
    my $graph=shift;
    my $species=shift;
    my $project=shift;
    my $plottype=shift;
    $graph=~s/.txt$//;
    my $lanes=`wc -l $graph.txt | gawk '{print \$1}'`;
    if ($lanes < 2) {
	print STDERR "Only one lane, barplot is not necessary\n";
	my $command=`touch $graph.ps`;
	system($command);
	return();
    }

    my $name_size=1 - (0.05 * $lanes);
    if ($name_size < 0.2) {
	$name_size=0.2;
    }
    # Create only a postscript graph
    my %figures=('ps' => "postscript(\"$graph.ps\")\n");

    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;

    # Read the data
    $r_string ="stats<-read.table(\"$graph.txt\",sep=\"\t\")\n";

    # Calculate the numbers
    $r_string.="matrix<-as.matrix(stats[,2:5])\n";
    $r_string.="stats2<-stats[,2:5] / stats[,2] * 100\n";
    $r_string.="matrix2<-as.matrix(stats2[,1:4])\n";
    $r_string.="stats3<-stats[,2:5] - c(stats[,3:5],0)\n";
    $r_string.="matrix3<-as.matrix(stats3[,1:4])\n";
    $r_string.="stats4<-stats2[,1:4] - c(stats2[,2:4],0)\n";
    $r_string.="matrix4<-as.matrix(stats4[,1:4])\n";

    # set some figure parameters
    $r_string.="cols=rainbow(5)\n";

    for my $type (keys %figures) {
	$r_string.=$figures{$type};
	$r_string.="par(oma=c(2,0,2,0))\n";
	$r_string.="layout(matrix(1:3,nrow=1,byrow=T),widths=c(3,1,3))\n";

	# Plot the barplots
	$r_string.="barplot(t(matrix),beside=T,main=\"Absolute counts\",names.arg=stats\$V6,cex.names=$name_size,col=cols[1:4])\n";
	$r_string.="par(mfg=c(1,2))\n";
	$r_string.="old.mar<-par(\"mar\")\n";
	$r_string.="par(mar=c(0,0,0,0))\n";
	$r_string.="legend('center',legend=c('TotalReads','MappedReads','UniqueReads','1:0:0 Reads','Unmapped'),fill=cols)\n";
	$r_string.="par(mar=old.mar)\n";
	$r_string.="matrix4<- matrix4[,ncol(matrix4):1]\n";
	$r_string.="barplot(t(matrix4),beside=F,main=\"Percentage distribution\",names.arg=stats\$V6,cex.names=$name_size,col=cols[c(4,3,2,5)])\n";
	$r_string.="title(main=\"$species $project $plottype mapping summary\",outer=T)\n";
	$r_string.="dev.off()\n";
    }

    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla < $execution_file";
    system($command);
    $command="rm $execution_file";
    system($command);
}    
