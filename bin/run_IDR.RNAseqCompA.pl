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
# Run the IDR analysis on the exons

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command check_table_existence);
use RNAseq_pipeline_settings3 qw(get_dbh read_config_file);
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $breakdown;
my $tabsuffix='gene_RPKM_pooled';
my $fraction='';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'limit|l=s' => \$fraction,
	   'breakdown|b' => \$breakdown);

if ($breakdown) {
    $tabsuffix='gene_RPKM';
}

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_gene_RPKM.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get pairs of samples to be used
my %pairs=%{get_sample_pairs($dbhcommon,
			     $project)};

foreach my $sample (keys %pairs) {
    print $sample,"\n";
    
    my @files;
    foreach my $bio (keys %{$pairs{$sample}}) {
	my $table_id=$pairs{$sample}{$bio};
	$table_id.='_gene_RPKM';
	my $table_ok=check_table_existence($dbh,
					   $table_id);

	if ($table_ok) {
	    my $file_id=$table_id.'.txt.gz';
	    write_gene_table($dbh,
			     $table_id,
			     $file_id);
	    push @files, $file_id;
	}
    }
    if (@files==2) {
	my $outputfile=$sample.'.IDR.out';
	print STDERR "I have two replicates for $sample\n";
	print STDERR "printing IDR to $outputfile\n";
	my $command="/users/rg/amerkel/bin/RAPexonRPKM_to_IDR_matchedPeaks.pl ";
	$command.="-i1 $files[0] ";
	$command.="-i2 $files[1] ";
	$command.="-o $outputfile";
	run_system_command($command);
	$command="rm $files[0] $files[1]";
	run_system_command($command);
	
	#run the IDR code
	$command="/users/rg/amerkel/bin/runR.sh /users/rg/amerkel/projects/ENCODE/gene_expression/IDR-discret.Qunha/batch-IDR-mappedInput-discrete-pRAP.r /users/rg/amerkel/projects/ENCODE/gene_expression/IDR-discret.Qunha/functions-09-12-2010.r $outputfile $sample T 10 0.95";
	run_system_command($command);
	
	#IDR code fix: reassign the element IDs back to the peaks (they are not kep in the original code)
	$command="/users/rg/amerkel/bin/Idrdiscrete_to_overlapPeaks.pl--infile *-categories.txt --ref *matchedPeaks.txt --outfile somemeaningfulname";
	run_system_command($command);

	last;
    }
}

# Remove any tables that do not exist
#check_tables($dbh,
#	     \%tables);

# For each of tables extract the RPKMs of interest and get for each of the
# tables the different samples present in them
#my %samples=%{get_samples(\%tables,
#			  $dbh,
#			  $breakdown)};
#my @experiments=sort {$a cmp $b} keys %samples;
#my @values;
#my %all_genes;
#foreach my $experiment (@experiments) {
#    my ($table,$sample)=split('_sample_',$experiment);
#    print $log_fh "Extracting $sample, data from $table\n";
#    my $data=get_RPKM_data($dbh,
#			   $table,
#			   \%all_genes,
#			   $sample,
#			   $breakdown);
#    if ($data) {
#	push @values, [$experiment,$data];
#    } else {
#	print STDERR "Skipping $experiment\n";
#    }
#}

# Get the human readable lables
#oreach my $experiment (@experiments) {
#   my $label;
#   if ($nolabels) {
#	$label=$samples{$experiment}->[1];
#   } else {
#	$label=get_labels($experiment);
#   }
#   if ($label) {
#	$experiment=$label;
#   }
#

# Print the expression values for each gene in each of the tables into a
# temporary file
my $tmpfn="Expression.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
#print $tmpfh join("\t",@experiments),"\n";
#foreach my $gene (keys %all_genes) {
#    my @row;
#    my $no_print=0;
#    foreach my $exp (@values) {
#	my $value=0;
#	if ($exp->[1] &&
#	    ($exp->[1]->{$gene})) {
#	    $value=$exp->[1]->{$gene};
#	} else {
#	    $no_print=1;
#	}
#	push @row,$value;
#    }
#    unless ($no_print) {
#	print $tmpfh join("\t",
#			  $gene,
#			  @row),"\n";
#    }
#}
#close($tmpfh);

# Plot the correlation graph if we have more than 2 samples
#if (@experiments > 2) {
#    my $outfn=$project.".clusters";
#    plot_graphs_R($tmpfn,
#		  $outfn);
#} else {
#    print STDERR "Not enough samples to cluster\n";
#}

exit;

sub get_sample_pairs {
    my $dbh=shift;
    my $project=shift;

    my $table='experiments';
    my ($sth,$query,$count);
    $query ='SELECT experiment_id,CellType,RNAType,Compartment,';
    $query.='Bioreplicate ';
    $query.="FROM $table ";
    $query.='WHERE project_id = ?';
    $query.='AND CellType is not NULL ';
    $query.='AND RNAType is not NULL ';
    $query.='AND Compartment is not NULL ';
    $query.='AND Bioreplicate is not NULL';
    $sth=$dbh->prepare($query);
    $sth->execute($project);

    my %pairs;
    my %remove;
    while (my ($exp,$cell,$RNA,$comp,$bio)=$sth->fetchrow_array()) {
	my $id=join('_',$cell,$RNA,$comp);
	$pairs{$id}{$bio}=$project.'_'.$exp;
    }
    $count=keys %pairs;
    print STDERR $count,"\tSamples are in the database for $project\n";

    foreach my $sample (keys %pairs) {
	my $replicates=keys %{$pairs{$sample}};
	unless ($replicates == 2) {
	    $remove{$sample}=1;
	}
    }

    foreach my $sample (keys %remove) {
	delete $pairs{$sample};
    }
    $count=keys %pairs;
    print STDERR $count,"\t$project samples have exactly 2 bioreplicates\n";

    return(\%pairs);
}

sub write_gene_table {
    my $dbh=shift;
    my $table=shift;
    my $file=shift;

    my $outfh=get_fh($file,1);

    my ($query,$sth,$count);
    $query ='SELECT gene_id, RPKM ';
    $query.="FROM $table ";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();
    
    if ($count && ($count > 1)) {
	print STDERR $count,"\tGenes are detected in $table\n";
    } else {
	die "No genes present in $table\n";
    }

    # get all the necessary tables
    while (my ($gene,$rpkm)=$sth->fetchrow_array()) {
	print $outfh join("\t",
			  $gene,$rpkm),"\n";
    }
    close($outfh);
}

# Build a postscript tree using the  canberra distance and the complete linkage
# for clustering
sub plot_graphs_R {
    my $statsfn=shift;
    my $outfile=shift;

    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;

    # Read the data
    $r_string.="rpkms<-read.table(\"$statsfn\",sep=\"\t\",header=T)\n";

    # Calculate the distance matrix
    $r_string.="genesdist<-dist(t(rpkms),method='canberra')\n";

    # Setup the figure
    $r_string.="pdf(\"$outfile.pdf\")\n";

    # Build the tree
    $r_string.="plot(hclust(genesdist))\n";
    $r_string.="dev.off()\n";
	
    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet < $execution_file";
    run_system_command($command);

    $command="rm $execution_file";
    run_system_command($command);
}

