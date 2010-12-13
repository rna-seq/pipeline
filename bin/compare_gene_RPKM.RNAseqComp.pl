#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take a project from the database and for each of the
# experiments belonging to it build a comparison table that can be run
# through R
# We have to be able to decide what table to get the data from (pooled or not)
# And also we have to be able to select genes that are expressed 

use RNAseq_pipeline3 qw(get_fh get_log_fh);
use RNAseq_pipeline_settings3 qw(get_dbh read_config_file);
use DGK_stats qw(log10);
use Getopt::Long;

# Declare som variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $breakdown;
my $tabsuffix='gene_RPKM_pooled';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
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

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables,
	     $log_fh);


# For each of tables extract the RPKMs of interest and get for each of the
# tables the different samples present in them
my %samples=%{get_samples(\%tables,
			  $dbh,
			  $breakdown)};
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_experiment_data($dbh,
				 $table,
				 \%all_genes,
				 $sample,
				 $breakdown);
    if ($data) {
	push @values, [$experiment,$data];
    } else {
	print STDERR "Skipping $experiment\n";
    }
}

# Get the human readable lables
unless ($nolabels) {
    foreach my $experiment (@experiments) {
	my $label=get_labels($experiment);
	if ($label) {
	    $experiment=$label;
	}
    }
}

# Print the expression values for each gene in each of the tables into a
# temporary file
my $tmpfn="Expression.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",@experiments),"\n";
foreach my $gene (keys %all_genes) {
    my @row;
    my $no_print=0;
    foreach my $exp (@values) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	} else {
#	    $no_print=1;
	}
	push @row,$value;
    }
    unless ($no_print) {
	print $tmpfh join("\t",
			  $gene,
			  @row),"\n";
    }
}
close($tmpfh);

# Plot the correlation graph if we have more than 2 samples
if (@experiments > 2) {
    my $outfn=$project.".clusters";
    plot_graphs_R($tmpfn,
		  $outfn);
} else {
    print STDERR "Not enough samples to cluster\n";
}

exit;

# This sub should take the information from the experiments table and where
# available build more human readable lables for the each of the datasets
sub get_labels_sub {
    my $dbh=shift;
    my $table='experiments';

    my %cache;
    my ($query,$sth,$count);

    $query ='SELECT CellType, RNAType, Compartment, Bioreplicate ';
    $query.="FROM $table ";
    $query.='WHERE project_id = ? AND experiment_id = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $table_id=shift;
	my ($proj_id,$exp_id)=split('_',$table_id);

	my $key=join('_',$proj_id,$exp_id);

	unless (exists $cache{$key} &&
		$cache{$key}) {
	    $count=$sth->execute($proj_id,$exp_id);
	    if ($count == 1) {
		my ($cell,$rna,$comp,$bio)=$sth->fetchrow_array();
		my @label;
		if ($cell) {
		    push @label, $cell;
		    if ($rna) {push @label, $rna};
		    if ($comp) {push @label, $comp};
		    if ($bio) {push @label, $bio};
		} else {
		    push @label,$exp_id;
		}
		$cache{$key}=join('.',@label);
	    } else {
		die "Wrong number of elements ($count) returned for $key\n";
	    }
	}
	return($cache{$key});
    };

    return($subroutine);
}

sub get_samples {
    my $tables=shift;
    my $dbh=shift;
    my $breakdown=shift;

    my %samples;

    foreach my $exp (keys %{$tables}) {
	my ($query,$sth);
	my $table=$tables->{$exp};
	
	if ($breakdown) {
	    $query ='SELECT distinct LaneName ';
	} else {
	    $query ='SELECT distinct sample ';
	}
	$query.="FROM $table";
	
	$sth = $dbh->prepare($query);
	my $count=$sth->execute();
	while (my ($sample)=$sth->fetchrow_array()) {
	    my $sample_id=join('_sample_',
			       $table,
			       $sample);
	    $samples{$sample_id}=$table;
	}	
    }
    return(\%samples);
}

sub check_tables {
    my $database=shift;
    my $tables=shift;
    my $log_fh=shift;
    
    print $log_fh "Checking for the presence of required tables...\n";

    my %remove;
    foreach my $exp (keys %{$tables}) {
	my ($query,$sth);
	my $table=$tables->{$exp};
	
	$query ='SELECT count(*) ';
	$query.="FROM $table";
	
	$sth = $dbh->table_info(undef,undef,$table,"TABLE");
	
	my $count=$sth->execute();
	my $results=$sth->fetchall_arrayref();
	
	if (@{$results}) {
	    # Print the table location for the junctions of this experiment
	    print $log_fh join("\t",
			       $table,
			       "Present"),"\n";
	} else {
	    # Continue, as the table must be created
	    print $log_fh $table,"\tIs not present\n";
	    $remove{$exp}=1;
	}
    }

    # Remove the missing tables from the analysis
    foreach my $exp (keys %remove) {
	delete $tables->{$exp};
    }
}

sub get_experiment_data {
    my $dbh=shift;
    my $table=shift;
    my $all=shift;
    my $sample=shift;
    my $breakdown=shift;

    my %expression;

    my ($query,$sth,$count);
    $query ='SELECT gene_id, RPKM ';
    $query.="FROM $table ";
    if ($breakdown) {
	$query.='WHERE LaneName = ?';
    } else {
	$query.='WHERE sample = ?';
    }
    $sth=$dbh->prepare($query);
    $count=$sth->execute($sample);
    
    if ($count && ($count > 1)) {
	print STDERR $count,"\tGenes are detected in $table\n";
    } else {
	die "No genes present in $table\n";
    }

    if ($count < 16000) {
	print STDERR "Too few genes detected for $sample\n";
#	return();
    }
    
    # get all the necessary tables
    while (my ($gene,$rpkm)=$sth->fetchrow_array()) {
	$expression{$gene}=$rpkm;
	$all->{$gene}=1;
    }

    return(\%expression);
}

sub get_tables {
    my $dbh=shift;
    my $project=shift;
    my $table_sufix=shift;

    my %tables;

    my ($query,$sth,$count);
    $query ='SELECT experiment_id ';
    $query.='FROM experiments ';
    $query.='WHERE project_id = ?';
#    $query.=' AND experiment_id like "Wold%" and experiment_id NOT like "%geneid%"';
#    $query.=' AND experiment_id like "Ging%" and experiment_id != "Ging003WC_spikeIN"';
    $query.=' AND experiment_id like "Short%"';
#    $query.=' limit 5';
    $sth=$dbh->prepare($query);

    $count=$sth->execute($project);
    if ($count && ($count > 1)) {
	print STDERR $count,"\tExperiments are present for $project\n";
    } else {
	die "No experiments present for project\n";
    }

    # get all the necessary tables
    while (my ($experiment)=$sth->fetchrow_array()) {
	my $table_id=join('_',
			  $project,
			  $experiment,
			  $table_sufix);
	$tables{$experiment}=$table_id;
    }
    return(\%tables);
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
    $r_string.="postscript(\"$outfile.ps\")\n";

    # Build the tree
    $r_string.="plot(hclust(genesdist))\n";
    $r_string.="dev.off()\n";
	
    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet < $execution_file";
    system($command);

    $command="rm $execution_file";
    system($command);

    # Clean up
#    $command="rm $statsfn";
#    system($command);
}

