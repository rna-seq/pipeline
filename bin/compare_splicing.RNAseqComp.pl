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
# This script should take a project from the database and for each of the
# experiments belonging to it build a comparison table that can be run
# through R This table will contain for each junction the number of reads
# detected.
# We have to be able to decide what table to get the data from (pooled or not)
# And also we have to be able to select genes that are expressed 

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 qw(get_dbh read_config_file);
use RNAseq_pipeline_stats3 qw(log10);
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_samples');
use Getopt::Long;

# Declare som variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $breakdown=0;
my $tabsuffix='all_junctions_class_pooled';
my $limit;

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'breakdown|b' => \$breakdown,
	   'limit|l=s' => \$limit);

if ($breakdown) {
    $tabsuffix='all_junctions_class';
}

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_splicing.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);
*junc2gene=get_gene_from_short_junc_sub($dbhcommon);

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			$limit)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# For each of tables extract the splice site support of interest and get for
# each of the tables the different samples present in them
my %samples=%{get_samples(\%tables,
			  $dbh,
			  $breakdown)};
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_splicing_data($dbh,
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

# Get the human readable labels
foreach my $experiment (@experiments) {
    my $label;
    if ($nolabels) {
	$label=$samples{$experiment}->[1];
    } else {
	$label=get_labels($experiment);
    }
    if ($label) {
	$experiment=$label;
    }
}

# Print the detected reads for each junction in each of the tables into a
# temporary file
my $tmpfn="Junctions.$project.txt";
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
	my $gene_id=join('_',@{junc2gene($gene)});
	print $tmpfh join("\t",
			  $gene_id.'_'.$gene,
			  @row),"\n";
    }
}
close($tmpfh);

exit;

sub get_splicing_data {
    my $dbh=shift;
    my $table=shift;
    my $all=shift;
    my $sample=shift;
    my $breakdown=shift;

    my %expression;

    my ($query,$sth,$count);
    $query ='SELECT chr1, start, chr2, end, support ';
    $query.="FROM $table ";
    if ($breakdown) {
	$query.='WHERE LaneName = ?';
    } else {
	$query.='WHERE sample = ?';
    }
    $query.=' AND junc_type not like "split%"';
    $sth=$dbh->prepare($query);
    $count=$sth->execute($sample);
    
    if ($count && ($count > 1)) {
	print STDERR $count,"\tGenes are detected in $table\n";
    } else {
	die "No genes present in $table\n";
    }

    # get all the necessary tables
    while (my ($chr1,$start,$chr2,$end,$support)=$sth->fetchrow_array()) {
	my $splice_id=join('_',
			   $chr1,$start,$end);
	if ($chr1 ne $chr2) {
	    $splice_id=join('_',
			    $chr1,$start,$chr2,$end);
	}
	$expression{$splice_id}=$support;
	$all->{$splice_id}=1;
    }

    return(\%expression);
}


# This sub should take a short junction_id and extract the gene it belongs to
# exon junction belongs
sub get_gene_from_short_junc_sub {
    my %options=%{read_config_file()};
    my $dbh=shift;
    my $table=$options{'JUNCTIONSTABLE'};

    # For saving time, as the junctions table is huge
    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT gene_id ';
    $query.="FROM $table ";
    $query.='WHERE chr = ? AND start = ? AND end = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $junc=shift;
	### TO DO fix for problematic chromosomes
	my ($chr,$start,$end)=split('_',$junc);

	unless ($cache{$junc}) {
	    $count=$sth->execute($chr,$start,$end);

	    if ($count == 0) {
		die "No gene in $table corresponds to $junc\n";
	    } else {
		while (my ($gene)=$sth->fetchrow_array()) {
		    push @{$cache{$junc}}, $gene;
		}
	    }
	}
	return($cache{$junc});
    };
    return($subroutine);
}
