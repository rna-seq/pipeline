package RNAseq_pipeline_comp3;

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('get_tables','check_tables','get_labels_sub','get_samples');

use strict;
use warnings;

# Load other modules required
use POSIX qw(uname);
use RNAseq_pipeline3 qw(check_table_existence);

# Get the tables that we will use from the expreiments table
sub get_tables {
    my $dbh=shift;
    my $project=shift;
    my $table_sufix=shift;
    my $fraction=shift;

    my %tables;

    my ($query,$sth,$count);
    $query ='SELECT experiment_id ';
    $query.='FROM experiments ';
    $query.='WHERE project_id = ?';
    if ($fraction) {
	$query.=" $fraction";
    }
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

sub check_tables {
    my $database=shift;
    my $tables=shift;
    
    my %remove;
    foreach my $exp (keys %{$tables}) {
	my ($query,$sth);
	my $table=$tables->{$exp};

	my $present=check_table_existence($database,
					  $table);
	
	unless($present) {
	    # Remove the table from the dataset
	    $remove{$exp}=1;
	}
    }

    # Remove the missing tables from the analysis
    foreach my $exp (keys %remove) {
	delete $tables->{$exp};
    }
}

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
	    $samples{$sample_id}=[$table,$exp];
	}	
    }
    return(\%samples);
}

1;
