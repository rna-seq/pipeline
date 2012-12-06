#!/soft/bin/perl
# DGK & A.Merkel

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

#    Author : David Gonzalez, Angelika Merkel david.gonzalez@crg.eu

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
my $bindir='';

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
$bindir=$options{'BINDIR'};

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
	print $log_fh "I have two replicates for $sample\n";
	print $log_fh "printing IDR to $outputfile\n";
#	my $command="/users/rg/amerkel/bin/RAPexonRPKM_to_IDR_matchedPeaks.pl ";
        my $command=$bindir."/RAPexonRPKM_to_IDR_matchedPeaks.pl";
	$command.="-i1 $files[0] ";
	$command.="-i2 $files[1] ";
	$command.="-o $outputfile";
	run_system_command($command);
	$command="rm $files[0] $files[1]";
	run_system_command($command);
	
	# Run the IDR code
#	$command="/users/rg/amerkel/bin/runR.sh /users/rg/amerkel/projects/ENCODE/gene_expression/IDR-discret.Qunha/batch-IDR-mappedInput-discrete-p.r $outputfile $sample T 10 0.95";
	$command=$bindir."/runR.sh ".$bindir."/batch-IDR-mappedInput-discrete-p.r $outputfile $sample T 10 0.95";
	run_system_command($command);
	
	# IDR code fix: reassign the element IDs back to the peaks (they are not
	# kept in the original code)
#	$command="/users/rg/amerkel/bin/Idrdiscrete_to_overlapPeaks.pl--infile *-categories.txt --ref *matchedPeaks.txt --outfile somemeaningfulname";
	$command=$bindir."/Idrdiscrete_to_overlapPeaks.pl--infile *-categories.txt --ref *matchedPeaks.txt --outfile idrout.txt";
	run_system_command($command);

	last;
    }
}

close($log_fh);

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

