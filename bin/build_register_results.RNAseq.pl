#!/soft/bin/perl

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

#    Author : David Gonzalez, david.gonzalez@crg.eu

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
# This script will store all the results that need to be kept in the results
# folder and gzip any file that is not gzipped

# Load some modules
use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command get_md5sum check_table_existence);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh',
			       'get_exp_info_sub');
use Bio::SeqIO;

# Get some options from the configuration file
my $project;
my $tmpdir;
my $samdir;
my $file_list;
my $debug=0;
my $proj_id;
my $exp_id;
my $prefix;

my %options=%{read_config_file()};
$project=$options{'PROJECT'};
$tmpdir=$options{'LOCALDIR'};
$samdir=$options{'PROJECT'}.'/SAM';
$file_list=$options{'FILELIST'};
$prefix=$options{'PREFIX'};
$proj_id=$options{'PROJECTID'};
$exp_id=$options{'EXPID'};

# Get some subs
*get_exp_info=get_exp_info_sub('md5sum');

# Get the files we are going to process
my %files=%{read_file_list()};

# Get a log file
my $log_fh=get_log_fh('build_register_results.RNAseq.log',
		      $debug);

# Connect to the database
my $dbh=get_dbh();

my @files;

# Get the merged mapping files
get_merged_mapping(\@files,
		   $prefix,
		   $project,
		   $dbh);

# Get the merged bam files
get_merged_bam(\@files,
	       $prefix,
	       $project,
	       $dbh);


# Get the gtf files 
get_gtf_files(\@files,
	      \%files,
	      $project);

# Get the bed files
get_bed_files(\@files,
	      \%files,
	      $project,
	      $tmpdir);

# Print out the results
my ($globalmd5)=@{get_exp_info($proj_id,
			       $exp_id)};

my $resultsdir=$project.'/results';
foreach my $file (@files) {
    $file->[1]=~s/\//./;
    $file->[1]=~s/^work\.//;
    my $source=$file->[3];
    my $target=$resultsdir.'/'.$file->[1];
    if (-r $target) {
	print STDERR $target,"\tpresent. Skipping...\n";
    } else {
	my $command="cp $source $target";
	run_system_command($command);
    }
    $file->[3]=$target;
    print join("\t",
	       @{$file},
	       $globalmd5),"\n";
}

exit;

sub get_merged_mapping {
    my $files=shift;
    my $prefix=shift;
    my $project=shift;
    my $dbh=shift;

    my $table=$prefix.'_merged_mapping';
    
    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT filename ';
	$query.="FROM $table ";
	$sth=$dbh->prepare($query);
	$sth->execute();

	while (my ($filename)=$sth->fetchrow_array()) {
	    my $filepath=$project.'/SAM/'.$filename;
	    my ($file,$md5sum)=get_md5sum($filepath);
	    push @{$files},['gem.mapping',$filename,$md5sum,$filepath];
	}
    }
}

sub get_merged_bam {
    my $files=shift;
    my $prefix=shift;
    my $project=shift;
    my $dbh=shift;

    my $table=$prefix.'_merged_SAM';
    
    my ($query,$sth);
    $query ='SELECT BAM_file ';
    $query.="FROM $table ";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($filename)=$sth->fetchrow_array()) {
	my $filepath=$project.'/SAM/'.$filename;
	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['bam',$filename,$md5sum,$filepath];
    }
}

sub get_gtf_files {
    my $files=shift;
    my $fileprefix=shift;
    my $project=shift;

    my @genomefiles=`ls genome/*.gtf.gz`;

    foreach my $filename (@genomefiles) {
	chomp($filename);
	my $filepath=$project.'/'.$filename;
	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['genome.gtf',$filename,$md5sum,$filepath];
    }

    my @transfiles=`ls transcriptome/*.gtf.gz`;

    foreach my $filename (@transfiles) {
	chomp($filename);
	my $filepath=$project.'/'.$filename;
	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['transcriptome.gtf',$filename,$md5sum,$filepath];
    }

    my @juncfiles=`ls junctions/*.gtf.gz`;

    foreach my $filename (@juncfiles) {
	chomp($filename);
	my $filepath=$project.'/'.$filename;
	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['junctions.gtf',$filename,$md5sum,$filepath];
    }

    my @splitfiles=`ls splitmapping/*.gtf.gz |grep -v 'breakdown'`;

    foreach my $filename (@splitfiles) {
	chomp($filename);
	my $filepath=$project.'/'.$filename;
	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['split.gtf',$filename,$md5sum,$filepath];
    }

    my @clusterfiles=`ls clusters/*.gtf.gz`;

    foreach my $filename (@clusterfiles) {
	chomp($filename);
	my $filepath=$project.'/'.$filename;
	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['cluster.gtf',$filename,$md5sum,$filepath];
    }
}

sub get_bed_files {
    my $files=shift;
    my $fileprefix=shift;
    my $project=shift;
    my $tmpdir=shift;

    my @bedfiles=`ls $tmpdir/*.bed`;

    foreach my $filename (@bedfiles) {
	chomp($filename);
	my $filepath=$project.'/'.$filename;
	$filepath =gzip_file($filepath);

	my ($file,$md5sum)=get_md5sum($filepath);
	push @{$files},['bed',$filename.'.gz',$md5sum,$filepath];
    }
}

sub gzip_file {
    my $infile=shift;
    my $level=shift || 7;

    unless (-r $infile) {
	print STDERR "Can't find $infile\n";
	return(0);
    }
    my $outfile=$infile.'.gz';

    if (-r $outfile) {
	print STDERR $outfile,"\tIs present. Skipping...\n";
    } else {
	my $command="gzip -$level -c $infile > $outfile";
	run_system_command($command);
    }

    return($outfile);
}
