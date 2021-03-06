#!/soft/bin/perl
# DGK

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
# This script should take a project from the database and for each of the
# experiments belonging to it build a file with the number of uniquely mapped
# reads and total mapped reads at each stage

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file','read_file_list');
use RNAseq_pipeline_comp3 ('get_tables','check_tables');
use Getopt::Long;
use File::Find;
use Cwd;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $header;
my @tabsuffixes=('read_stats');
my $fraction='';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'limit|l=s' => \$fraction,
	   'header|h' =>\$header);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('extract_read_stats_info.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

my %samples;
foreach my $tabsuffix (@tabsuffixes) {
    # Get the tables belonging to the project
    my %tables=%{get_tables($dbhcommon,
			    $project,
			    $tabsuffix,
			    $fraction)};

    # Remove any tables that do not exist
    check_tables($dbh,
		 \%tables);
    
    # For each of tables extract the RPKMs of interest and get for each of the
    # tables the different samples present in them
    $samples{$tabsuffix}=get_read_stats_info(\%tables,
					     $dbh);

}


foreach my $tab (@tabsuffixes) {
    foreach my $exp (keys %{$samples{$tab}}) {
	my ($exp_id_all)=split('_sample_',$exp);
	my ($proj_id,$exp_id)=split('_',$exp_id_all,2);
	$exp_id=~s/_read_stats//o;
	
	print join("\t",
		   $proj_id,
		   $exp_id,
		   $tab,
		   @{$samples{$tab}{$exp}}),"\n";
    }
}
exit;

sub get_link_source {
    my $dir=shift;
    my $realpath='?';
    find sub {
	my @right = split /\//, $File::Find::name;
	
	my @left = do {
	    @right && ($right[0] eq "") ?
		shift @right :            # quick way
		split /\//, $dir;
	};
	
	while (@right) {
	    my $item = shift @right;
	    next if $item eq "." or $item eq "";
	    if ($item eq "..") {
		pop @left if @left > 1;
		next;
	    }
	    my $link = readlink (join "/", @left, $item);
	    if (defined $link) {
		my @parts = split /\//, $link;
		if (@parts && ($parts[0] eq "")) { # absolute
		    @left = shift @parts;   # quick way
		}
		unshift @right, @parts;
		next;
		
	    } else {
		push @left, $item;
		next;
	    }
	    
	}
	$realpath=join("/", @left);
	if ($dir ne $realpath) {
	    print STDERR "$File::Find::name is really $realpath\n";
	}
    }, ($dir);
    return($realpath)
}

sub get_read_list_info {
    my $projdir=shift;
    my $lane=shift;
    my $readlist=$projdir.'/read.list.txt';

    my $pair='?';

    if (-r $readlist) {
	my %files=%{read_file_list($readlist)};
	
	foreach my $file (keys %files) {
	    if ($files{$file}->[1] eq $lane) {
		$pair=$files{$file}->[0];
	    }
	}
    } else {
	warn "Can't find $readlist. Skipping...\n";
    }

    return($pair);
}

sub get_read_stats_info {
    my $tables=shift;
    my $dbh=shift;

    my %samples;

    foreach my $exp (keys %{$tables}) {
	my ($query,$sth);
	my $table=$tables->{$exp};
	
	$query ='SELECT sum(TotalReads),sum(NoAmbiguousBases),';
	$query.='sum(AmbiguousBases),sum(UniqueReads) ';
	$query.="FROM $table ";
	
	$sth = $dbh->prepare($query);
	my $count=$sth->execute();
	my ($total,$noamb,$amb,$unique)=$sth->fetchrow_array();
	$samples{$table}=[$total,$noamb,$amb,$unique];
    }
    return(\%samples);
}


