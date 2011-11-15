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
# This script should extract from the database the total number of reads
# existing in the initial file and the total number of mapped reads from the
# mapping results table.
# It will compare these numbers for each of the lanes and if in any case they
# do not coincide it will print out the cases that do not coincide and die.
# It could maybe make the pipeline run again, but I don't know how stable this
# would be (if it starts calling itself recursively...)

use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_dbh);
use Getopt::Long;

my $file_list;
my $prefix;
my $type;
my $dbh;

GetOptions('type=s' => \$type);

my %options=%{read_config_file()};
$file_list=$options{'FILELIST'};
$prefix=$options{'PREFIX'};

my $readstab=$prefix.'_read_stats';
my $mappingstab=$type;

# Connect to the database
$dbh=get_dbh();

# get the numbers from the reads table
my %reads=%{get_reads($readstab,
		      $dbh)};

# get the reads from the mapping table
my %mapped=%{get_mapped($mappingstab,
			$dbh)};

# Close the dbh connection
$dbh->disconnect();

# Check if there are any differences
my @differences;
foreach my $lane (keys %reads) {
    unless ($mapped{$lane}) {
	push @differences, [$lane,'absent'];
	next;
    }
    unless ($mapped{$lane} == $reads{$lane}) {
	push @differences, [$lane,
			    "mapped: $mapped{$lane}",
			    "existing: $reads{$lane}"];
	next;
    }
}

if (@differences) {
    foreach my $diff (@differences) {
	print join("\t",
		   @{$diff}),"\n";
    }
    die "INCOMPLETE MAPPING\n";
}

exit;

sub check_table {
    my $table=shift;

    my %options=%{read_config_file()};

    my $dbh=get_dbh();
    my $present=0;

    my ($query,$sth);
    $query ='SELECT count(*) ';
    $query.="FROM $table";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    my $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    $sth->finish();
    $dbh->disconnect();
    
    if (@{$results}) {
	$present=1;
    } 

    return($present);
}

sub get_reads {
    my $table=shift;
    my $dbh=shift;

    print STDERR "Extracting info from  $table\n";

    my %reads;
    my $present=check_table($table);

    if ($present) {
	my ($query,$sth);
	$query ='SELECT LaneName, TotalReads ';
	$query.="FROM $table";
	$sth=$dbh->prepare($query);
	$sth->execute();

	while (my ($lane,$reads) =$sth->fetchrow_array()) {
	    $reads{$lane}=$reads;
	}
	$sth->finish();
    } else {
	die "$table is not present in the database\n";
    }

    return(\%reads);
 }   

sub get_mapped {
    my $table=shift;
    my $dbh=shift;

    print STDERR "Extracting info from  $table\n";

    my %mapped;
    my $present=check_table($table);

    if ($present) {
	my ($query,$sth);
	$query ='SELECT LaneName, TotalReads ';
	$query.="FROM $table";
	$sth=$dbh->prepare($query);
	$sth->execute();

	while (my ($lane,$reads) =$sth->fetchrow_array()) {
	    $mapped{$lane}=$reads;
#	    print STDERR join("\t",
#			      $lane,$reads),"\n";
	}
	$sth->finish();
    } else {
	die "$table is not present in the database\n";
    }

    return(\%mapped);
 }  
