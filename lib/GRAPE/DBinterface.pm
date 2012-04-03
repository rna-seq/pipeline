package GRAPE::DBinterface;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export subroutines to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
push @EXPORT_OK,('MySQL_DB_Connect');


use strict;
use warnings;

use DBI;
use DBD::mysql;

sub MySQL_DB_Connect {
    my $database = shift;

    my $datasource = "DBI:mysql:$database";
    my $dbh;
    my $cnf_file='.my.cnf';

    # Check for .my.cnf as this will contain all required information for the
    # database connection
    if (-e "$ENV{HOME}/$cnf_file") {
	$datasource .= ";mysql_read_default_file=$ENV{HOME}/$cnf_file";
	$dbh = DBI->connect($datasource, undef, undef, {RaiseError => 1});
    } else {
	print STDERR "Unable to find $ENV{HOME}/.my.cnf\n";
    }
    return $dbh;
}

sub check_key_value {
    my $dbh=shift;
    my $table=shift;
    my $field=shift;
    my $key=shift;
    
    my $field_value='';

    my ($query,$sth,$count);
    $query ="SELECT $field ";
    $query.="FROM $table ";
    $query.="WHERE $key = ?";
    
    $sth=$dbh->prepare($query);
    $count=$sth->execute($key);
    
    if ($count > 1) {
	# There is a problem
	die "Entry $key is present more than once in $table\n";
    } else {
	($field_value)=$sth->fetchrow_array();
    }
    return($field_value);
}

sub set_species_field {
    my $dbh=shift;
    my $key=shift;
    my $table=shift; # species_info
    my $field=shift; # species_id

    my $value=check_field_value($dbh,
				$table,
				$field,
				$key);
    
    unless ($value) {
	# The entry is absent and we must set it
	# Get the genus and abbreviation
	my ($genus,$specific)=split(/\s+/,$key,2);
	my $abbreviation=join('',
			      substr($genus,0,1),
			      substr($specific,0,3));
	
	# Insert the info into the database
	my $query;
	$query ="INSERT INTO $table ";
	$query.='SET species = ? , genus = ? , abbreviation = ?, sp_alias = ? ';
	print STDERR "Executing: $query\n";
	my $sth2=$dbh->prepare($query);
	$sth2->execute($key,$genus,$abbreviation,'-');
    }
    return($value);
}

1;
