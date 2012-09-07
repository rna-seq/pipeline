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
push @EXPORT_OK,('MySQL_DB_Connect',
		 'check_field_value','check_key_value',
		 'check_table_existence',
		 'create_MySQL_table');


use strict;
use warnings;

use DBI;
use DBD::mysql;

# Load GRAPE modules
use GRAPE::Base ('get_fh','run_system_command');

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

# Create a table having the 
sub create_MySQL_table {
    my $database=shift;
    my $table=shift;
    my $tablebuild=shift;

    print STDERR "Creating $table...";
    my $file_name=$table.'.sql';
    my $table_fh=get_fh($file_name,1);
    print $table_fh "SET FOREIGN_KEY_CHECKS=0;\n";
    print $table_fh $tablebuild,"\n";
    print $table_fh "SET FOREIGN_KEY_CHECKS=1;\n";
    close($table_fh);
    
    my $command="mysql $database < $file_name";
    run_system_command($command);
    $command="rm $file_name";
    run_system_command($command);
    print STDERR join("\t",
		      $table,
		      "Generated"),"\n"
}

# Check if a table exists in the database
sub check_table_existence {
    my $dbh=shift;
    my $table=shift;

    print STDERR "Checking database for $table...";

    my ($query,$sth);
    $query ='SELECT count(*) ';
    $query.="FROM $table";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    my $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    my $present=0;
    
    if (@{$results}) {
	# Print the table location for the junctions of this experiment
	print STDERR "$table is present\n";
	$present=1;
    } else {
	print STDERR "$table is absent\n";
    }
    return($present);
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

# Get the value of a certain filed to see if it is defined
sub check_field_value {
    my $dbh=shift;
    my $reference=shift;
    my $table=shift;
    my $field=shift;

    # Check if the table actually exists
    my $present=check_table_existence($dbh,
				      $table);
    die "Unable to search a nonexistant table\n" unless $present;

    my ($key,$value)=@{$reference};
    my ($query,$sth,$count);
    $query ="SELECT $field ";
    $query.="FROM $table ";
    $query.="WHERE $key = ?";
    $sth=$dbh->prepare($query);

    # A sth is returned here instead of the results as this allows for the
    # query to be executed again with less code
    return($sth);
}

1;
