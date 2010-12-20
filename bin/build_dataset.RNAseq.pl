#!/soft/bin/perl

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
# This script should read the file_list file, and it will create a table
# with the contents of this file. It will also do some quality checks over it
# to make sure all the formatting, etc.. is correct.

# Load some modules
use RNAseq_pipeline3 qw(get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_dbh);

# Get some options from the configuration file
my $file_list;
my $prefix;
my $debug=1;

my %options=%{read_config_file()};
$file_list=$options{'FILELIST'};
$prefix=$options{'PREFIX'};

my $log_fh=get_log_fh('build_dataset.RNAseq.log',
		      $debug);
my %files=%{read_file_list()};
my %pairs;

print $log_fh "Checking file $file_list\n";
foreach my $readfile (keys %files) {
    # Check if the number of fields is correct;
    my @fields=@{$files{$readfile}};
    if (@fields < 2) {
	die "Fields missing in $file_list\n";
    } elsif (@fields > 3) {
	die "Too many fields in $file_list\n";
    } else {
	push @{$pairs{$fields[0]}},$fields[1];
    }
}

# Get from the database the number of reads in each lane
my $readstab=$prefix.'_read_stats';
my $dbh=get_dbh();
my %reads=%{get_read_number($readstab,
			    $dbh)};
$dbh->disconnect();

# Check the pairing is correct
foreach my $pair (keys %pairs) {
    my @lanes=@{$pairs{$pair}};
    if (@lanes == 1) {
	print $log_fh "$pair reads identified as single\n";
	unless ($pair eq $lanes[0]) {
	    die "In single reads the pair id and lane Id should be the same\n";
	}
    } elsif (@lanes == 2) {
	print $log_fh "$pair reads identified as paired\n";
	unless ($reads{$lanes[0]}==$reads{$lanes[1]}) {
	    die "In paired reads both read files should have the same number of entries.\nThis is not the case for pair $pair\n";
	}
    } else {
	die "The format of $file_list is not correct\n";
    }
}

# Everything seems OK, so we will write a file
foreach my $file (keys %files) {
    my $group=$files{$file}->[2] || 'All';
    my $pair=$files{$file}->[0];
    my $lane=$files{$file}->[1];
    my $lane_id=$pair;
    $lane_id=~s/[^\d\w]//g;
    print join("\t",
	       $lane_id,
	       $lane,
	       $pair,
	       $group),"\n";
}
close($log_fh);

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

sub get_read_number {
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
