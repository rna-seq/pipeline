#!/usr/bin/perl
#
# DGK 2010
#

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
    unshift @INC, "$libdir";
}

use Getopt::Long;
use DBI;
use Cwd;
use RNAseq_pipeline_start3 qw(get_tables_hash);
use RNAseq_pipeline3 qw( MySQL_DB_Connect);

# Set the global options
# Read any command line arguments
my $pipfile = 'RNAseq_pipeline.txt';	#  -f   control file name

my @targets;				# command line args -- names of stages to "build"

# Set some other global variables:
my $dbh;				# database handle
my %rulelist; # Hash structured list of rules
my %phony;				# rule names that aren't database tables
my %symbols;				# symbol table

# Initialize the variables and set any command line options
initialize();

# Setup some environmental variables
$ENV{'PATH'} = $symbols{'BIN'}. ':' . $ENV{'PATH'};
$ENV{'PERL5LIB'} = $symbols{'LIB'}.':'.$ENV{'PERL5LIB'};
$ENV{'DB'}=$symbols{'DB'};

# Once all the variables have been initialized process the pipfile
processFile($pipfile,
	    \%rulelist);

# This must be done after the initialization step and the file has been
# processed to set all the variables correctly
#my @root_rules=@{get_root_rules(\%rulelist)};

# Get the database handle
$dbh=MySQL_DB_Connect($symbols{'DB'});
# Get table last modification dates
Get_table_Dates(\%rulelist);

# Get the actual tables we generate
my %tables=%{get_tables_hash($symbols{'PREFIX'})};

foreach my $rule (keys %tables) {
    my $rule_id=$tables{$rule};
    unless ($rulelist{$rule_id}) {
	next;
    }
    my $status=check_rule(\%rulelist,
			  $rule_id);
    print join("\t",
	       $rule_id,
	       $status),"\n";
}

exit;

# Check if the rule is up to date or not
sub check_rule {
    my $list=shift;
    my $rule=shift;

    # get the date at which the rule was completed
    my $date=$list->{$rule}->{'timestamp'};
    my $status='Finished';
    if ($date == 19000101010101) {
	$status='Pending';
    } else {
	my @requisites=@{$list->{$rule}->{'precursor'}};
	foreach my $req (@requisites) {
	    my $reqdate=$list->{$req}->{'timestamp'};
	    if ($reqdate &&
		($reqdate > $date)) {
		$status='Outdated';
	    }
	}
    }
 #   $list->{$req}->{'status'}=$status;

    return($status);
}

######################################################################
#
# Rule processing
#

# Return the current time as a MySQL DATETIME value.

sub now {
    my $time = `date +"%Y-%m-%d %H:%M:%S"`;
    chomp($time);
    return $time;
}

######################################################################
#
# Setup
#
# Get the modification date of a table.  If no table name is supplied,
# get the dates for all tables, in which case the return value is a
# hash associating a table name with its modification time.
sub getDate {
    my $tbl = shift;

    my $query = "SHOW TABLE STATUS";
    if ($tbl) {
	$query .= " LIKE '$tbl'";
    }

    my $dh = {};

    # Modified in order to be able to access the database after long queries that
    # will kill the $dbh
    my $sth;

    # Check if the database is still responding
    eval {
	$sth= $dbh->prepare("SHOW TABLE STATUS") or
	    die "PIP: error at prepare: $DBI::errstr\n";
	$sth->execute or
	    die "PIP: error at execute: $DBI::errstr\n";
    };

    # If not disconnect and reconnect
    if ($@) {
	$dbh->disconnect();
	print STDERR 'Reconnecting to the database...';
	$dbh=MySQL_DB_Connect($symbols{'DB'});
	print STDERR "done\n";

	# Try again and this time if the query failed it is true
	$sth= $dbh->prepare("SHOW TABLE STATUS") or
	    die "PIP: query failed: $DBI::errstr\n";
	$sth->execute or
	    die "PIP: query failed: $DBI::errstr\n";
    }

    while (my $href = $sth->fetchrow_hashref() ) {
	$dh->{$href->{'Name'}} = $href->{'Update_time'};
    }

    if ($tbl) {
	return $dh->{$tbl};
    }
    else {
 	return $dh;
    }
}

# Each target corresponds to a table in the project database.  Look up
# the modification date of each table and save it in the rule set.  If
# a table doesn't exist yet it will get an empty timestamp.
sub Get_table_Dates {
    my $rules = shift;

    my $dates = getDate();

    foreach my $rule (keys(%{$rules})) {
	next if $phony{$rule};
	$rules->{$rule}->{"timestamp"} = $dates->{$rule} || "1900-01-01 01:01:01";
	$rules->{$rule}->{"timestamp"}=~s/[^\d]//g;
    }
}

sub rules {
    my $line=shift;
    my $rules=shift;
    my $lines=shift;
    my ($target,$prec)=split(':',$line);
    $prec=process_line($prec);
    my @prec=split(/\s+/,$prec);
    $rules->{$target}->{'precursor'}=[@prec];
    
    while (my $command=shift(@{$lines})) {
	# Check if it is actually a command and if not exit loop
	# Basically if the line starts with anything except empty space
	# or a comment the line is a new rule
	if ($command !~ /^[\s\#]/) {
	    unshift @{$lines},$command;
	    last;
	}
    }
}

sub substitute_values {
    my $string=shift;

    my @parts=split(/\$/,$string);
    for (my $i = 1; $i < @parts; $i++) {
	my ($var) = ($parts[$i]=~/(\w+)/);
	if ((exists $symbols{$var}) && 
	    (defined $symbols{$var})) {
	    $parts[$i] =~s/\w+/$symbols{$var}/;
	} else {
	    $parts[$i] = "\$".$parts[$i];
	}
    }
    return(join('',@parts));
}

sub symbol {
    my $line=shift;
    # Get the varible name and value and also remove preceding and trailing whitespace
    my ($var,$val)=($line=~/\s*(\w+)\s*=\s*(.*)\s*$/);
    my $substituted=substitute_values($val);

    $symbols{$var}=$substituted;
}

# Get the next noncomment line from the array of lines.  Some details:
# * comments and empty lines are skipped, and comments are stripped
#   from the ends of lines
sub process_line {
    my $line=shift;

    # Remove comment lines
    $line=~s/\#.*//;
    # Remove trailing whitespace
    $line=~s/\s*$//;
    # Remove leading whitespace
    $line=~s/^\s*//;

    return($line);
}

# Read the command file and save a list of rules and variables.  The
# three types of file entries are:
# PIP directives (lines starting with periods)
# Variable definitions (single lines of the form "X = Y")
# Rule definitions.
sub processFile {
    my $filename=shift;
    my $rules=shift;
    my $rules_fh;

    # Get the information out of the pipfile
    open($rules_fh,$filename) ||
	die "Can't open $filename:$!\n";
    my @lines=<$rules_fh>;
    close($rules_fh);
    
    while (my $line=shift(@lines)) {
	chomp($line);

	my $line=process_line($line);
	# Skip lines that are empty after processing
	if ($line=~/^\s*$/) {next;}
	elsif ($line=~/^\./) {
	    # Get directives (lines starting with dot)
	    next;
	} elsif ($line=~/\w+\s*=/) {
	    # Symbol definitions
	    symbol($line);
	} elsif ($line=~/^\w+\s*:/) {
	    # Rule definition
	    rules($line,$rules,\@lines);
	} else {
	    die "Pip syntax error: $line\n";
	}
    }    
}

# Initialize the variables
sub initialize {
    GetOptions(
	"file:s"  => \$pipfile,
    );
    # Check for option errors
    exit if $Getopt::Long::error;

    # Get the rules to execute
    @targets = @ARGV;

    # Set the defaults for the variables in the pipfile
    chomp(my $projdir = `pwd`);		# main project directory
    my $bindir = "$projdir/bin";	# project applications
    my $sqldir = "$projdir/mysql/table_build";	# project tables
    my $datadir= "$projdir/data"; # Bulk data
    my $logsdir= "$projdir/logs"; # Logs from the scripts
    my $libdir= "$projdir/lib"; # modules required by the project scripts
    my $tabdatdir="$projdir/mysql/table_data"; # Compressed copy of the database tables
    my $graphsdir="$projdir/graphs"; # Any pictures or plots from the project
    my $database=(split(/\//,$projdir))[-1];

    %symbols = (
		"PROJECT" => $projdir,
		"PATH"    => $bindir, # In case there is a problem 
		"BIN"    => $bindir,
		"TABLES"  => $sqldir,
		"DB"      => $database,
		'DATA' => $datadir,
		'LOGS' => $logsdir,
		'LIB' => $libdir,
		'TAB_DAT' => $tabdatdir,
		'GRAPHS' => $graphsdir,
		'DEBUG' => 0
	);
}

