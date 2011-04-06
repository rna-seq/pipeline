#!/usr/bin/perl
#
#    execute_RNAseq_pipeline.pl
#    Copyright (C) 2009-2011  DGK (Based on pip by John Conery)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

use IO::Handle;
use Getopt::Long;
use Pod::Usage;
use DBI;
use Cwd;
use RNAseq_pipeline3 qw(MySQL_DB_Connect get_log_fh);

# Set some general variables that will be used through the script:

# Set a signal handler to kill any process started by the script if it
# receives the kill signal.

# I don't know if this will actually work but its a pain when all the children
# keep executing after the main program dies
my ($kill_all_sub,$add_proc_sub)=\&set_process_management;
$SIG{'INT'}=\&{$kill_all_sub};

# Set the global options
# Read any command line arguments
my $help;				#  -h   print help and exit
my $pipfile = 'RNAseq_pipeline.txt';	#  -f   control file name
my $nomake;				#  -n   print commands, don't execute
my $silent;				#  -s   suppress log messages
my $debug;				#  -d   print trace (debug mode)
my $force;				#  -r   restart (force execution of
                                        #       stage)
### TO DO
my $startpoint;                         #  -b   base Set the base rule from
                                        #       which to start

my $log_fh=get_log_fh('execute_pipeline.err');
unless($debug) {
    STDERR->fdopen( $log_fh,  'w' ) or die $!;
}

my $success=GetOptions(
    "file:s"  => \$pipfile,
    'n'       => \$nomake,
    "help"    => \$help,
    "silent"  => \$silent,
    "debug"   => \$debug,
    'restart' => \$force,
    'startpoint|b=s' => \$startpoint
    );
# Check for option errors
pod2usage(1) unless $success;

# Set some variables to override each other
$silent = 0 if $debug;		# -debug overrides -silent
$silent = 0 if $nomake;		# -nomake overrides silent

# If silent, set STDOUT to /dev/null
if ($silent) {
    open(TRASH,"> /dev/null");
    *STDOUT = *TRASH;
}

my @targets = @ARGV;			# command line args: rules to run

# Set some other global variables:
my $dbh;				# database handle
my %rulelist; # Hash structured list of rules
my %phony;				# rule names that aren't database tables
my %symbols;				# symbol table

# If help option is used print the help
pod2usage(1) if $help;

# Set the defaults for the variables in the config file
chomp(my $projdir = `pwd`);		   # main project directory
my $bindir = "$projdir/bin";	           # script directory
my $sqldir = "$projdir/mysql/table_build"; # project tables
my $logsdir= "$projdir/logs";              # Logs from the scripts
my $tabdatdir="$projdir/mysql/table_data"; # Compressed copy of the database tables
my $graphsdir="$projdir/graphs";        # Any pictures or plots from the project
my $database=(split(/\//,$projdir))[-1];

# Fill in the symbols hash
%symbols = (
    "PROJECT" => $projdir,
    "PATH"    => $bindir, # In case there is a problem 
    "BIN"    => $bindir,
    "TABLES"  => $sqldir,
    "DB"      => $database,
    'LOGS' => $logsdir,
    'TAB_DAT' => $tabdatdir,
    'GRAPHS' => $graphsdir,
    'DEBUG' => 0
    );

# Set the debugging options
if ($debug) {
    $symbols{'DEBUG'}=1;
    print "Config:\n";
    foreach my $var (keys %symbols) {
	print join("\t",
		   $var,
		   $symbols{$var}),"\n";
    }
}

# Setup some environmental variables
$ENV{'PATH'} = $symbols{'BIN'}. ':' . $ENV{'PATH'};
$ENV{'DB'}=$symbols{'DB'};

# Once all the variables have been initialized process the pipfile
processFile($pipfile,
	    \%rulelist);

# Get a list of common subroutines we will need.
# This must be done after the initialization step and the file has been
# processed to set all the variables correctly
my $common_subs=create_common_subs();
my @root_rules=@{get_root_rules(\%rulelist)};

# Get the database handle
$dbh=MySQL_DB_Connect($symbols{'DB'},
		      $symbols{'HOST'});

# Get table last modification dates
Get_table_Dates(\%rulelist);

# If no target rule is given set the targets to the root rules
if (@root_rules && (@targets==0)) {
    @targets=@root_rules;
}

# Print a message for the user
print "RNAseq analysis workflow started ", `date`;
print "cwd: ", `pwd`;
print "Command file: $pipfile\n";

# Execute the rules
foreach my $rule (@targets) {
    execute_rule($rule,0);				# 0 means top level call
}

# Print an ending message
print "\nRNAseq analysis workflow completed ", `date`;

close($log_fh);

exit;

######################################################################
#
# Rule processing
#

# Execute a rule of the pipeline.  The first param is the name of the
# rule, the second is a call level used in pretty-printing status
# messages.  The return value is the new timestamp for the table.
sub execute_rule {
    my $rule=shift;
    my $level=shift;

    print "\ncall: $rule\n" if $debug;

    # Check if the rule actually exists
    unless ($rulelist{$rule}) {
	print STDERR "*** pip: no rule for $rule\n";
	return 0;
    }

    my $prec = $rulelist{$rule}->{"prec"};
    my $time = $rulelist{$rule}->{"timestamp"} || 0;
    my $body = $rulelist{$rule}->{"body"};

    my $ok = 1;
    my $execute = 0;
    # Set default return value to the current date
    my $return_value = $time;

    # always execute the body if either:
    #   the precedent list is empty
    #   the -r flag was given on the command line
    # otherwise skip the body if the (first) precedent is "@"
    # otherwise call make recursively on each precedent, and execute
    # the body if any of:
    #   this target is a phony rule
    #   the target is not phony but the table doesn't exist
    #   any precedent is newer than this target
    # note: don't check phony rules in the first test because phony
    #   rules can have precedents and we want to call those rules....
    # note: a rule will either be phony (always run it) or have an @
    #   (only run if forced)
    # note: when a rule is forced its precedents are not checked

    if (@$prec == 0 || $force) {
	$execute = 1;
    } elsif ($$prec[0] eq "@") {
	$execute = ($time eq 0);
    } else {
	$execute = ($phony{$rule} || $time eq 0);
	foreach my $prec (@$prec) {
	    my $time_prec = execute_rule($prec,$level+1);
	    if ($time_prec gt $time) {
		$execute = 1;
	    }
	    print "$rule [$time] vs $prec [$time_prec]: $execute\n" if $debug;
	}
    }

    if ($execute) {
	print "\nExecuting $rule\n";
	$ok = process_rule($rule,$body,$level+1);
	$return_value = update_time_stamp($rule);
    }    else {
	print "\n$rule is up to date\n" if ($level == 0 || $debug);
    }

    # Check if the execution has been correct
    if ($ok) {
	return $return_value;
    } else {
	print STDERR "*** ERROR: rule for $rule failed\n";
	exit 0;
    }
}

# Return the current time as a MySQL DATETIME value.
sub now {
    my $time = `date +"%Y-%m-%d %H:%M:%S"`;
    chomp($time);
    return $time;
}

# Update the timestamp for a rule X.  If X is a phony rule, or if
# the "nomake" flag is set, simulate the update by getting the current
# time, otherwise access the table's real timestamp in the database.
sub update_time_stamp {
    my $rule = shift;;

    if ($phony{$rule}) {
	return now();
    }
    elsif ($nomake) {
	return ($rulelist{$rule}->{"timestamp"} = now());
    }
    else {
	return ($rulelist{$rule}->{"timestamp"} = getDate($rule));
    }
}

# Execute the commands in the body of a rule.  The argument is a
# pointer to a list of commands.  Return 0 if the command fails
sub process_rule {
    my ($rule,$body,$level) = @_;
    my $ok = 1;

    print "\n$rule:\n";

    foreach my $line (@$body) {
	my $command = substitute_values($line);
 	print "Executing $command\n" unless ($command =~ /^echo/ );
	# Check if it is a key form the %common_subs hash
	if (exists $common_subs->{$command}) {
	    print "Substituting $command\n";
	    $common_subs->{$command}->(\$command,$rule);
	} 
	if (ref($command) eq 'ARRAY') {
	    foreach my $comm (@$command) {
		# bitwise adition, if both are 1 eq 1 if not 0I  think
		print "\t",$comm,"\n";
		next if $nomake;
		$ok &= (system($comm) == 0);
	    	last unless $ok;
	    }
	} else {
	    next if $nomake;
	    $ok &= (system($command) == 0);
	}
	last unless $ok;
    }

    return $ok;
}

# Create a hash with commonly used subs
sub create_common_subs {
    my %common_subs;
    my $database=$symbols{'DB'};
    my $table_build_dir=$symbols{'TABLES'};
    my $table_store_dir=$symbols{'TAB_DAT'};
    my $log_dir=$symbols{'LOGS'};
    $common_subs{'table_create'}= sub {
	my $command=shift;
	my $rule=shift;
	$$command="mysql $database < $table_build_dir/$rule.sql";
    };
    $common_subs{'table_import'}=sub {
	my $command=shift;
	my $rule=shift;
	$$command="mysqlimport -L $database $rule.txt";
    };
    $common_subs{'table_zip'}=sub {
	my $command=shift;
	my $rule=shift;
	$$command="gzip -9 $rule.txt";
    };
    $common_subs{'table_store'}=sub {
	my $command=shift;
	my $rule=shift;
	$$command="mv $rule.txt.gz $table_store_dir";
    };
    $common_subs{'store_log'}=sub {
	my $command=shift;
	my $rule=shift;
	my $short=$rule;
	$short=~s/[^_]+_//;
	$$command="mv build_${short}.log $log_dir";
    };
    $common_subs{'all_standard_table'}=sub {
	my $command=shift;
	my $rule=shift;
	my $short=$rule;
	$short=~s/[^_]+_//;
	$$command=["mysql $database < $table_build_dir/$rule.sql",
		   "mysqlimport -L $database $rule.txt",
		   "gzip -9 $rule.txt",
		   "mv $rule.txt.gz $table_store_dir",
		   "mv build_${short}.log $log_dir"];
    };
    $common_subs{'all_no_log_table'}=sub {
	my $command=shift;
	my $rule=shift;
	$$command=["mysql $database < $table_build_dir/$rule.sql",
		   "mysqlimport -L $database $rule.txt",
		   "gzip -9 $rule.txt",
		   "mv $rule.txt.gz $table_store_dir"]
    };

    return(\%common_subs);
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

    # Modified in order to be able to access the database after long queries
    # that will kill the $dbh
    my $sth;

    # Check if the database is still responding
    eval {
	$sth= $dbh->prepare("SHOW TABLE STATUS") or
	    die "Error at prepare: $DBI::errstr\n";
	$sth->execute or
	    die "Error at execute: $DBI::errstr\n";
    };

    # If not disconnect and reconnect
    if ($@) {
	$dbh->disconnect();
	print STDERR 'Reconnecting to the database...';
	$dbh=MySQL_DB_Connect($symbols{'DB'});
	print STDERR "done\n";

	# Try again and this time if the query failed it is true
	$sth= $dbh->prepare("SHOW TABLE STATUS") or
	    die "Error: query failed: $DBI::errstr\n";
	$sth->execute or
	    die "Error: query failed: $DBI::errstr\n";
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
	$rules->{$rule}->{"timestamp"} = $dates->{$rule};
	print "date($rule) = $$dates{$rule}\n" if $debug;
    }
}

sub get_root_rules {
    my $rules=shift;
    my @root_rules;

    foreach my $rule (keys %{$rules}) {
	unless (@{$rules->{$rule}->{'prec'}}) {
	    push @root_rules,$rule;
	}
    }
    return(\@root_rules);
}

sub rules {
    my $line=shift;
    my $rules=shift;
    my $lines=shift;
    my ($target,$prec)=split(':',$line);
    $prec=process_line($prec);
    my @prec=split(/\s+/,$prec);
    $rules->{$target}->{'prec'}=[@prec];
    
    while (my $command=shift(@{$lines})) {
	# Check if it is actually a command and if not exit loop
	# Basically if the line starts with anything except empty space
	# or a comment the line is a new rule
	if ($command !~ /^[\s\#]/) {
	    unshift @{$lines},$command;
	    last;
	}
	$command=process_line($command);
	# remove any left empty lines
	if ($command) {
	    push @{$rules->{$target}->{'body'}},$command;
	}
    }
}

# Replace in the string any variables that are in the symbol
# hash
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

# Save the name and value of the variable in the symbol hash
# after interpolating any variable contained in the value
sub symbol {
    my $line=shift;
    # Get the varible name and value and also remove preceding and trailing whitespace
    my ($var,$val)=($line=~/\s*(\w+)\s*=\s*(.*)\s*$/);
    my $substituted=substitute_values($val);

    $symbols{$var}=$substituted;

    if ($debug) {
	print join("\t",
		   $var,
		   $substituted),"\n";
    }
}

sub directive {
    my $line=shift;
    if ($line =~ /\.phony/i) {
	$line =~ s/^.*:\s*//;
	foreach my $x (split(/\s+/,$line)) {
	    $phony{$x} = 1;
	}
    } else {
	print "ERROR: syntax error:  unknown directive $line\n";
	exit 1;
    }
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
# Pipeline directives (lines starting with periods)
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
	    directive($line);
	} elsif ($line=~/\w+\s*=/) {
	    # Symbol definitions
	    symbol($line);
	} elsif ($line=~/^\w+\s*:/) {
	    # Rule definition
	    rules($line,$rules,\@lines);
	} else {
	    die "Template syntax error: $line\n";
	}
    }    
}

sub set_process_management {
    my @processes;

    my $add_proc_sub=sub {
	my $proc=shift;
	push @processes, $proc;
    };

    my $kill_all_sub=sub {
	foreach my $proc (@processes) {
	    print STDERR $proc,"\n";
	}
	kill 9,@processes;
	exit;
    };

    return($kill_all_sub,$add_proc_sub);
}

__END__

=head1 NAME
    
    execute_RNAseq_pipeline.pl [options] rulename
    
=head1 SYNOPSIS
    
sample [options] 
    
  Options:
    -help|h            brief help message
    -man|m             full documentation
    -file|f            Control file name. Defaults to RNAseq_pipeline.txt
    -n                 Print commands but don't execute them
    -silent|s          Silent
    -debug|d           Debug
    -restart|r         Restart at the specific rule (force execution of rule)
        
=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.

=item B<-n>
    
    Print the commands that would be executed to the command line without
    actually doing anything.

=item B<-restart|r>
    
    Restart at the specific rule. As it is restarting it will force the
    execution of the rule regardless of the presence of the prerequisites. As
    it does not chech for these prerequisites make sure they are there, if not
    the rule will fail.
    
=back
    
=head1 DESCRIPTION
    
    This program is not documented yet, sorry

=cut
