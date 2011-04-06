#!/usr/bin/perl

use strict;
use warnings;

# Objective
# This script should take a file containing a list of directories and a script
# name. The script mus be run without any arguments (so everything should be
# self contained, maybe using a wrapperif arguments are necessary)
# The script will be run using the available resources by queuing up the
# different instances in as many parallel runs as are specified by the -threads
# option

# Load modules
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil);
use Cwd;
use RNAseq_pipeline3 qw(get_fh run_system_command);

# Daclare variables
my $help;
my $man;
my $threads=2;
my $directories;
my $script;
my $debug=0;

# Get input
GetOptions('tabs|threads|t=s' => \$threads,
	   'directories|f=s' => \$directories,
	   'script|s=s' => \$script,
	   'help|h' => \$help,
	   'debug|d' => \$debug,
	   'man|m' => \$man);
 
# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage("No list of project directories suplied") unless ($directories);
pod2usage("No script suplied") unless ($script);

# Get available clients
### TO DO
# Here we should check the local load before starting anything in case we
# cannot submit as many simultaneous jobs as we intended

# Check if the number of threads requested is too large. Here we will have the
# number of available clients at some stage, but currently 8 is fine.
if ($threads > 8) {
    print STDERR "Are you sure you want to use $threads threads?(y/N)";
    my $yes=<STDIN>;
    chomp($yes);
    if ($yes=~/^y/i) {
	print STDERR "Well...as you wish, but don't you say I didn't warn you...\n";
    } else {
	print STDERR "Ok, in that case I'll let you reconsider...\n";
	exit;
    }
}

# Get directories from file
my @dirs=();
my $dirfh=get_fh($directories);
while (my $dir=<$dirfh>) {
    chomp($dir);
    push @dirs, $dir;
    if ($debug &&
	(@dirs >=6)) {
	last;
    }
}
close($dirfh);
my $dir_number=@dirs;

# Count the number of jobs
print STDERR $dir_number,"\n";

# If we have more threads than jobs reduce the threads to the job number
if ($threads > $dir_number) {
    print STDERR "Reducing threads to $dir_number\n";
    $threads=$dir_number;
}

# Determine the number of sequences per job
my $dirs_per_thread=ceil($dir_number/$threads);
print STDERR $dirs_per_thread, "\tJobs will be sent to each thread\n";

my %tmp_jobs;
for (my $i=0;$i < $threads;$i++) {
    @{$tmp_jobs{$i}}=splice(@dirs,0,$dirs_per_thread);
}

# Start the jobs
my $cwd=getcwd();
run_script($script,\%tmp_jobs,$cwd);

print STDERR "All jobs submitted\n";

exit;

# Subroutines
sub run_script {
    my $script=shift;
    my $jobs=shift;
    my $cwd=shift;
    my $exit_status=0;

    my $command;
    $command ="gnome-terminal ";
    foreach my $terminal (keys %{$jobs}) {
	$command.=" --tab -e 'bash -c \"echo $terminal;";
	my @directories=@{$jobs->{$terminal}};
	foreach my $dir (@directories) {
	    $dir=$cwd.'/'.$dir;
	    $command.="cd $dir;echo \$PWD;nice $script;";
	}
	$command.="echo done;exec bash\"'";
    }
    run_system_command($command);

    return($exit_status);
}

=head1 NAME
    
    multirun.lp run a script/command on all directories listed in a file in a 
    new window. The new window will be opened with as many tabs as indicated
    in the tabs option with one thread running in each of them
    
=head1 SYNOPSIS
    
multirun.pl -f directory.list.txt -s 'command/script' -t threads
    
  Options:
    -help            brief help message
    -man             full documentation
    -tabs|threads|t  Number of tabs to use. Each tab will correspond to one
                     thread running in parallel. Default 2.
    -directories|f   List of directories on which to run the provided script.
                     (mandatory)
    -script|s        Script or command to run in each directory.(mandatory)
    -debug|d         Debugging
        
=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.

=item B<-tabs|threads|t>
    
    A new gnome terminal window will be open containing as many tabs as
    specified by the t option. The jobs will be distributed between these tabs
    which correspond to one thread each, and run sequentially within them using
    nice with the default priority. This should prevent  clogging up things too
    much if a very large number of threads is specified. In any case a warning
    will be issued if more then 8 threads are requested.

=item B<-directories|f>
    
    This should be a file containing a list of directories one in each line.
    Currently these directories are expected to stem from the directory where
    the script is run and no checking is made for their existence before trying
    to use them

=item B<-script|s>
    
    The script or command to be run, this can be more than one word if given
    within quotes
    
=back
    
=head1 DESCRIPTION
    
    The script given is run in each of the directories listed by specifically
    changing the working directory to that directory and executing the script
    which means that it can use files in these directories and any other things
    that may have standard naming.

=cut
