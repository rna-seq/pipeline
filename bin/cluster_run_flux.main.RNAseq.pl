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
# This script will take as an argument the name of a script. This must be a
# single word, and it will run that script on the cluster

use RNAseq_pipeline3 qw(get_fh run_system_command);
use RNAseq_pipeline_settings3 qw(read_config_file send2cluster);

# Declare some variables
my $bindir;
my $script;
my $queue;
my $projdir;
my $logsdir;

# get some options from the config file
my %options=%{read_config_file()};
$bindir=$options{'BIN'};
$projdir=$options{'PROJECT'};
$queue=$options{'CLUSTER'};
$logsdir=$options{'LOGS'};

$script=shift;

unless($script) {
    die "No input script supplied\n";
}

if ($queue) {
    my ($subfile,$jobname)=build_script_submission($script,
						   $bindir);    
    send2cluster($subfile,
		 $queue,
		 $jobname);

    # Clean up
    my $command="rm $subfile";
    run_system_command($command);
} else {
    print STDERR "Unable to run this on the cluster. Running locally...\n";
    my $command=$bindir.'/'.$script;
    run_system_command($command);
}

exit;

sub build_script_submission {
    my $script=shift;
    my $bindir=shift;

    unless ($bindir) {
	die "I don't know where to find the binaries\n";
    }
     
    unless($script) {
	die "No input supplied\n";
    }

    my $jobname=$script;
    $jobname=~s/(.pl)?$//;
    $jobname=~s/_//g;
    
    # Print the submission file
    my $subfile="subfile.$$.job";
    my $outfh=get_fh($subfile,1);
    
    print STDERR "Building submission file for $jobname...";

    print $outfh <<FORMEND;
# Get the job name
#\$ -N $jobname
    
# Write in to the current working directory
#\$ -cwd 
export PATH=\$PATH:/soft/bin
export LC_ALL='C'
#\$ -l h_vmem=32G

$bindir/$script
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile,$jobname);
}

__END__

=head1 NAME
    
    ?
    
=head1 SYNOPSIS
    
sample [options] 
    
  Options:
    -help            brief help message
    -man             full documentation
        
=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.
    
=back
    
=head1 DESCRIPTION
    
    This program is not documented yet, sorry

=cut
