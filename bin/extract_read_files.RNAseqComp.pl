#!/soft/bin/perl
# DGK

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
# experiments belonging to it build a table with four columns, the file
# location for each of the different read files in the experiment, the lane name
# the pair name or lane name if unpaired, and the experiment name

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 qw(get_dbh read_config_file read_file_list);
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
my $tabsuffix='read_stats';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('extract_read_files.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			'',
			1)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# For each of tables extract the RPKMs of interest and get for each of the
# tables the different samples present in them
my %samples=%{get_read_stats_info(\%tables,
				  $dbh)};

foreach my $exp (keys %samples) {
    # get the readlist info
    my $pair=get_read_list_info($samples{$exp}->[2],
				$samples{$exp}->[3]);

    # get the link info
    my $readfile=$samples{$exp}->[2].'/readData/'.$samples{$exp}->[0].'.gz';
    my $source=get_link_source($readfile);
    print join("\t",
	       $exp,
	       @{$samples{$exp}},
	       $pair,
	       $readfile,
	       $source),"\n";
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
	
	$query ='SELECT Filename,Project,LaneName ';
	$query.="FROM $table";
	
	$sth = $dbh->prepare($query);
	my $count=$sth->execute();
	while (my ($file,$projdir,$lane)=$sth->fetchrow_array()) {
	    my $sample_id=join('_sample_',
			       $table,
			       $lane);
	    $samples{$sample_id}=[$file,$exp,$projdir,$lane];
	}	
    }
    return(\%samples);
}


