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
# This script will check if the required scripts are present in the path

# Load some modules
use RNAseq_pipeline3 qw(get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file);

# Get some options from the configuration file
my $prefix;
my $bindir;
my $debug;

my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$bindir=$options{'BIN'};

$ENV{'PATH'}.=":$bindir";

my $log_fh=get_log_fh('check_prerequisites.RNAseq.log',
		      $debug);

my @prerequisites=('gem-mapper',
		   'gem-split-mapper',
		   'gem-do-index',
		   'overlap',
		   'flux.sh');

print $log_fh "Checking prerequisites\n";
foreach my $req (@prerequisites) {
    my $absent=check_prerequisite($req);

    if ($absent) {
	die "I cannot find required executable $req\n";
    } else {
	print $log_fh join("\t",
			   $req,
			   'present'),"\n";
    }
}

exit;

sub check_prerequisite {
    my $program=shift;

    my $command="which $program > /dev/null";
    my $ok=system($command);

    print STDERR join("\t",
		      $program,
		      $ok),"\n";

    return($ok);
}
