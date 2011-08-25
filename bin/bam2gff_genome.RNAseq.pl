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
# This script should parse a BAM file and create a gff filefrom it
# It will take the read.list file in order to know which are the files to
# process

# After this for each genome mapped file it will take the corresponding junction
# mapped file and filter out those possible pseudogene hits (those cases where
# a read hits both genome and junctions)

# It will print out the hits in gtf format separating them into two files, the
# unique mappings and the multimappings.

use Getopt::Long;
use Bio::DB::Sam;
use RNAseq_pipeline3 ('get_fh','print_gff');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh');
use RNAseq_GEM3 qw(coords2gtf);
use Tools::Bam qw(bam2coords);

# Declare variables & Get command line options
my $mismatches;
my $file_list;
my $tmpdir;
my $genomedir;

my %options=%{read_config_file()};
$mismatches=$options{'MISMATCHES'};
$file_list=$options{'FILELIST'};
$tmpdir=$options{'LOCALDIR'};
$genomedir=$options{'GENOMEDIR'};

# Get a database connection
my $dbh=get_dbh();

# First ge the files we are going to analyze
my %files=%{read_file_list($file_list)};


# Process the files
foreach my $file (keys %files) {
    print STDERR "Processing $file\n";

    # Get the name for the output file
    my $pair=$files{$file}->[0];
    my $outfileunique=$genomedir.'/'.$pair.'.single.unique.gtf.gz';
    my $outfilemulti=$genomedir.'/'.$pair.'.single.multi.gtf.gz';
    my $outuniquefh=get_fh($outfileunique,1);
    my $outmultifh=get_fh($outfilemulti,1);
    my $sam = Bio::DB::Sam->new(-bam  => $tmpdir.'/'.$file,
				-expand_flags => 1);

    my $all_alignments=$sam->features(-iterator =>1);

    while (my $a=$all_alignments->next_seq()) {
	my $coords=bam2coords($a);
	foreach my $coord (@{$coords}) {
	    my $gtf=coords2gtf($coord,
			       'BAM');
	    if ($coord->{'unique'}) {
		print $outuniquefh $gtf,"\n";
	    } else {
		print $outmultifh $gtf,"\n";
	    }
	}
    }
    close($outuniquefh);
    close($outmultifh);
}

exit;
