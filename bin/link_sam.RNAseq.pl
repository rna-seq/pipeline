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
# This script should take the files in bam format in the readData directory
# and link them to the correct location in the SAM directory

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub run_system_command);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list');
use Bio::SeqIO;

my $genomefile;
my $prefix;
my $tmpdir;
my $samdir;
my $file_list;
my $bindir;
my $readdir;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$tmpdir=$options{'LOCALDIR'};
$genomefile=$options{'GENOMESEQ'};
$samdir=$options{'SAMDIR'};
$bindir=$options{'BIN'};
$readdir=$options{'READDIR'};

# Add the bindir to the path
$ENV{'PATH'}.=":$bindir";

# Connect to the database
my $dbh=get_dbh();

# First ge the files we are going to analyze
my %files=%{read_file_list()};

# Process the files
foreach my $file (keys %files) {
    my $samfn=$readdir.'/'.$file;
    my $pair=$files{$file}[0];
    my $bamfn=$pair.'.merged';
    print STDERR "Indexing $pair\n";

    # check if the bam file is already present and skip if it is
    if (-r $samdir.'/'.$bamfn.'.bam') {
	print STDERR "$bamfn.bam exists already. Skipping...\n";
	next;
    } else {
	generate_bam_index($samfn,
			   $samdir.'/'.$bamfn);
    }

    print join("\t",
	       $pair,
	       $bamfn.'.bam'),"\n";
}


exit;

sub generate_bam_index {
    my $infile=shift;
    my $bamfile=shift;

    my $command='samtools sort ';
    $command.="$infile $bamfile";
    run_system_command($command);

    $command='samtools index ';
    $command.="$bamfile.bam";
    run_system_command($command);
}
