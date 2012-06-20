#!/soft/bin/perl

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
# This script will take the output of cufflinks and create a 
# table to load in the database with this information

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_log_fh parse_gff_line run_system_command);
use RNAseq_pipeline_settings3 qw(read_file_list read_config_file);

# declare some variables
my $debug=0;
my $samdir;
my $mysqldir;
my $table;

my %options=%{read_config_file()};
$samdir=$options{'SAMDIR'};
$mysqldir=$options{'TABLES'};
$table=$options{'PREFIX'}.'_novel_transcripts';

# Get the log_fh
my $log_fh=get_log_fh('build_novel_transcripts.RNAseq.log',
		      $debug);

# First get a list of the bed files we are going to process
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# read the results from the files named after the lane ids
foreach my $lane (keys %lanes) {
    my $filename=$samdir.'/'.$lane.'.transcripts.gtf';

    if (-r $filename) {
	process_file($filename,
		     $lane,
		     $log_fh);
    } else {
	die "I can't find $filename\n";
    }
}

close($log_fh);

make_db_table($table,
	      $mysqldir);

exit;

sub make_db_table {
    my $table=shift;
    my $mysqldir=shift;
    my $tablefh=get_fh("$table.sql",1);
    print $tablefh "DROP TABLE IF EXISTS $table;
CREATE TABLE $table (
       cuff_id varchar(50) NOT NULL,
       trans_id varchar(50) NOT NULL,
       FPKM int unsigned NOT NULL,
       Sample varchar(50) NOT NULL,
       index idx_Sample (Sample),
       index idx_cuff_id (cuff_id),
       index idx_trans (trans_id)
);\n";
    close($tablefh);

    my $command="mv $table.sql $mysqldir";
    run_system_command($command);
}

sub process_file {
    my $infile=shift;
    my $lane=shift;
    my $log_fh=shift;

    my $infh=get_fh($infile);

    while (my $line=<$infh>) {
	my %line=%{parse_gff_line($line,
				  $log_fh)};

	my $type=$line{'type'};

	if ($type=~/transcript/) {
	    my $cuff_id=$line{'feature'}{'gene_id'};
	    my $trans_id=$line{'feature'}{'transcript_id'};
	    my $fpkm=$line{'feature'}{'FPKM'};
	    if ($fpkm > 0) {
		print join("\t",
			   $cuff_id,
			   $trans_id,
			   $fpkm,
			   $lane),"\n";
	    }
	}
    }
    close($infh);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	push @{$lanes{$files->{$file}->[0]}},$files->{$file}->[1];
    }

    return(\%lanes);
}
