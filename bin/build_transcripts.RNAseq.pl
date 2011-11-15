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
# This script should take a gtf/gff file containing the annotation
# and it will extract information for the different transcripts

use RNAseq_pipeline3 ('get_fh','get_log_fh','check_table_existence',
		      'get_annotation_from_gtf','run_system_command');
use RNAseq_pipeline_settings3 qw(read_config_file get_dbh);

my $annotation_file;
my $prefix;
my $table;
my $commondb;
my $outtable;
my $debug=1;

my %options=%{read_config_file()};
$annotation_file=$options{'ANNOTATION'};
$prefix=$options{'PREFIX'};
$table=$options{'TRANSCLASSTABLE'} || die "No TRANSCLASS table found\n";
$commondb=$options{'COMMONDB'};
$outtable=$prefix.'_transcripts.txt';

# Get a log file
my $log_fh=get_log_fh("build_transcripts.RNAseq.log",
		      $debug);
print $log_fh "Extracting all transcripts from $annotation_file\n";

# Check if the table is present in the database and contains data
my $dbh=get_dbh(1);
my $present=check_table_existence($dbh,
				  $table);

my $tablefh=get_fh($outtable,1);
if ($present) {
    # Print the table location for the transcripts of this experiment
    print $log_fh "$table file is present. Skipping\n";
    print $tablefh join("\t",
			$outtable,
			"Present"),"\n";
    close ($tablefh);
    exit;
} else {
    # Continue, as the table must be created
    print $log_fh $table,"\tIs not present\n";
    print $tablefh join("\t",
			$outtable,
			"Generated"),"\n";
    close ($tablefh);
    make_db_table($table)
}

# Check for the presence of the annotation file
unless ($annotation_file &&
	-r $annotation_file) {
    die "Annotation file $annotation_file is not readable\n";
}

# Print out all transcripts
my $transfn=$table.'.txt';
if (-r $transfn) {
    # Check if the file is already present
    print $log_fh $transfn,"\tPresent\n";
} else {
    # Generate all the necessary files
    print $log_fh "Extracting information for each transcript\n";
    print $log_fh "Processing $annotation_file\n";

    # Extract all the transcript coordinates belonging to each of the genes
    my %trans=%{get_annotation_from_gtf($annotation_file,
					$log_fh,
					'transcript')};

    print $log_fh "Printing transcript information in $transfn\n";
    my $transfh=get_fh($transfn,1);
    foreach my $transcript (keys %trans) {
	my $string=build_trans_info($trans{$transcript},
				   $log_fh,
				   'transcript');
	print $transfh $string,"\n";
    }
    close($transfh);
    print $log_fh "done\n";
}

# Load the tables in the database
my $command;
if (-r $table.'.sql') {
    $command="mysql $commondb < $table.sql";
    run_system_command($command,
		       $log_fh);

    $command="mysqlimport -L $commondb $transfn";
    run_system_command($command,
		       $log_fh);

    $command="rm $table.sql $transfn";
    run_system_command($command,
		       $log_fh);
} else {
    die "Can't find $table.sql\n";
}

close($log_fh);

exit;

sub make_db_table {
    my $table=shift;
    my $tablefh=get_fh("$table.sql",1);
    print $tablefh "DROP TABLE IF EXISTS $table;
CREATE TABLE $table (
       gene_id varchar(100) NOT NULL,
       transcript_id varchar(100) NOT NULL,
       type varchar(50) NOT NULL,
       status varchar(20) NOT NULL,
       name varchar(50) NULL,
       description text NULL,
       chr varchar(50) NOT NULL,
       index idx_gene (gene_id),
       index idx_trans (transcript_id),
       index idx_type (type),
       index idx_status (status),
       index idx_name (name),
       index idx_chr (chr)
);\n";
    close($tablefh);
}

sub build_trans_info {
    my $trans=shift;
    my $log_fh=shift;
    my $string;

    my %transcripts;

    if ($trans->{'type'}!~/transcript/o) {
	return();
    }
    my $gene_id=$trans->{'gene_id'} || die "No gene id...\n";
    my $trans_id=$trans->{'transcript_id'} || die "No transcript id...\n";
    my $trans_type=$trans->{'feature'}{'transcript_type'} || 'gene';
    my $trans_status=$trans->{'feature'}{'transcript_status'} || 'unknown';
    my $trans_name=$trans->{'feature'}{'transcript_name'} || '\N';
    my $trans_desc='\N';
    my $trans_chr=$trans->{'chr'} || '\N';

    $string=join("\t",
		 $gene_id,
		 $trans_id,
		 $trans_type,
		 $trans_status,
		 $trans_name,
		 $trans_desc,
		 $trans_chr);

    return($string);
}
