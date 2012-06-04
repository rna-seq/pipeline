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

# Author : David Gonzalez, david.gonzalez@crg.eu

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
# This script should take a gtf/gff or genpred file containing exons with
# gene names and it will build a set of all possible nonoverlapping exon-exon
# junctions for each of the genes.
# The introns, regardless of the strand will alway be ordered from the 5'
# coordinate to the 3'

use RNAseq_pipeline3 ('get_fh','get_log_fh',
		      'get_annotation_from_gtf',
		      'build_annotated_junctions',
		      'run_system_command',
		      'check_table_existence');
use RNAseq_pipeline_settings3 qw(read_config_file check_db get_dbh);
use Bio::Range;

my $debug=0;

# Get the pipeline options
my %options=%{read_config_file()};
# Set a couple of variables
my $annotation_file=$options{'ANNOTATION'};
my $table=$options{'JUNCTIONSTABLE'};

my $dbh=get_dbh(1);

# Get a log file
my $log_fh=get_log_fh('build_exon_junctions.RNAseq.log',
		      $debug);

# Check if the table is present in the database and contains data
my $present=check_db($table,
		     $options{'PREFIX'}.'_junctions.txt');

if ($present) {
    print $log_fh "$table is present. Skipping\n";
    exit;
} else {
    # Prepare for table generation
    print $log_fh "$table is absent. Generating...\n";
    make_db_table($table)
}

unless ($annotation_file &&
	-r $annotation_file) {
    die "Annotation file $annotation_file is not readable\n";
}
print $log_fh "Building all junctions combinations for each gene\n";
print $log_fh "Processing $annotation_file\n";

# Extract all the exons coordinates belonging to each of the genes
my %genes=%{get_annotation_from_gtf($annotation_file,
				    $log_fh)};


# get all the annotated junctions
my %annotated_junctions=%{build_annotated_junctions(\%genes,
						    $log_fh)};


# Get all the possible junctions    
my %junctions=%{build_all_junctions(\%genes,
				    \%annotated_junctions,
				    $log_fh)};

# Print out all exon junctions
my $juncfn=$options{'JUNCTIONSTABLE'}.'.txt';
if (-r $juncfn) {
    print $log_fh $juncfn,"\tPresent\n";
} else {
    print $log_fh $juncfn,"\tNot Present. building...\n";
    my $junctionfh=get_fh($juncfn,1);
    foreach my $junction (keys %junctions) {
	my ($exon1,$exon2)=split('_splice_',$junction);
	my $end1=(split('_',$exon1))[-2];
	my $start2=(split('_',$exon2))[-3];

	# Build the chromosome id taking into account that the chromosome id may
	# contain '_' which is the character we are using to split
	my @chr=split('_',$exon1);
	splice(@chr,-3,3);
	my $chr=join('_',@chr);
	my $key=$junctions{$junction}[1];
	print $junctionfh join ("\t",
				$junction,
				@{$junctions{$junction}}[0,1],
				$exon1,
				$exon2,
				$chr,
				$end1,
				$start2),"\n";
    }
    close($junctionfh);
}

# Load the tables in the database
my $command;
if (-r $table.'.sql') {
    my $database=$options{'COMMONDB'};
    $command="mysql $database < $table.sql";
    print $log_fh "Executing:\t",$command,"\n";
    print $log_fh `$command`,"\n";
    $command="mysqlimport -L $database $juncfn";
    print $log_fh "Executing:\t",$command,"\n";
    print $log_fh `$command`,"\n";
    sleep(1);
    $command="rm $table.sql $juncfn";
    run_system_command($command,
		       $log_fh);
} else {
    die "Can't find $table.sql\n";
}

# Check if the table has been created correctly
$present=check_table_existence($dbh,
			       $table);

unless ($present) {
    die "ERROR: $table does not seem to be present\n";
}

close($log_fh);

exit;

sub make_db_table {
    my $table=shift;
    my $tablefh=get_fh("$table.sql",1);
    print $tablefh "DROP TABLE IF EXISTS $table;
CREATE TABLE $table (
       junction_id varchar(100) NOT NULL,
       type varchar(20) NOT NULL,
       gene_id varchar(50) NOT NULL,
       exon1_id varchar(50) NOT NULL,
       exon2_id varchar(50) NOT NULL,
       chr varchar(50) NOT NULL,
       start int unsigned NOT NULL,
       end int unsigned NOT NULL,
       index idx_junction (junction_id),
       index idx_gene (gene_id),
       index idx_exon1 (exon1_id),
       index idx_exon2 (exon2_id),
       index idx_chr (chr),
       index idx_start (start),
       index idx_end (end)
);\n";
    close($tablefh);
}

# get the coordinates of each possible exon junction
sub build_all_junctions {
    my $genes=shift;
    my $annotated=shift;
    my $log_fh=shift;

    my %exons;
    my %exon_list;
    my %junctions;
    my $type;
    my $skipped=0;
    
    print $log_fh "Extracting all exons...";
    
    foreach my $gene (keys %{$genes}) {
	foreach my $exon ($genes->{$gene}->{'gene'}->exons()) {
	    my $exon_id=$exon->display_name();
	    
	    $exons{$gene}{$exon_id}=[$exon->strand(),
				     $exon->start(),
				     $exon->end()];
	    $exon_list{$exon_id}{$gene}=1;
	}
    }
    print $log_fh "done\n";

    my $count=keys %exon_list;
    print $log_fh $count,"\tAnnotated exons retrieved\n";

    ### QC
    # Print a list of those exons assigned to more than one gene
    my $exonerror_fh=get_fh('exon.strange.list',1);
    my $strange=0;
    my %strange;
    foreach my $exon_id (keys %exon_list) {
	my $genes=keys %{$exon_list{$exon_id}};
	if ($genes !=1) {
	    foreach my $gene_id (keys %{$exon_list{$exon_id}}) {
		print $exonerror_fh join("\t",
					 $exon_id,
					 $gene_id),"\n";
	    }
	    $strange{$exon_id}=1;
	    $strange++;
	}
    }
    close($exonerror_fh);
    print $log_fh $strange,"\tExons were assigned to more then one gene?\n";

    print $log_fh "Building all gene exon junction combinations\n";
    foreach my $key (keys %exons) {
	my @exons= sort {$exons{$key}->{$a}->[1] <=>
			     $exons{$key}->{$b}->[1]} keys %{$exons{$key}};
	
	while (my $exon1_id=shift(@exons)) {
	    my $exon1_strand=$exons{$key}->{$exon1_id}->[0];
	    foreach my $exon2_id (@exons) {
		# Skip junctions that would be built from exons that overlap
		if (($exons{$key}->{$exon1_id}->[2] >= 
		     $exons{$key}->{$exon2_id}->[1]) ||
		    ($exons{$key}->{$exon1_id}->[1] == 
		     $exons{$key}->{$exon2_id}->[1])){
		    $skipped++;
		    next;
		}

		# Build the junctions
		my $exon2_strand=$exons{$key}->{$exon2_id}->[0];
		my $junction_id=$exon1_id.'_splice_'.$exon2_id;

		### QC
		# Skip exon combinations from different strands
		unless ($exon1_strand eq $exon2_strand) {
		    # This happens in the case of some unprocessed pseudogenes
		    next;
		}

		# Skip exon combinations that include an exon present in more
		# than one gene
		if ($strange{$exon1_id} || $strange{$exon2_id}) {
		    next;
		}

		# Determine if the junction is known or novel
		if (exists $annotated->{$junction_id}){
		    $type='known';
		} else {
		    $type='novel';
		}
		$junctions{$junction_id}->[0]=$type;
		$junctions{$junction_id}->[1]=$key;
	    }
	}
    }
    $count=keys %junctions;
    print $log_fh $count,"\tJunctions retrieved\n";
    print $log_fh $skipped,"\tOverlapping exons skipped\n";
    return(\%junctions);
}
