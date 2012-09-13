#!/usr/bin/perl

use strict;
use warnings;

#    GRAPE
#    Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
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

# Objective
# This script should prepare everything for the execution of the RNAseq
# analysis pipeline
# It will also do some basic checks on the provided input.
# - This is check the format of the gft annotation file if ok add to the DB
# - Check the chromosome entries in the genome file
# - Insert the genome annotation into the 
# It can be used also to provide information to the database such as the 
# - Takes one species and project identification as arguments and creates a
# code to identify them (DmeImdisc)
# - Checks for the presence of all the necessary directories and if they
# are not present they are created
# -Check for the presence of all the necessary scripts (not implemented yet)
# -Check for the presence of all the necessary sql table creation files and
# if necessary create them
# -Check for the presence of a paired file if paired end reads are provided
# (not implemented yet)

### TO DO:
# Add an option to insert an experiment description or a project description
# Add umask(0115)

# Add the path to the current library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

use Pod::Usage;
use GRAPE::Start ('base_table_build',
		  'get_existing_data_subs',
		  'get_tables_hash', 'build_table_files',
		  'create_directory_structure','build_file_list',
		  'print_config_file','print_pipeline_file',
		  'add_project','add_proj_info','add_exp_info',
		  'add_experiment',
		  'clear_experiment');

use RNAseq_pipeline3 ('MySQL_DB_Connect','get_fh','get_log_fh');
use GRAPE::Logs;
use GRAPE::Options;
use Cwd;

my $tables=base_table_build();
my $options=GRAPE::Options->new($tables);

# Behaviour
my $debug=$options->get_debug();

# Print help and exit if required
pod2usage(1) if $options->get_help();
pod2usage(-verbose => 2) if $options->get_man();

if ($debug) {
    $options->print_options(1);
}
# Declare the variables
# naming variables
my $species=$options->get_species() || die "No species\n";
my $project=$options->get_project() || die "No project\n";
my $experiment=$options->get_experiment() || die "No experiment\n";

# Database information
# Most of the connection information will be read from the .my.cnf file, so it
# is not necessary to supply this
my $host=$options->get_host() || die "Unknown host\n";
my $commondb=$options->get_commondb() || die "Unknown database Common\n";
my $database=$options->get_database() || die "Unknown database\n";

if ($options->get_clean()) {
    clear_experiment($options);
    exit;
}

# Cluster specification
my $cluster=$options->get_cluster() || die "No cluster\n";

# Basic files
my $template=$options->get_template() || die "No template\n";
my $annotation_file=$options->get_annotation_file();
my $genome_file=$options->get_genome_file();
my $file_list='read.list'; # this file will list all the files in the 
                           # directory to map as
                           # well as if they are paired, and/or stranded.

# Common tables information
my $junctionstable=$options->get_junctionstable() || die "Unknown junctions table\n";
my $geneclasstable=$options->get_geneclass() || die "Unknown geneclass table";
my $transclasstable=$options->get_transclass() || die "Unknown transclass table";
my $junctionsclasstable=$options->get_junctionsclass() || die "Unknown junctionsclass table";
my $exonsclasstable=$options->get_exonsclass() || die "Unknown exonsclass table";

# Mapping settings
my $mismatches=$options->get_mismatches() || die "Unknown mismatches\n";
my $paired=$options->get_paired() || 0;
my $stranded=$options->get_stranded() || 0;
my $readlength=$options->get_readlength() || warn "No read length supplied\n";
my $mapper=$options->get_mapper() || warn "Unknown mapper\n";
my $threads=$options->get_threads() || die "Unknown threads\n";
my $qualities=$options->get_qualities() || die "Unknown qualities\n";

# Flux settings
my $fluxmem=$options->get_fluxmem() || warn "No fluxmem set\n";

# Overlap settings

# Indices
my $genomeindex;
my $transcriptomeindex;
my $junctionsindex;
my $exonindex;
my $exonsfasta;
my $junctionsfasta;
my $transcriptomefasta;
my $exclusionfile;

# Meta data
my $cellline;
my $compartment;
my $proj_description;
my $run_description;
my $rnafraction;
my $bioreplicate;

# Preprocessing
my $preprocess=$options->get_preprocess();
my $preprocess_trim_length=$options->get_preprocess_trim_length();

# Get the command line for the log file
my $program_name=$0;
$program_name=~s/.*\///;
my $command_line=join(' ',
		      $program_name,
		      @ARGV);

# Get the project directory where the program will be run
my $project_dir=getcwd();
chomp($project_dir);

# Set the temporary directory
my $localdir=$options->get_localdir() || die "Unable to set localdir\n";

# Get the options from the options object
my %options=(
	     'genome' => \$genome_file,
	     'files' => \$file_list,
	     'host' => \$host,
	     'commondb' => \$commondb,
	     'database' => \$database,
	     'geneclass' => \$geneclasstable,
	     'transclass' => \$transclasstable,
	     'junctionstable' => \$junctionstable,
	     'junctionsclass' => \$junctionsclasstable,
	     'exonsclass' => \$exonsclasstable,
	     'mismatches' => \$mismatches,
	     'paired' => \$paired,
	     'stranded' => \$stranded,
	     'readlength' => \$readlength,
	     'mapper' => \$mapper,
	     'threads' => \$threads,
	     'localdir' => \$localdir,
	     'genomeindex' => \$genomeindex,
	     'transcriptomeindex' => \$transcriptomeindex,
	     'junctionsindex' => \$junctionsindex,
	     'exonsfasta' => \$exonsfasta,
	     'junctionsfasta' => \$junctionsfasta,
	     'transcriptomefasta' => \$transcriptomefasta,
	     'exclusionfile' => \$exclusionfile,
	     'cluster' => \$cluster,
	     'qualities' => \$qualities,
	     'cellline' => \$cellline,
	     'compartment' => \$compartment,
	     'run_description' => \$run_description,
	     'projdesc' => \$proj_description,
	     'rnafrac' => \$rnafraction,
	     'bioreplicate' => \$bioreplicate,
	     'preprocess' => \$preprocess,
	     'preprocess_trim_length' => \$preprocess_trim_length,
	     'fluxmem' => \$fluxmem
    );

# Additional options should be files to be run
my @files=sort @{$options->get_filelist()};
my $quickrun=0;
if (@files) {
    # Check that the reads are fastq
    foreach my $file (@files) {
	unless ($file=~/fastq(.gz)?$/) {
	    die "Currently quickrun is available only for fastq files\n";
	}
    }
    $quickrun=1;
}

# If this is a quickrun we need to link the files in the readData
# directory and the read.list.txt file
if ($quickrun) {
    create_read_list_file($options);
    link_files($options);
}

# Get a log file
my $log_fh=GRAPE::Logs->new('start_RNAseq_pipeline.log',
			    $debug);
$log_fh->printlog("Executing: $command_line");


# All the necessary information is available so we will start the pipeline
$log_fh->printlog("All necessary information seems present.");
$log_fh->printlog("Building pipeline for:");
$log_fh->printlog(join("\n",
		   $species,
		   $project.':'.$experiment));

# Get the identifier string we want:
my $prefix=$options->get_prefix();
my $id_string=$prefix;
$log_fh->printlog("Id string: $prefix");

# Get the output table fh
my $table_fh=get_fh($prefix.'_start.txt',1);
$log_fh->printlog("Table for this rule is: ${prefix}_start.txt");

# Get the necessary tables
my %tables=%{get_tables_hash($prefix)};
$log_fh->printlog("I will try to build the following tables");
foreach my $table (keys %tables) {
    $log_fh->printlog($prefix.$table);
}

# Create the necessary directory structure
my %vars=%{create_directory_structure($log_fh,
				      $options,
				      $table_fh)};

my @table_files=build_table_files($options,
				  $vars{'TABLES'});
# print out pipfile
print_pipeline_file(\%tables,
		    \%vars,
		    $template,
		    $log_fh);

# Print out the configuration file
print_config_file(\%vars,
		  $log_fh);

# close the log file
#close($log_fh);

exit;

sub create_read_list_file {
    my $self=shift;

    my @files=sort @{$options->get_filelist()};
    my $paired=$self->get_paired();

    my $fh=get_fh('read.list.txt',1);


    foreach my $file (@files) {
	$file=~s/.*\///;
	$file=~s/.gz$//;
    }

    if ($paired) {
	print $fh join("\t",
		       $files[0],
		       'Quick',
		       'Quick.1',
		       'Quick'),"\n";
	print $fh join("\t",
		       $files[1],
		       'Quick',
		       'Quick.2',
		       'Quick'),"\n";
    } else {
	print $fh join("\t",
		       $files[0],
		       'Quick',
		       'Quick',
		       'Quick'),"\n";
    }
    close($fh);
}

sub link_files {
    my $self=shift;

    my @files=sort @{$options->get_filelist()};

    foreach my $file (@files) {
	my $source=$file;
	$file=~s/.*\///;
	my $dir.='readData/'.$file;
	my $symlink_ok = eval {symlink($source,$dir);1};
	die "Unable to create link\n" unless $symlink_ok;
    }
}

__END__

=head1 NAME
    
    start_RNAseq_pipeline.pl
    
=head1 SYNOPSIS
    
    start_pair_pipeline.3.0.pl -species ... -genome ... -annotation ... -readlength...
    
  Help options:
    -help:           brief help message
    -man:            full documentation
    -debug:

  Behavior options
    -clean:          Remove all the tables corresponding to the project and
                     experiment as well as all directories

  Mandatory options:
    -species:        Species for which the pipeline is run. Mandatory
    -genome:         File with the genomic sequence. Mandatory
    -annotation:     File with the annotation to use. Mandatory
    -project:        The project to which the experiment will be added.
                     Recommended
    -experiment:     The set of reads to be added. Recommended
    -template:       File containing the commands that will be executed.
                     Mandatory unless present in the same directory as
                     template.txt
    -readlength:     Nucleotide length of the reads. Mandatory
    -qualities:      Encoding of the qualities in fastq format (solexa|phred|
                     ignore). The none option will perform the mapping ignoring
                     the quality information

  Mapping Options:
    -mapper:         Mapper to be used.
                     Defaults to GEM which is the only one supported currently
    -mismatches:     Number of mismatches with which the mapping will be done.
                     Default 2.
    -stranded:       Reads are stranded. Default unstranded
    -threads:        Number of threads to use for those parts of the pipeline
                     that can benefit from multithreading. Default 2
    
  Advanced Options:
    -database:       Sets the database to use for the experiment tables.
                     Recommended
    -commondb:       Sets the database to use for the common tables.
                     Recommended
    -localdir:       Directory in which to store the temporary files generated
                     during the pipeline execution. It is advisable to set it
                     to a local drive, particularly in NFS mounted systems where
                     it could cause significant network traffic if not. The
                     default is a directory named work within the project
                     directory

  Optional
    -cellline:      Sets the cell line on which the experiment was performed.
    -compartment:   Sets the compartment on which the experiment was performed.
    -run_description:   Experiment description
    -projdesc:      Project description
    -rnafrac:       RNA fraction on which the experiment was performed.
    -bioreplicate:  Bioreplicate (if the experiment is a bioreplicate)
    -preprocess:    Preprocessing script to be run on each of the read files
                    before anything else
    -preprocess_trim_length: Length by which the reads are trimmed. Must be an
                             integer N or in the form <=N where N is the maximum
                             number of nucleotides that could be trimmed (as
                             when an adaptor is trimmed and the whole adaptor
                             is found).
    -fluxmem:       Amount of memory to use for the Flux Capacitor during
                    transcript quantification. Default 16G
    -genomeassembly:    Assembly version of the reference genome.
    -genomesource:  Source of the reference genome
    -gender:        Gender of the individual from which the RNA was obtained.
    -annotationversion: Version of the reference annotation
    -annotationsource:  Source of the reference annotation
    -protocolinfo:  Information on the protocol used
    -protocolshort: Abbreviation for the protocol used
        
=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.

=item B<-clean>
    
    This option will remove all the tables from the database as well as removing
    all the directories in the project directory with the exception of the bin,
    results, sequences and readData directories.
    
=item B<-preprocess>

    The script or command line supplied here must take as an input a fasta/fastq
    file and output fasta/fastq format

=back
    
=head1 DESCRIPTION
    
    This script creates the configuration files and directory structure that are
    necessary for running the GRAPE RNAseq analysis pipeline. Ideally each run
    should be configured via buildout using an accessions and a profile files
    that will be kept for future reference in order to store the complete
    settings used.

    There is also a minimal QuickRun mode, which requires the species, genome,
    annotation and read length to be provided on the command line as well as
    one (single) or two (paired) read files in fastq(.gz) format. In this case
    all additional options are set to reasonable defaults and the database used
    is the test database that is usually available to all users in most standard
    MySQL installations. For these runs the template file if not provided as an
    argument should be in the project directory and named template.txt. The
    qualities are not taken into account in the mapping.

=cut
