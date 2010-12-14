#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script should prepare everything for the execution of the RNAseq
# analysis pipeline
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

use Pod::Usage;
use Getopt::Long;
use RNAseq_pipeline_start3 ('defaults','get_species_hash','base_table_build',
			    'check_base_tables','get_existing_data_subs',
			    'get_tables_hash', 'build_table_files',
			    'create_directory_structure','build_file_list',
			    'print_config_file','print_pipeline_file',
			    'add_project','add_experiment');
use RNAseq_pipeline3 ('MySQL_DB_Connect','get_fh','get_log_fh');
use Cwd;

# Declare the variables
# naming variables
my $species;
my $proj_prefix;
my $exp_prefix;

# Cluster specification
my $cluster;

# Behaviour
my $clean;
my $help;
my $man;
my $debug;

# Basic files
my $template;
my $annotation_file;
my $genome_file;
my $file_list; # this file will list all the files in the directory to map as
               # well as if they are paired, and/or stranded.

# Database information
my $host;
my $commondb='RNAseqPipelineCommon';
my $database='RNAseqPipeline';
my $junctionstable;
my $junctionsclasstable;
my $exonsclasstable;

# Mapping settings
my $mismatches;
my $paired;
my $stranded;
my $readlength;
my $mapper;
my $threads;
my $parallel;
my $qualities;

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
my $exp_description;

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
my $localdir=$project_dir."/work";

# Get the options from the command line
my %options=('species' => \$species,
	     'project' => \$proj_prefix,
	     'experiment' => \$exp_prefix,
	     'clean' => \$clean,
	     'help' => \$help,
	     'man' => \$man,
	     'debug' => \$debug,
	     'template' => \$template,
	     'annotation' => \$annotation_file,
	     'genome' => \$genome_file,
	     'files' => \$file_list,
	     'host' => \$host,
	     'commondb' => \$commondb,
	     'database' => \$database,
	     'junctionstable' => \$junctionstable,
	     'junctionsclass' => \$junctionsclasstable,
	     'exonsclass' => \$exonsclasstable,
	     'mismatches' => \$mismatches,
	     'paired' => \$paired,
	     'stranded' => \$stranded,
	     'readlength' => \$readlength,
	     'mapper' => \$mapper,
	     'threads' => \$threads,
	     'parallel' => \$parallel,
	     'localdir' => \$localdir,
	     'genomeindex' => \$genomeindex,
	     'transcriptomeindex' => \$transcriptomeindex,
	     'junctionsindex' => \$junctionsindex,
	     'exonindex' => \$exonindex,
	     'exonsfasta' => \$exonsfasta,
	     'junctionsfasta' => \$junctionsfasta,
	     'transcriptomefasta' => \$transcriptomefasta,
	     'exclusionfile' => \$exclusionfile,
	     'cluster' => \$cluster,
	     'qualities' => \$qualities
    );
GetOptions(\%options,
	   'species|s=s',
	   'project=s',
	   'experiment=s',
	   'clean',
	   'help|h',
	   'man',
	   'debug|d+',
	   'template|t=s',
	   'annotation|a=s',
	   'genome|g=s',
	   'files=s',
	   'host=s',
	   'commondb=s',
	   'database=s',
	   'junctionstable=s',
	   'junctionsclass=s',
	   'exonsclass=s',
	   'mismatches=i',
	   'paired',
	   'stranded',
	   'readlength=i',
	   'mapper=s',
	   'threads=i',
	   'parallel',
	   'genomeindex=s',
	   'transcriptomeindex=s',
	   'junctionsindex=s',
	   'exonindex=s',
	   'exclusionfile=s',
	   'localdir=s',
	   'cluster=s',
	   'qualities=s'
    );

# Print help and exit if required
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Print message if in clean mode
if ($clean) {
    print STDERR "WARNING: Clean mode.\n";
    print STDERR "This will remove all tables for this experiment as well as
the common tables\n";
    print STDERR "Do you want to continue? (y/N)\n";
    my $reply=<STDIN>;
    unless ($reply=~/y/i) {
	die "Aborting\n";
    }
}

# Print message if debuging
my $log_fh=get_log_fh('start_RNAseq_pipeline.log',
		      $debug);
print $log_fh "Executing: $command_line\n";

### TO DO
# check all the options for the presence of dashes or unauthorized characters
# that may cause problems with the database if used for naming tables
check_options(\%options,
	      $log_fh);

# Get some useful subs
*check_default=defaults();

# Check the input
my ($missing,$guessmaster)=check_option_values(\%options);
pod2usage("WARNING: $missing information is necessary but missing") if $missing;

# Check for the existence of the tables that will contain the annotation files
# and the other necessary tables and if they are not present create them
# If clean option is used this step will also delete the tables
check_base_tables($host,
		  $commondb,
		  $log_fh,
		  $clean,);
# Now if we have used calena we have to exit, because the database is empty
if ($clean) {exit;};

# For those values that have not been provided do a bit of guess work
# First get the required subroutines. This MUST be done here, after using check
# values if the default database is to be used
*guess=get_existing_data_subs($genome_file,
			      $annotation_file,
			      $readlength,
			      $mismatches,
			      $species,
			      $project_dir,
			      $commondb,
			      $host,
			      $log_fh);

if (keys %{$guessmaster}) {
    guess_values($guessmaster,
		 \%options);
}

# Check that the template file that has been supplied is readable
pod2usage("I need a readable template file") unless (-r $template);

# All the necessary information is available so we will start the pipeline
print $log_fh "All necessary information seems present.\n";
print $log_fh "Building pipeline for:\n";
print $log_fh join("\n",
		   $species,
		   $proj_prefix.':'.$exp_prefix),"\n";

# Get the species information
my %species=%{get_species_hash($commondb,
			       $host)};

unless ($species{$species}) {
    die "Something has happened and I cannot find $species in $commondb\n";
}
print $log_fh "Species Abbreviation will be: $species{$species}\n";

# Get the identifier string we want:
my $id_string=$proj_prefix.'_'.$exp_prefix;
print $log_fh "Id string: $id_string\n";

# Get the output table fh
my $table_fh=get_fh($id_string.'_start.txt',1);
print $log_fh "Table for this rule is: ${id_string}_start.txt\n";

# Get the necessary tables
my %tables=%{get_tables_hash($id_string)};
print $log_fh "I will try to build the following tables\n";
foreach my $table (keys %tables) {
    print $log_fh $id_string.$table,"\n";
}

# Create the necessary directory structure
my %vars=%{create_directory_structure($id_string,
				      $log_fh,
				      \%options,
				      $table_fh)};

my @table_files=build_table_files($id_string,
				  $vars{'TABLES'});
# print out pipfile
print_pipeline_file(\%tables,
		    \%vars,
		    $template,
		    $log_fh);

# print out the configuration file
print_config_file(\%vars,
		  $log_fh);

# Build the file list
build_file_list(\%vars,
		$log_fh);

# execute the first step of the pipeline
### THE FIRST STEP OF THE PIPELINE should be added here

# Add the project and experiment entries to the respective tables
add_project($host,
	    $commondb,
	    $proj_prefix,
	    $species,
	    $log_fh);

add_experiment(\%options,
	       $log_fh);

# close the log file
close($log_fh);

exit;

sub check_table {
    my $dbh=shift;
    my $table=shift;

    print STDERR "Checking database for $table...\n";

    my ($query,$sth);
    $query ='SELECT count(*) ';
    $query.="FROM $table";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    my $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    my $present=0;
    
    if (@{$results}) {
	# Print the table location for the junctions of this experiment
	print STDERR join("\t",
			  $table,
			  "Present"),"\n";
	$present=1;
    }
    return($present);
}

# This script should for each of the necessary files
sub guess_values {
    my $guessmaster=shift;
    my $options=shift;

    print STDERR "Some values missing; time for some guesswork...";
    foreach my $key (keys %{$guessmaster}) {
	print $log_fh "Guessing $key ...\n";
	my $value=guess($key);
	${$options->{$key}}=$value;
	print $log_fh "done\n";
    }
    print STDERR "done\n";
}


# Check if the options are set and if not set the defaults.
# If any option remains undefiuned it means we have no default and the user
# has to supply it
sub check_option_values {
    my $options=shift;

    my $missing;
    my %guessmaster;
    print STDERR "Checking arguments...\n";

    foreach my $key (sort keys %{$options}) {
	my $value=$options->{$key};
	# Skip the behaviour switches
	if ($key=~/(clean|help|debug|man)/) {
	    next;
	}
	if (check_default($key,
			  $value,
			  \%guessmaster)) {
	    $missing=$key;
	}
    }
    return($missing,\%guessmaster);
}

sub check_options {
    my $options=shift;
    my $log_fh=shift;

    print $log_fh "Checking options for unadvisable characters...\n";
    foreach my $option (keys %{$options}) {
	my $value=${$options->{$option}} || '';
	if ($value=~/([^\w_\/\. ])/) {
	    my $char=$1;
	    die "$value contains an invalid character: '$char'\n";
	}
    }
    print $log_fh "done\n";

    print $log_fh "Print checking options for common problems...\n";
    # Check the species name to see if it has 2 words
    my $species=${$options->{'species'}};
    my @species=split(/\s+/,$species);
    if (@species < 2) {
	die "Species name ($species) does not look right. Should have a Genus and species name at least\n";
    }

    print $log_fh "done\n";
}

__END__

=head1 NAME
    
    start_RNAseq_pipeline.3.0.pl
    
=head1 SYNOPSIS
    
    start_pair_pipeline.3.0.pl -species ... -genome ... -annotation ... -project ... -experiment ... -template ... -readlength...
    
  Help options:
    -help:           brief help message
    -man:            full documentation
    -debug:
  Mandatory options:
    -species:        Species for which the pipeline is run.
    -genome:         File with the genomic sequence.
    -annotation:     File with the annotation to use.
    -project:        The project to which the experiment will be added.
    -experiment:     The set of reads to be added.
    -template:       File containing the commands that will be executed.
    -readlength:     Nucleotide length of the reads.
  Mapping Options:
    -mapper:         Mapper to be used.
                      Defaults to GEM which is the only one supported currently
    -mismatches:     Number of mismatches with which the mapping will be done.
                      Default 2.
    -paired:         Reads are paired.
    -stranded:       Reads are stranded.
    -threads:        Number of threads to use for mapping.
                      (if mapper allows multiple threads)
    
  Advanced Options:
    -database:       Sets the database to use for the experiment tables.
    -commondb:       Sets the database to use for the common tables.
    -host:           Sets the host where the databases are located.
    -localdir:       Directory in which to store the temporary files generated during the pipeline execution. It is recommandable to set it to a local drive.

        
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
