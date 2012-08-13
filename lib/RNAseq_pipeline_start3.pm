package RNAseq_pipeline_start3;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export subroutines to the caller namespace
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('check_option_characters','defaults','check_base_tables',
	    'get_species_hash','base_table_build',
	    'get_existing_data_subs');
push @EXPORT_OK,('get_tables_hash','build_table_files');
push @EXPORT_OK,('create_directory_structure');
push @EXPORT_OK,('print_config_file','print_pipeline_file','build_file_list');
push @EXPORT_OK,('add_project','add_experiment','add_proj_info','add_exp_info');
push @EXPORT_OK,('clear_tables','clear_dirs','clear_common_tables',
		 'clear_files');

# Set strict and warnings
use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# Load other modules required
## POSIX 
use POSIX qw(uname);
use Cwd;
use Cwd 'abs_path';

## Batabase interface
use DBI;
use DBD::mysql;

## Bio::Perl modules
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;

## Pipeline specific modules
use RNAseq_pipeline3 ('get_fh','MySQL_DB_Connect','check_file_existence',
		      'check_gff_file','check_fasta_file','run_system_command',
		      'check_table_existence','check_field_value');

### Subroutines

### Subroutines for checking the options as well as what is present already in
# the database
# Check the options that are provided to the script in order to decide if there 
# are any unadvisable characters
## OK
sub check_option_characters {
    my $options=shift;
    my $log_fh=shift;
    my $problems='';

    print $log_fh "Checking options for unadvisable characters...\n";
    foreach my $option (keys %{$options}) {
	my $value=${$options->{$option}} || '';
	# Experiment id should only contain alphanumeric and underscore
	if ($option eq 'experiment') {
	    if ($value=~/([^\w_])/o) {
		my $char=$1;
		$problems.="WARNING: Value $value corresponding to $option contains an invalid character: '$char'\n";
	    }
	    next;
	} elsif ($option eq 'project') {
	    # Project id should only be alphanumeric
	    if ($value=~/([^\w])/o) {
		my $char=$1;
		$problems.="WARNING: Value $value corresponding to $option contains an invalid character: '$char'\n";
	    }
	    next;
	} elsif ($option eq 'qualities') {
	    # Qualities can only be solexa, phred or ignore
	    unless ($value=~/^(solexa|phred|ignore)$/o) {
		my $char=$1;
		$problems.="WARNING: Value for qualities is set to $value, but the only valid options are 'solexa', 'phred' or 'ignore'\n";
	    }
	    next;
	} elsif ($option eq 'preprocess_trim_length') {
	    if ($value) {
		unless ($value=~/^(<?=?(\d)+)$/o) {
		    my $char=$1;
		    $problems.="WARNING: Value for preprocess_trim_length is set to $value, which is not a valid option\n";
		}
	    }
	    next;
	} elsif ($value=~/([^\w_\/\. -])/o) {
	    my $char=$1;
	    $problems.="WARNING: Value $value corresponding to $option contains an invalid character: '$char'\n";
	} else {
	    print $log_fh "Option $option => $value OK\n";
	}
    }
    print $log_fh "done\n";

    print $log_fh "Print checking options for common problems...\n";
    # Check the species name to see if it has 2 words
    my $species=${$options->{'species'}};
    my @species=split(/\s+/,$species);
    if (@species < 2) {
	$problems.="WARNING: Species name ($species) does not look right. Should have a Genus and species name at least\n";
    }
    print $log_fh "done\n";

    return($problems);
}

# This contains the default values for each of the variables and contains a
# subroutine that will set them or change them
## OK
sub defaults {
    my $log_fh=shift;
    my %defaults=('species' => undef,
		  'project' => undef,
		  'experiment' => undef,
		  'template' => 'template.txt',
		  'annotation' => undef,
		  'genome' => undef,
		  'files' => 'read.list.txt',
		  'commondb' => 'RNAseqPipelineCommon',
		  'database' => 'RNAseqPipeline',
		  'mismatches' => 2,
		  'paired' => 0,
		  'stranded' => 0,
		  'readlength' => undef,
		  'mapper' => 'GEM',
		  'threads' => 2,
		  'localdir' => undef,
		  'cluster' => '-', # Should be defined to the local cluster
		  'qualities' => undef,
		  'preprocess' => 'zcat',
		  'preprocess_trim_length' => 0,
		  'fluxmem' => '16G'
	);
    my $check_default=sub {
	my $variable=shift;
	my $value=shift;
	my $guessmaster=shift;

	my $missing=0;
	print STDERR $variable,":\t";

	# Check if a value has been supplied
	if (defined $$value) {
	    print STDERR "Value supplied as: $$value\n";

	    # If the value is one of the genome or annotation files check that
	    # it actually exists
	    if (($variable eq 'annotation') ||
		($variable eq 'genome')){
		check_file_existence($$value,
				     $log_fh);
	    }
	} elsif (exists $defaults{$variable}) {
	    print STDERR "No value supplied. ";
	    if (defined $defaults{$variable}) {
		my $default=$defaults{$variable};
		$$value=$default;
		print STDERR "Setting to $default\n";
	    } else {
		print STDERR "WARNING: No default for $variable please specify a value\n";
		$missing=1;
	    }
	} else {
	    print STDERR "No value supplied. ";
	    print STDERR "I'll try to guess...\n";
	    $guessmaster->{$variable}=1;
	}
	return($missing);
    };
    return($check_default);
}

## Subroutines used for building the tables for the common database and
# populating it
# Build the tables that will be shared between the projects and experiments
sub check_base_tables {
    my $database=shift;
    my $log_fh=shift;
    my $clean=shift;
    
    my $dbh=MySQL_DB_Connect($database);
    
    my %tables=%{base_table_build()};
    
    print $log_fh "Checking the common tables...\n";
    
    foreach my $table (keys %tables) {
	my $present=check_table_existence($dbh,
					  $table);
	
	if ($present) {
	    # Print the table location for the junctions of this experiment
	    print $log_fh join("\t",
			      $table,
			      "Present"),"\n";
                print $log_fh $tables{$table},"\n";
	    if ($clean) {
		print STDERR "Deleting $table\n";
		my $query="drop table $table";
		my $sth=$dbh->prepare($query);
		print STDERR "Executing: $query\n";
		$sth->execute();
	    }
	} else {
	    unless ($clean) {
		# Continue, as the table must be created
		print $log_fh $table,"\tIs not present\n";
		my $file_name=$table.'.sql';
		my $table_fh=get_fh($file_name,1);
		print $table_fh "SET FOREIGN_KEY_CHECKS=0;\n";
		print $table_fh $tables{$table},"\n";
		print $table_fh "SET FOREIGN_KEY_CHECKS=1;\n";
		close($table_fh);
		my $command="mysql $database < $file_name";
		run_system_command($command,
				   $log_fh);
		$command="rm $file_name";
		run_system_command($command);
		print $log_fh join("\t",
				   $table,
				   "Generated"),"\n";
	    }
	}
    }
}

# Actual table syntax
sub base_table_build {
    my %tables=('projects' => 'CREATE TABLE IF NOT EXISTS projects (
       project_id varchar(50) NOT NULL,
       species varchar(50) NOT NULL,
       proj_description mediumtext,
       secret varchar(50) NULL,
       PRIMARY KEY (project_id)
);',
		'experiments' => 'CREATE TABLE IF NOT EXISTS experiments (
       experiment_id varchar(50) NOT NULL,
       project_id varchar(50) NOT NULL REFERENCES projects(project_id),
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       template_file varchar(200) NOT NULL,
       read_length smallint unsigned NULL,
       mismatches smallint unsigned NOT NULL,
       exp_description mediumtext,
       expDate date NULL,
       CellType varchar(50) NULL,
       RNAType varchar(50) NULL,
       Compartment varchar(50) NOT NULL DEFAULT "CELL",
       Bioreplicate varchar(10) DEFAULT 1,
       partition varchar(50) NULL,
       md5sum varchar(45) NULL,
       Preprocessing varchar(45) NULL,
       annotation_version varchar(45) NULL,
       genome_assembly varchar(45) NULL,
       paired bit(1) DEFAULT 1, 
       lab varchar(45) NULL,
       protocol_id mediumint unsigned DEFAULT 1,
       PRIMARY KEY (project_id,experiment_id)
);',
		'annotation_files' => 'CREATE TABLE IF NOT EXISTS annotation_files (
       annotation_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       annotation varchar(60) NOT NULL,
       location varchar(200) NOT NULL,
       version varchar(50) NULL,
       source mediumtext NULL,
       index idx_annotation (annotation),
       PRIMARY KEY (annotation_id)
) ENGINE=INNODB;',
		'genome_files' => 'CREATE TABLE IF NOT EXISTS genome_files (
       genome_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       genome varchar(60) NOT NULL,
       location varchar(200) NOT NULL,
       assembly varchar(50) NULL,
       source mediumtext NULL,
       gender varchar(20) NULL,
       index idx_genome (genome),
       PRIMARY KEY (genome_id)
);',
		'indices' => 'CREATE TABLE IF NOT EXISTS indices (
       index_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       type varchar(100) NOT NULL,
       location varchar(200) NOT NULL,
       PRIMARY KEY (index_id)
) ENGINE=INNODB;',	
		'exclusion_files' => 'CREATE TABLE IF NOT EXISTS exclusion_files (
       exclusion_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       location varchar(200) NOT NULL,
       PRIMARY KEY (exclusion_id)
);',
		'annotation_tables' => 'CREATE TABLE IF NOT EXISTS annotation_tables (
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       type varchar(20) NOT NULL,
       table_name varchar(100) NOT NULL,
       index idx_genome (genome_id),
       index idx_annotation (annotation_id),
       index idx_type (type)
) ENGINE=INNODB;',
		'fasta_files' => 'CREATE TABLE IF NOT EXISTS fasta_files (
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       filetype varchar(20) NOT NULL,
       file_name varchar(200) NOT NULL,
       index idx_genome (genome_id),
       index idx_annotation (annotation_id),
       index idx_type (filetype),
       PRIMARY KEY (genome_id,annotation_id,filetype)
);',
		'species_info' => 'CREATE TABLE IF NOT EXISTS species_info (
       species_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species varchar(200) NOT NULL,
       genus varchar(100) NOT NULL,
       sp_alias varchar(100) NOT NULL,
       abbreviation varchar(20) NOT NULL,
       PRIMARY KEY (species_id),
       index idx_alias (sp_alias)
);',
		'protocol_info' => 'CREATE TABLE IF NOT EXISTS protocol_info (
       protocol_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       abbreviation varchar(20) NOT NULL,
       description varchar(100) NOT NULL,
       PRIMARY KEY (protocol_id)
);'
);
    return(\%tables);
}

## Species information. This should be added to a table with the same
# subroutine as an interface, as that woudl allow us to add a species
sub get_species_hash {
    my $database=shift;
    my $dbh=MySQL_DB_Connect($database);
    my %species;
    my $table='species_info';

    my ($query,$sth);
    $query ='SELECT species, abbreviation ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($species,$abbrev)=$sth->fetchrow_array()) {
	if ($species{$species}) {
	    print STDERR $species,"Seems to appear more than once in $table\n";
	} else {
	    $species{$species}=$abbrev;
	}
    }
    $dbh->disconnect();

    return(\%species);
}

sub gset_field_id {
    my $dbh=shift;
    my $value=shift;
    my $table=shift; # species_info
    my $field=shift; # species_id
    my $type=shift; # species
    my $species=shift;
    my $log_fh=shift;

    my $field_id;
    my $value2=$value;
    $value2=~s/.*\///;
    my $sth=check_field_value($dbh,
			      [$type,$value2],
			      $table,
			      $field);

    my $count=$sth->execute($value2);

    print $log_fh "Found $count entries in $table for $value\n";
    if ($count == 1) {
	# The entry is present
	($field_id)=$sth->fetchrow_array();
    } elsif ($count > 1) {
	# There is a problem
	die "Entry $value is present more than once in $table\n";
    } else {
	if ($type eq 'species') {
	    # The entry is absent and we must set it
	    # Get the genus and abbreviation
	    my ($genus,$specific)=split(/\s+/,$value,2);
	    my $abbreviation=join('',
				  substr($genus,0,1),
				  substr($specific,0,3));

	    # Insert the info into the database
	    my $query;
	    $query ="INSERT INTO $table ";
	    $query.='SET species = ? , genus = ? , abbreviation = ?, sp_alias = ? ';
	    print $log_fh "Executing: $query\n";
	    my $sth2=$dbh->prepare($query);
	    $sth2->execute($value,$genus,$abbreviation,'-');
	} elsif ($type eq 'genome') {
	    # First check if the fasta files is ok
	    my $fileok=check_fasta_file($value);
	    my ($version,$source,$gender)=('\N','\N','\N');

	    unless ($fileok) {
		die "There is a problem with the genome file, so I'm quitting. Please check it\n";
	    }

	    # The entry is absent and we must set it
	    my $query;
	    $query ="INSERT INTO $table ";
	    $query.='SET genome= ? , species_id= ? , location= ? ';
	    print $log_fh "Executing: $query\n";
	    my $sth2=$dbh->prepare($query);
	    $sth2->execute($value2,$species,$value);
	} elsif ($type eq 'annotation') {
	    # first we will check the file to see if it conforms to the gff
	    # format
	    my $fileok=check_gff_file($value);
	    my ($version,$source)=('\N','\N');

	    # The entry is absent and we must set it
	    my $query;
	    $query ="INSERT INTO $table ";
	    $query.='SET annotation= ? , species_id= ? , location= ? ';
	    my $sth2=$dbh->prepare($query);
	    $sth2->execute($value2,$species,$value);
	} else {
	    die "Unknown field $value\n";
	}
	# Get the Id of the value we just inserted
	$sth->execute($value2);
	($field_id)=$sth->fetchrow_array();
    }
    print $log_fh "$type: $field_id obtained from $table for $value\n";
    return($field_id);
}

# This subroutine will check the tables for available files using the
# following information:
# genome, annotation, project, experiment, mismatches and readlength
### TO DO rewrite to reduce code, some is very redundant
sub get_existing_data_subs {
    my $genome=shift;
    my $annotation=shift;
    my $readlength=shift;
    my $mismatches=shift;
    my $species=shift;
    my $project_dir=shift;
    my $database=shift;
    my $log_fh=shift;
    my %guessing_subs;

    # Connect to the database as we will need the info in
    my $dbh=MySQL_DB_Connect($database);

    # Get the genome unique id and the annotation unique id
    my $species_id=gset_field_id($dbh,
				 $species,
				 'species_info',
				 'species_id',
				 'species',
				 '',
				 $log_fh);
    my $genome_id=gset_field_id($dbh,
				$genome,
				'genome_files',
				'genome_id',
				'genome',
				$species_id,
				$log_fh);

    my $annotation_id=gset_field_id($dbh,
				    $annotation,
				    'annotation_files',
				    'annotation_id',
				    'annotation',
				    $species_id,
				    $log_fh);

    # Get the tables where the files and tables dependent on the genome
    # and annotation will be stored
    my $table1='indices';
    my $table2='exclusion_files';
    my $table3='annotation_tables';
    my $table5='fasta_files';

    $genome=~s/.*\///;
    $annotation=~s/.*\///;
    $genome=~s/\/$//;
    $annotation=~s/\/$//;

    # Prepare the queries
    # Index query: depends on genome and annotation
    my ($query1,$sth1);
    $query1 ='SELECT location ';
    $query1.="FROM $table1 ";
    $query1.='WHERE genome_id = ? AND annotation_id = ? AND type = ?';
    $sth1=$dbh->prepare($query1);

    # Exclusion query: depends only on annotation
    my $sth2=check_field_value($dbh,
			       ['annotation_id',$annotation],
			       $table2,
			       'location');

    # Annotation tables: Depend only on annotation
    my ($query3,$sth3);
    $query3 ='SELECT table_name ';
    $query3.="FROM $table3 ";
    $query3.='WHERE annotation_id = ? AND type = ?';
    $sth3=$dbh->prepare($query3);

    # Fasta files for building indices: Depend on annotation and genome
    my ($query5,$sth5);
    $query5 ='SELECT file_name ';
    $query5.="FROM $table5 ";
    $query5.='WHERE annotation_id = ? AND genome_id = ? AND filetype = ?';
    $sth5=$dbh->prepare($query5);

    # Get the indices
    # Prepare an sth for the indices table
    my ($ins_query1,$ins_sth1);
    $ins_query1 ="INSERT INTO $table1 ";
    $ins_query1.='SET genome_id = ?,annotation_id = ?,species_id = ?,';
    $ins_query1.='type = ?,location = ?';
    $ins_sth1=$dbh->prepare($ins_query1);

    $guessing_subs{'genomeindex'}=sub {
	my $filetype=shift;
	my $present=0;
	my $fasta_name=$genome;
	$fasta_name=~s/\.fa$//;
	my $location=$project_dir.'/GEMIndices/'.$fasta_name;
	
	my $count=$sth1->execute($genome_id,0,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth1->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query1 for $location\n";
	    $ins_sth1->execute($genome_id,0,$species_id,
			       $filetype,$location);
	}
	return($location);
    };
    $guessing_subs{'transcriptomeindex'}=sub {
	my $filetype=shift;
	my $present=0;
	my $fasta_name=$annotation.'.'.$genome;
	$fasta_name=~s/\.g[tf]f//g;
	$fasta_name=~s/\.fa(stq)?$//;
	$fasta_name.='.transcriptome';
	my $location=$project_dir.'/GEMIndices/'.$fasta_name;

	my $count=$sth1->execute($genome_id,$annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth1->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome $annotation\n";
		$present=1
	    }
	} elsif (defined $count) {
	    print $log_fh "$filetype absent\n";
	} else {
	    print $log_fh "Failed at:\t",join("\t",$genome_id,$annotation_id,$filetype),"\n";
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query1 for $location\n";
	    $ins_sth1->execute($genome_id,$annotation_id,$species_id,
			       $filetype,$location);
	}
	return($location);
    };
    $guessing_subs{'junctionsindex'}=sub {
	my $filetype=shift;
	my $present=0;
	my $fasta_name=$annotation.'.'.$genome;
	$fasta_name=~s/\.g[tf]f//g;
	$fasta_name=~s/\.fa(stq)?$//;
	$fasta_name.='.junctions';
	my $location=$project_dir.'/GEMIndices/'.$fasta_name;

	my $count=$sth1->execute($genome_id,$annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth1->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome $annotation\n";
		$present=1
	    }
	} elsif (defined $count) {
	    print $log_fh "$filetype absent\n";
	} else {
	    print $log_fh "Failed at:\t",join("\t",$genome_id,$annotation_id,$filetype),"\n";
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query1\n";
	    $ins_sth1->execute($genome_id,$annotation_id,$species_id,
			       $filetype,$location);
	}
	return($location);
    };

    # Get the exclusion file if present or set the location to build it if
    # absent
    $guessing_subs{'exclusionfile'}=sub {
	my $filetype=shift;
	my $present=0;
	my $exclusion_name=$annotation;
	$exclusion_name=~s/\.g[tf]f//;
	$exclusion_name.='.exclusion';
	my $location=$project_dir.'/GEMIndices/'.$exclusion_name;

	my $count=$sth2->execute($annotation_id);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth2->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome $annotation\n";
		$present=1;
	    }
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    my ($ins_query,$ins_sth);
	    $ins_query ="INSERT INTO $table2 ";
	    $ins_query.='SET annotation_id = ?,species_id = ?,location = ?';
	    print $log_fh "Executing: $ins_query\n";
	    $ins_sth=$dbh->prepare($ins_query);
	    $ins_sth->execute($annotation_id,$species_id,$location);
	    $ins_sth->finish();
	}
	return($location);
    };

    # Set or get the information on the naming of the annotation tables that
    # correspond to the genome and annotation we are looking at.
    # Get the common sth for this table
    my ($ins_query3,$ins_sth3);
    $ins_query3 ="INSERT INTO $table3 ";
    $ins_query3.='SET genome_id = ?,annotation_id = ?,';
    $ins_query3.='type = ?,table_name = ?';
    $ins_sth3=$dbh->prepare($ins_query3);

    $guessing_subs{'junctionstable'}=sub {
	my $filetype=shift;
	my $present=0;
	my $annot_string=substr($annotation,0,30);
	my $table_name=$annotation_id.'_'.$annot_string;
	$table_name.='.junct';
	$table_name=~s/\./_/g;
	my $location=$table_name;

	# The table naming is based on the unique indices of the genome and
	# annotation so we don't need to check anything.
	my $count=$sth3->execute($annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth3->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $annotation\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "Table present at $location\n";
	} else {
	    print $log_fh "Table not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query3\n";
	    $ins_sth3->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };

    $guessing_subs{'junctionsclass'}=sub {
	my $filetype=shift;
	my $present=0;
	my $annot_string=substr($annotation,0,30);
	my $table_name=$annotation_id.'_'.$annot_string;
	$table_name.='.junctclass';
	$table_name=~s/\./_/g;
	my $location=$table_name;

	# The table naming is based on the unique indices of the genome and
	# annotation so we don't need to check anything.
	my $count=$sth3->execute($annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth3->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $annotation\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "Table present at $location\n";
	} else {
	    print $log_fh "Table not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query3\n";
	    $ins_sth3->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };

    $guessing_subs{'geneclass'}=sub {
	my $filetype=shift;
	my $present=0;
	my $annot_string=substr($annotation,0,30);
	my $table_name=$annotation_id.'_'.$annot_string;
	$table_name.='.geneclass';
	$table_name=~s/\./_/g;
	my $location=$table_name;

	# The table naming is based on the unique indices of the genome and
	# annotation so we don't need to check anything.
	my $count=$sth3->execute($annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth3->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $annotation\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "Table present at $location\n";
	} else {
	    print $log_fh "Table not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query3\n";
	    $ins_sth3->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };

    $guessing_subs{'transclass'}=sub {
	my $filetype=shift;
	my $present=0;
	my $annot_string=substr($annotation,0,30);
	my $table_name=$annotation_id.'_'.$annot_string;
	$table_name.='.transclass';
	$table_name=~s/\./_/g;
	my $location=$table_name;

	# The table naming is based on the unique indices of the genome and
	# annotation so we don't need to check anything.
	my $count=$sth3->execute($annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth3->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $annotation\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "Table present at $location\n";
	} else {
	    print $log_fh "Table not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query3\n";
	    $ins_sth3->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };

    $guessing_subs{'exonsclass'}=sub {
	my $filetype=shift;
	my $present=0;
	my $annot_string=substr($annotation,0,30);
	my $table_name=$annotation_id.'_'.$annot_string;
	$table_name.='.exclass';
	$table_name=~s/\./_/g;
	my $location=$table_name;

	# The table naming is based on the unique indices of the genome and
	# annotation so we don't need to check anything.
	my $count=$sth3->execute($annotation_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth3->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome $annotation\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "Table present at $location\n";
	} else {
	    print $log_fh "Table not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query3\n";
	    $ins_sth3->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };

    # Get the locations of the files for building the indices
    # Get a comon stk for the fasta files
    my ($ins_query5,$ins_sth5);
    $ins_query5 ="INSERT INTO $table5 ";
    $ins_query5.='SET genome_id = ?,annotation_id = ?,';
    $ins_query5.='filetype = ?,file_name = ?';
    $ins_sth5=$dbh->prepare($ins_query5);

    $guessing_subs{'exonsfasta'}=sub {
	my $filetype=shift;
	my $present=0;
	my $fasta_name=$annotation.'.'.$genome;
	$fasta_name=~s/\.g[tf]f//g;
	$fasta_name=~s/\.fa(stq)?$//;
	$fasta_name.='.exons.fa';
	my $location=$project_dir.'/GEMIndices/'.$fasta_name;
	
	my $count=$sth5->execute($annotation_id,$genome_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth5->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query5 for $location\n";
	    $ins_sth5->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };
    $guessing_subs{'transcriptomefasta'}=sub {
	my $filetype=shift;
	my $present=0;
	my $fasta_name=$annotation.'.'.$genome;
	$fasta_name=~s/\.g[tf]f//g;
	$fasta_name=~s/\.fa(stq)?$//;
	$fasta_name.='.transcripts.fa';
	my $location=$project_dir.'/GEMIndices/'.$fasta_name;
	
	my $count=$sth5->execute($annotation_id,$genome_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth5->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query5\n";
	    $ins_sth5->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };
    $guessing_subs{'junctionsfasta'}=sub {
	my $filetype=shift;
	my $present=0;
	my $fasta_name=$annotation.'.'.$genome;
	$fasta_name=~s/\.g[tf]f//g;
	$fasta_name=~s/\.fa(stq)?$//;
	$fasta_name.='.junctions.fa';
	my $location=$project_dir.'/GEMIndices/'.$fasta_name;
	
	my $count=$sth5->execute($annotation_id,$genome_id,$filetype);
	if ($count > 0) {
	    if ($count == 1) {
		($location)=$sth5->fetchrow_array();
		$present=1;
	    } else {
		warn "More than one entry retrieved for $genome\n";
		$present=1
	    }
	}
	if ($present) {
	    print $log_fh "File present at $location\n";
	} else {
	    print $log_fh "File not present. Will be built at $location\n";
	    # The entry is absent and we must set it
	    print $log_fh "Executing: $ins_query5\n";
	    $ins_sth5->execute($genome_id,$annotation_id,
			       $filetype,$location);
	}
	return($location);
    };

    my $guess=sub {
	my $key=shift;
	if ($guessing_subs{$key}) {
	    my $value=$guessing_subs{$key}->($key);
	} else {
	    print $log_fh "Haven't heard of $key\n";
	    return();
	}
    };

    return($guess);
}

# This script should for each of the necessary files
sub guess_values {
    my $guessmaster=shift;
    my $options=shift;

    print STDERR "Some values missing; time for some guesswork...\n";
    foreach my $key (keys %{$guessmaster}) {
	print STDERR "Guessing $key ...";
	my $value=guess($key);
	${$options->{$key}}=$value;
	print STDERR "done\n";
    }


    print STDERR "done\n";
}

## The following subs will build the tables for each of the steps in the
# pipeline
sub get_tables_hash {
    my $prefix=shift;
    # Build a hash with the table root as the key and the name as the value
    my %tables=('_start' => '',
		'_read_stats' => '',
		'_dataset' => '',
		'_junctions' => '',
		'_genes' => '',
		'_detected_genes' => '',
		'_transcripts' => '',
		'_detected_transcripts' => '',
		'_exon_seqs' => '',
		'_junction_seqs' => '',
		'_transcript_seqs' => '',
		'_indices' => '',
		'_genome_mapping' => '',
		'_transcriptome_mapping' => '',
		'_junctions_mapping' => '',
		'_unmapped_reads' => '',
		'_split_mapping' => '',
		'_recursive_mapping' => '',
		'_unique_maps_genome' => '',
		'_unique_maps_transcripts' => '',
		'_unique_maps_junctions' => '',
		'_unique_maps_split' => '',
		'_initial_clusters' => '',
		'_split_mapping_breakdown' => '',
		'_extract_annotation' => '',
		'_gene_coverage' => '',
		'_proj_coverage' => '',
		'_exon_coverage' => '',
		'_junction_coverage' => '',
		'_splicing_summary' => '',
		'_exon_RPKM_pooled' => '',
		'_gene_readcount_pooled' => '',
		'_gene_RPKM_pooled' => '',
		'_gene_RPKM_dist' => '',
		'_EJEI' => '',
		'_store_reads' => '',
		'_exon_classification' => '',
		'_junction_classification' => '',
		'_novel_junctions_summary' => '',
		'_all_junctions_class' => '',
		'_all_junctions_class_pooled' => '',
#		'_transcript_expression_levels' => '',
		'_transcript_expression_levels_pooled' => '',
#		'_exon_inclusion' => '',
		'_exon_inclusion_reads' => '',
		'_exon_inclusion_pooled' => '',
		'_inclusion_correlation' => '',
		'_inclusion_dist' => '',
		'_merged_SAM' => '',
		'_bed_files' => '',
		'_fusion_transcripts' => '',
		'_fusion_transcripts_support' => '',
		'_read_classification' => '',
		'_qualitiespos' => '',
		'_ambiguous' => '',
		'_read_dist_transcripts' => '',
		'_merged_mapping' => '',
		'_junction_maps_class' => '',
		'_completion_status' => '',
		'_summaries' => '',
		'_register_results' => '',
		'_gene_RPKM_pooled_flux' => ''
		);

    # Add the species prefix to the table
    foreach my $table (keys %tables) {
	$tables{$table}=$prefix.$table;
    }
    return(\%tables);
}

# This sub will create the mysql command files necessary for building the tables
sub build_table_files {
    my $prefix=shift;
    my $mysql_dir=shift;

    # Description of the tables to be created
    my %tables=('_start' => "DROP TABLE IF EXISTS ${prefix}_start;
CREATE TABLE ${prefix}_start (
    filename varchar(100),
    presence varchar(20)
    );",
		'_summaries' => "DROP TABLE IF EXISTS ${prefix}_summaries;
CREATE TABLE ${prefix}_summaries (
    summary varchar(120)
    );",
		'_completion_status' => "DROP TABLE IF EXISTS ${prefix}_completion_status;
CREATE TABLE ${prefix}_completion_status (
    step varchar(100),
    status varchar(20)
    );",
		'_inclusion_correlation' => "DROP TABLE IF EXISTS ${prefix}_inclusion_correlation;
CREATE TABLE ${prefix}_inclusion_correlation (
    LaneName1 varchar(50) NOT NULL,
    LaneName2 varchar(50) NOT NULL,
    Lane1Inc smallint unsigned NOT NULL,
    Lane2Inc smallint unsigned NOT NULL
    );",
		'_dataset' => "DROP TABLE IF EXISTS ${prefix}_dataset;
CREATE TABLE ${prefix}_dataset (
    lane_idx varchar(50) NOT NULL,
    lane_id varchar(50) NOT NULL,
    pair_id varchar(50) NOT NULL,
    sample_id varchar(50) NOT NULL,
    INDEX index_lane (lane_idx)
    );",
		'_fusion_transcripts' => "DROP TABLE IF EXISTS ${prefix}_fusion_transcripts;
CREATE TABLE ${prefix}_fusion_transcripts (
    readId varchar(50) not null,
    readlength mediumint unsigned not null,
    type varchar(10) not null,
    trans1mapstart varchar(20) not null,
    trans1 varchar(50) not null,
    chr1 varchar(20) not null,
    genomic1mapstart varchar(20) not null,
    gene1 varchar(50) not null,
    gene1rpkm float unsigned not null,
    trans2mapstart varchar(20) not null,
    trans2 varchar(50) not null,
    chr2 varchar(20) not null,
    genomic2mapstart varchar(20) not null,
    gene2 varchar(50) not null,
    gene2rpkm float unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_trans1 (trans1),
    INDEX idx_trans2 (trans2),
    INDEX idx_gene1 (gene1),
    INDEX idx_gene2 (gene2),
    INDEX idx_lane (LaneName)
    );",
		'_fusion_transcripts_support' => "DROP TABLE IF EXISTS ${prefix}_fusion_transcripts_support;
CREATE TABLE ${prefix}_fusion_transcripts_support (
    type varchar(10),
    location varchar(100),
    number mediumint unsigned not null,
    start_chr varchar(50) NULL,
    start int unsigned NULL,
    end_chr varchar(50) NULL,
    end int unsigned NULL,
    reference varchar(10) NOT NULL,
    LaneName varchar(20) NOT NULL,
    INDEX idx_location (location),
    INDEX idx_lane (LaneName)
    );",
		'_read_classification' => "DROP TABLE IF EXISTS ${prefix}_read_classification;
CREATE TABLE ${prefix}_read_classification (
       lane_id varchar(50) NOT NULL,
       Exonic double unsigned NOT NULL,
       ExFrac double unsigned NOT NULL,
       Intronic double unsigned NOT NULL,
       IntrFrac double unsigned NOT NULL,
       Intergenic double unsigned NOT  NULL,
       InterFrac double unsigned NOT  NULL,
       Total double unsigned NOT NULL,
       TptFrac double unsigned NOT NULL,
       index idx_lane (lane_id)
);",
#		'_exon_inclusion' => "DROP TABLE IF EXISTS ${prefix}_exon_inclusion;
#CREATE TABLE ${prefix}_exon_inclusion (
#       gene_id varchar(50) NOT NULL,
#       exon_id varchar(50) NOT NULL,
#       ExIncl double unsigned NOT NULL,
#       JuncInc double unsigned NOT NULL,
#       JuncExc double unsigned NOT  NULL,
#       inc_rate double unsigned NULL,
#       lane_id varchar(100) NOT NULL,
#       index idx_gene (gene_id),
#       index idx_exon (exon_id),
#       index idx_lane (lane_id)
#);",
		'_exon_inclusion_reads' => "DROP TABLE IF EXISTS ${prefix}_exon_inclusion_reads;
CREATE TABLE ${prefix}_exon_inclusion_reads (
       exon_id varchar(50) NOT NULL,
       chr varchar(50) NOT NULL,
       start int unsigned NOT NULL,
       end int unsigned NOT NULL,
       ExIncl int unsigned NOT NULL,
       JuncInc int unsigned NOT NULL,
       JuncExc int unsigned NOT  NULL,
       inc_rate double unsigned NULL,
       sample_id varchar(100) NOT NULL,
       index idx_exon (exon_id),
       index idx_sample (sample_id)
);",
		'_exon_inclusion_pooled' => "DROP TABLE IF EXISTS ${prefix}_exon_inclusion_pooled;
CREATE TABLE ${prefix}_exon_inclusion_pooled (
       gene_id varchar(50) NOT NULL,
       exon_id varchar(50) NOT NULL,
       ExIncl double unsigned NOT NULL,
       JuncInc double unsigned NOT NULL,
       JuncExc double unsigned NOT  NULL,
       inc_rate double unsigned NULL,
       sample varchar(50) NOT NULL,
       index idx_gene (gene_id),
       index idx_exon (exon_id),
       index idx_sample (sample)
);",
		'_inclusion_dist' => "DROP TABLE IF EXISTS ${prefix}_inclusion_dist;
CREATE TABLE ${prefix}_inclusion_dist (
    support mediumint unsigned not null,
    incl_percent smallint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_lane (LaneName)
);",
#		'_transcript_expression_levels' => "DROP TABLE IF EXISTS ${prefix}_transcript_expression_levels;
#CREATE TABLE ${prefix}_transcript_expression_levels (
#       gene_id varchar(50) NOT NULL,
#       transcript_id varchar(50) NOT NULL,
#       flux_locus varchar(50) NOT NULL,
#       rpkm double unsigned NULL,
#       lane_id varchar(100) NOT NULL,
#       index idx_gene (gene_id),
#       index idx_transcript (transcript_id),
#       index idx_flux_locus (flux_locus),
#       index idx_lane (lane_id)
#);",
		'_transcript_expression_levels_pooled' => "DROP TABLE IF EXISTS ${prefix}_transcript_expression_levels_pooled;
CREATE TABLE ${prefix}_transcript_expression_levels_pooled (
       gene_id varchar(50) NOT NULL,
       transcript_id varchar(50) NOT NULL,
       flux_locus varchar(50) NOT NULL,
       rpkm double unsigned NULL,
       sample varchar(100) NOT NULL,
       index idx_gene (gene_id),
       index idx_transcript (transcript_id),
       index idx_flux_locus (flux_locus),
       index idx_lane (sample)
);",
		'_junction_classification' => "DROP TABLE IF EXISTS ${prefix}_junction_classification;
CREATE TABLE ${prefix}_junction_classification (
       table_id varchar(100) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_exon_classification' => "DROP TABLE IF EXISTS ${prefix}_exon_classification;
CREATE TABLE ${prefix}_exon_classification (
       table_id varchar(100) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_store_reads' => "DROP TABLE IF EXISTS ${prefix}_store_reads;
CREATE TABLE ${prefix}_store_reads (
    filename varchar(100),
    uncompressed varchar(100) not null,
    md5uncomp char(32) not null,
    compressed varchar(100) not null,
    md5comp char(32) not null,
    md5global char(32) not null
    );",
		'_gene_readcount_pooled' => "DROP TABLE IF EXISTS ${prefix}_gene_readcount_pooled;
CREATE TABLE ${prefix}_gene_readcount_pooled (
    gene_id varchar(100) not null,
    readcount int unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_lane (LaneName)
    );",
		'_gene_RPKM_pooled' => "DROP TABLE IF EXISTS ${prefix}_gene_RPKM_pooled;
CREATE TABLE ${prefix}_gene_RPKM_pooled (
    gene_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_sample (sample)
    );",
		'_gene_RPKM_pooled_flux' => "DROP TABLE IF EXISTS ${prefix}_gene_RPKM_pooled_flux;
CREATE TABLE ${prefix}_gene_RPKM_pooled_flux (
    gene_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_sample (sample)
    );",
		'_gene_RPKM_dist' => "DROP TABLE IF EXISTS ${prefix}_gene_RPKM_dist;
CREATE TABLE ${prefix}_gene_RPKM_dist (
    support varchar(100) not null,
    RPKM mediumint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_lane (LaneName)
    );",
		'_exon_RPKM_pooled' => "DROP TABLE IF EXISTS ${prefix}_exon_RPKM_pooled;
CREATE TABLE ${prefix}_exon_RPKM_pooled (
    exon_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_exon (exon_id),
    INDEX idx_sample (sample)
    );",
		'_EJEI' => "DROP TABLE IF EXISTS ${prefix}_EJEI;
CREATE TABLE ${prefix}_EJEI (
    gene_id varchar(100) not null,
    junction_id varchar(100) not null,
    EJEI double unsigned not null,
    total int not null,
    LaneName varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_junction (junction_id),
    INDEX idx_lane (LaneName)
    );",
		'_merged_SAM' => "DROP TABLE IF EXISTS ${prefix}_merged_SAM;
CREATE TABLE ${prefix}_merged_SAM (
    LaneName varchar(100) not null,
    BAM_file varchar(100) not null,
    INDEX idx_lane (LaneName)
    );",

		'_junction_coverage' => "DROP TABLE IF EXISTS ${prefix}_junction_coverage;
CREATE TABLE ${prefix}_junction_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_splicing_summary' => "DROP TABLE IF EXISTS ${prefix}_splicing_summary;
CREATE TABLE ${prefix}_splicing_summary (
    junc_type varchar(50) not null,
    detected mediumint unsigned not null,
    total mediumint unsigned NULL
    );",
		'_novel_junctions_summary' => "DROP TABLE IF EXISTS ${prefix}_novel_junctions_summary;
CREATE TABLE ${prefix}_novel_junctions_summary (
    junction_id varchar(50) not null,
    chr varchar(50) not null,
    start int unsigned not null,
    end int unsigned not null,
    junc_type varchar(50) not null,
    support int unsigned not null,
    sample varchar(50) not null
    );",
		'_all_junctions_class' => "DROP TABLE IF EXISTS ${prefix}_all_junctions_class;
CREATE TABLE ${prefix}_all_junctions_class (
    chr1 varchar(50) not null,
    start int unsigned not null,
    chr2 varchar(50) not null,
    end int unsigned not null,
    junc_type varchar(50) not null,
    support int unsigned not null,
    exons1 text,
    exons2 text,
    LaneName varchar(50) not null
    );",
		'_all_junctions_class_pooled' => "DROP TABLE IF EXISTS ${prefix}_all_junctions_class_pooled;
CREATE TABLE ${prefix}_all_junctions_class_pooled (
    chr1 varchar(50) not null,
    start int unsigned not null,
    chr2 varchar(50) not null,
    end int unsigned not null,
    junc_type varchar(50) not null,
    support int unsigned not null,
    exons1 text,
    exons2 text,
    sample varchar(50) not null
    );",
		'_exon_coverage' => "DROP TABLE IF EXISTS ${prefix}_exon_coverage;
CREATE TABLE ${prefix}_exon_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_gene_coverage' => "DROP TABLE IF EXISTS ${prefix}_gene_coverage;
CREATE TABLE ${prefix}_gene_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_bed_files' => "DROP TABLE IF EXISTS ${prefix}_bed_files;
CREATE TABLE ${prefix}_bed_files (
    filename varchar(100),
    features int unsigned not null
    );",
		'_proj_coverage' => "DROP TABLE IF EXISTS ${prefix}_proj_coverage;
CREATE TABLE ${prefix}_proj_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_extract_annotation' => "DROP TABLE IF EXISTS ${prefix}_extract_annotation;
CREATE TABLE ${prefix}_extract_annotation (
    filename varchar(40),
    features mediumint unsigned not null
    );",
		'_split_mapping_breakdown' => "DROP TABLE IF EXISTS ${prefix}_split_mapping_breakdown;
CREATE TABLE ${prefix}_split_mapping_breakdown (
    filename varchar(100),
    type varchar(10),
    location varchar(100),
    number mediumint unsigned not null,
    start_chr varchar(50) NULL,
    start int unsigned NULL,
    end_chr varchar(50) NULL,
    end int unsigned NULL,
    INDEX idx_location (location)
    );",
		'_initial_clusters' => "DROP TABLE IF EXISTS ${prefix}_initial_clusters;
CREATE TABLE ${prefix}_initial_clusters (
    ClustertType varchar(100),
    chromosome varchar(40),
    ClusterNumber  mediumint unsigned not null
    );",
		'_unique_maps_split' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_split;
CREATE TABLE ${prefix}_unique_maps_split (
    filename varchar(100),
    type varchar(30),
    number mediumint unsigned not null
    );",
		'_unique_maps_junctions' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_junctions;
CREATE TABLE ${prefix}_unique_maps_junctions (
    filename varchar(100),
    type varchar(10),
    number mediumint unsigned not null
    );",
		'_read_dist_transcripts' => "DROP TABLE IF EXISTS ${prefix}_read_dist_transcripts;
CREATE TABLE ${prefix}_read_dist_transcripts (
    length_cat varchar(50) not null,
    start smallint unsigned not null,
    position tinyint unsigned not null,
    hits mediumint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_lane (LaneName),
    INDEX idx_length (length_cat)
    );",
		'_unique_maps_transcripts' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_transcripts;
CREATE TABLE ${prefix}_unique_maps_transcripts (
    filename varchar(100),
    transcript_id varchar(60),
    uniqueMaps mediumint unsigned NOT NULL,
    INDEX idx_trans_id (transcript_id)
    );",
		'_unique_maps_genome' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_genome;
CREATE TABLE ${prefix}_unique_maps_genome (
    filename varchar(100),
    chromosome varchar(40),
    uniqueMaps mediumint unsigned not null
    );",
		'_junction_maps_class' => "DROP TABLE IF EXISTS ${prefix}_junction_maps_class;
CREATE TABLE ${prefix}_junction_maps_class (
    junction_id varchar(100),
    type varchar(20),
    support mediumint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_junc (junction_id),
    INDEX idx_type (type),
    INDEX idx_lane (LaneName)
    );",
		'_split_mapping' => "DROP TABLE IF EXISTS ${prefix}_split_mapping;
CREATE TABLE ${prefix}_split_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_unmapped_reads' => "DROP TABLE IF EXISTS ${prefix}_unmapped_reads;
CREATE TABLE ${prefix}_unmapped_reads (
    filename varchar(100),
    unmappedReads mediumint unsigned not null
    );",
		'_merged_mapping' => "DROP TABLE IF EXISTS ${prefix}_merged_mapping;
CREATE TABLE ${prefix}_merged_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_junctions_mapping' => "DROP TABLE IF EXISTS ${prefix}_junctions_mapping;
CREATE TABLE ${prefix}_junctions_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_transcriptome_mapping' => "DROP TABLE IF EXISTS ${prefix}_transcriptome_mapping;
CREATE TABLE ${prefix}_transcriptome_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_genome_mapping' => "DROP TABLE IF EXISTS ${prefix}_genome_mapping;
CREATE TABLE ${prefix}_genome_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_recursive_mapping' => "DROP TABLE IF EXISTS ${prefix}_recursive_mapping;
CREATE TABLE ${prefix}_recursive_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_indices' => "DROP TABLE IF EXISTS ${prefix}_indices;
CREATE TABLE ${prefix}_indices (
    filename varchar(200) NOT NULL,
    type varchar(20) NOT NULL,
    created varchar(20) NOT NULL
    );",
		'_transcript_seqs' => "DROP TABLE IF EXISTS ${prefix}_transcript_seqs;
CREATE TABLE ${prefix}_transcript_seqs (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_junction_seqs' => "DROP TABLE IF EXISTS ${prefix}_junction_seqs;
CREATE TABLE ${prefix}_junction_seqs (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_genes' => "DROP TABLE IF EXISTS ${prefix}_genes;
CREATE TABLE ${prefix}_genes (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_detected_genes' => "DROP TABLE IF EXISTS ${prefix}_detected_genes;
CREATE TABLE ${prefix}_detected_genes (
       type varchar(200) NOT NULL,
       reliability varchar(10) NOT NULL,
       sample varchar(50) NOT NULL,
       detected int NOT NULL
);",
		'_transcripts' => "DROP TABLE IF EXISTS ${prefix}_transcripts;
CREATE TABLE ${prefix}_transcripts (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_detected_transcripts' => "DROP TABLE IF EXISTS ${prefix}_detected_transcripts;
CREATE TABLE ${prefix}_detected_transcripts (
       type varchar(200) NOT NULL,
       reliability varchar(10) NOT NULL,
       sample varchar(50) NOT NULL,
       detected int NOT NULL
);",
		'_exon_seqs' => "DROP TABLE IF EXISTS ${prefix}_exon_seqs;
CREATE TABLE ${prefix}_exon_seqs (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_junctions' => "DROP TABLE IF EXISTS ${prefix}_junctions;
CREATE TABLE ${prefix}_junctions (
       table_id varchar(100) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_qualitiespos' => "DROP TABLE IF EXISTS ${prefix}_qualitiespos;
CREATE TABLE ${prefix}_qualitiespos (
       LaneName varchar(100) NOT NULL,
       position smallint UNSIGNED NOT NULL,
       mean double UNSIGNED NOT NULL
);",
		'_ambiguous' => "DROP TABLE IF EXISTS ${prefix}_ambiguous;
CREATE TABLE ${prefix}_ambiguous (
       LaneName varchar(100) NOT NULL,
       position smallint UNSIGNED NOT NULL,
       ambiguous int UNSIGNED NOT NULL,
       A int UNSIGNED NOT NULL,
       C int UNSIGNED NOT NULL,
       G int UNSIGNED NOT NULL,
       T int UNSIGNED NOT NULL
);",
		'_read_stats' => "DROP TABLE IF EXISTS ${prefix}_read_stats;
CREATE TABLE ${prefix}_read_stats (
       Filename varchar(100)NOT NULL,
       Species varchar(50) NOT NULL,
       Project varchar(100) NOT NULL,
       ReadLength SMALLINT UNSIGNED NOT NULL,
       TotalReads INT UNSIGNED NOT NULL,
       NoAmbiguousBases INT unsigned NOT NULL,
       AmbiguousBases INT unsigned NOT NULL,
       UniqueReads INT unsigned NOT NULL,
       LaneName varchar(50) NOT NULL,
       pair_id varchar(50) NOT NULL
);",
		'_register_results' => "DROP TABLE IF EXISTS ${prefix}_register_results;
CREATE TABLE ${prefix}_register_results (
    filetype varchar(100),
    filename varchar(100),
    filemd5sum char(32) not null,
    registeredfile varchar(200),
    runmd5sum char(32) not null
    );",
	);

    my $working_dir=getcwd();
    unless ($mysql_dir) {
	$mysql_dir=$working_dir.'/mysql/table_build';
    }
    chdir($mysql_dir) ||	
	die "Could not chdir to $mysql_dir from $working_dir\n";
    foreach my $table (sort keys %tables) {	
	my $file_name=$prefix.$table.'.sql';
	my $table_fh=get_fh($file_name,1);
	print $table_fh $tables{$table},"\n";
	close($table_fh);
    }
    chdir($working_dir) ||
	die "Could not chdir from $mysql_dir to $working_dir\n";
    return(keys %tables);
}

# Create the directory structure necessary for the running of the pipeline
sub create_directory_structure {
    my $prefix=shift;
    my $log_fh=shift;
    my $options=shift;
    my $table_fh=shift;

    # First check the directory structure and create any directories that are
    # not there
    print $log_fh "Checking directory structure and attempting creation of missing dirs...\n";

    $options->{'localdir'}=~s/\/$//;
    my $localdir=${$options->{'localdir'}};
    print $log_fh "Setting localdir as: $localdir\n";
    my @directories=('mysql',
		     'mysql/table_build',
		     'mysql/table_data',
#		     'bin',
#		     'GEMIndices',
		     'genome',
		     'sequence',
		     'transcriptome',
		     'junctions',
		     'splitmapping',
		     'recursivemap',
		     'exons',
		     'readData',
		     'trimmedReads',
		     'graphs',
		     'logs',
		     "$localdir",
		     'clusters',
		     'exclusion',
#		     'results',
		     'SAM',
		     'work');
    # The -L option is necessary to print the same path for all nodes
    my $project_dir=getcwd();
    print $log_fh "Using dir: ".$project_dir,"\n";

    my %vars=('COMMONDB' => ${$options->{'commondb'}},
	      'DB'=> ${$options->{'database'}},
	      'PROJECT' => $project_dir,
	      'TABLES' => 'mysql/table_build',
	      'TAB_DAT' => 'mysql/table_data',
	      'BIN' => 'bin',
	      'LOGS' => 'logs',
	      'PAIRED' => ${$options->{'paired'}},
	      'READDIR' => 'readData',
	      'SEQUENCEDIR' => 'sequence', # (Raw data if not fa or fastq)
	      'GEMINDICES' => 'GEMIndices',
	      'SPLITMAPDIR' => 'splitmapping',
	      'RECMAPDIR' => 'recursivemap',
	      'GRAPHS' => 'graphs',
	      'SPECIES' => ${$options->{'species'}},
	      'PROJECTID' => ${$options->{'project'}},
	      'EXPID' => ${$options->{'experiment'}},
	      'PREFIX' => $prefix,
	      'ANNOTATION' => (abs_path(${$options->{'annotation'}}) || die "SORRY, I can't find the specified annotation file, please check the path"),
	      'GENOMESEQ' => (abs_path(${$options->{'genome'}}) || die "SORRY, I can't find the specified genome file, please check the path\n"),
	      'TRANSDIR' => 'transcriptome',
	      'EXONDIR' => 'exons',
	      'JUNCTIONSDIR' => 'junctions',
	      'MISMATCHES' => ${$options->{'mismatches'}},
	      'LOCALDIR' => ${$options->{'localdir'}},
	      'GENOMEINDEX' => ${$options->{'genomeindex'}},
	      'TRANSCRIPTOMEINDEX' => ${$options->{'transcriptomeindex'}},
	      'JUNCTIONSINDEX' => ${$options->{'junctionsindex'}},
	      'EXCLUSIONFILE' => ${$options->{'exclusionfile'}},
	      'EXONSFASTA' => ${$options->{'exonsfasta'}},
	      'JUNCTIONSFASTA' => ${$options->{'junctionsfasta'}},
	      'TRANSCRIPTOMEFASTA' => ${$options->{'transcriptomefasta'}},
	      'GENOMEDIR' => 'genome',
	      'STRANDED' => ${$options->{'stranded'}},
	      'READLENGTH' => ${$options->{'readlength'}},
	      'GENECLASSTABLE' => ${$options->{'geneclass'}},
	      'TRANSCLASSTABLE' => ${$options->{'transclass'}},
	      'JUNCTIONSTABLE' => ${$options->{'junctionstable'}},
	      'JUNCTIONSCLASSTABLE' => ${$options->{'junctionsclass'}},
	      'EXONSCLASSTABLE' => ${$options->{'exonsclass'}},
	      'FILELIST' => ${$options->{'files'}},
	      'THREADS' => ${$options->{'threads'}},
	      'MAPPER' => ${$options->{'mapper'}},
	      'CLUSTERDIR' => 'clusters',
	      'CLUSTER' => ${$options->{'cluster'}},
	      'SAMDIR' => 'SAM',
	      'LOCALPARALLEL' => 'work',
	      'QUALITIES' =>${$options->{'qualities'}},
	      'PREPROCESS' => ${$options->{'preprocess'}},
	      'RESULTS' => 'results',
	      'FLUXMEM' => ${$options->{'fluxmem'}}
);

    foreach my $dir (@directories) {
	print $log_fh "Checking $dir...\n";
	check_dir($dir,
		  $log_fh,
		  $table_fh);
	my $command="chmod g+x $dir";
	run_system_command($command,
			   $log_fh);
    }
    print $log_fh "done\n";
    return(\%vars);
}

sub check_dir {
    my $dir=shift;
    my $log_fh=shift;
    my $table_fh=shift;
    if (-e $dir) {
	print $log_fh "$dir already exists\n";
	print $log_fh "Skipping\n";
    } elsif ($dir=~/^bin/) {
	warn "THIS SHOULD NOT HAPPEN\n";
    } else {
	mkdir($dir) &&
	    print $log_fh "Creating dir $dir\n";
    }

    print $table_fh join ("\t",
			  $dir,
			  'present'),"\n";
    return();
}

## Print the necessary files for running the pipeline:
# These are the configuration file
sub print_config_file {
    my $vars=shift;
    my $log_fh=shift;
    my $conf_file='project_config.txt';
    my $conf_fh=get_fh($conf_file,1);

    print $log_fh "Creating the configuration file at $conf_file...";

    # Print the keys sorted as it makes it easier to check the file if necessary
    foreach my $var (sort keys %{$vars}) {
	my $value=$vars->{$var};
	if (defined $value) {
	    print $conf_fh join("\t",
				$var,
				$value),"\n";
	} else {
	    die "No value found for $var. Please check\n";
	}
    }
    close($conf_fh);
    print $log_fh "done\n";
}

# Print the pipeline file
sub print_pipeline_file {
    my $tables=shift;
    my $vars=shift;
    my $template=shift;
    my $log_fh=shift;
    my $pipfile='RNAseq_pipeline.txt';

    print $log_fh "Creating the master file at $pipfile...";

    my $template_fh=get_fh($template);
    my $new_pip_fh=get_fh($pipfile,1);

    while (my $line=<$template_fh>) {
	chomp($line);
	# Set the variable lines
	# Ignore comment lines
	if ($line!~/^\#/) {
	    # Set the variable lines
	    if ($line=~/=/) {
		my($variable,$assignment)=split(/\s*=\s*/,$line);
		my $value=$vars->{$variable};
		if ($assignment) {
		    $value=join('/',
				$vars->{'PROJECT'},
				$vars->{$variable});
		}
		print $log_fh join("\t",
				   $variable,
				   $value),"\n";
		$line=join(' = ',
			   $variable,
			   $value);
	    } elsif ($line=~/(^|\s+)(_\w+)/) {
		if (exists $tables->{$2}) {
		    # Table is present in the module. This has to be removed
		    # moving the tables gradually to the rules
		} elsif ($vars->{'PREFIX'}) {
		    # The table does not exist in the module, and has to be
		    # creted by the rule. When all the tables are migrated this
		    # warning should be removed
		    warn "WARNING: $2 Does not exist as a table\n";
		    my $prefix=$vars->{'PREFIX'};
		    $tables->{$2}=$prefix.$2;
		} else {
		    die "Something whent very wrong...\n";
		}
		$line=~s/(^|\s+)(_\w+)/$1$tables->{$2}/g;
	    }
	}
	print $new_pip_fh $line,"\n";
    }
    print $log_fh "done\n";
}

# Add the project and experiment information to the database
sub add_project {
    my $database=shift;
    my $project=shift;
    my $species=shift;
    my $log_fh=shift;

    my $dbh=MySQL_DB_Connect($database);

    print $log_fh "Checking for the presence of $project in the database...\n";
    my $table='projects';
    my ($query,$sth,$count);

    $query ='SELECT project_id ';
    $query.="FROM $table ";
    $query.="WHERE project_id = ?";
    $sth=$dbh->prepare($query);
    $count=$sth->execute($project);

    if ($count == 1) {
	print $log_fh "$project is present in the database\n";
    } elsif ($count > 1) {
	die "Project ID is present many times. This should not happen\n";
    } else {
	print $log_fh "$project is not present in the database\n";
	print $log_fh "Adding a new entry\n";

	# Insert the info into the database
	$query ="INSERT INTO $table ";
	$query.='SET species = ? , project_id = ?';
	print $log_fh "Executing: $query\n";
	$sth=$dbh->prepare($query);
	$sth->execute($species,$project);
    }
    return();
}

# Add the experiment information making sure the project id is present in the
# projects table
sub add_experiment {
    my $opts=shift;
    my $database=${$opts->{'commondb'}};
    my $project_id=${$opts->{'project'}};
    my $exp_id=${$opts->{'experiment'}};
    my $species=${$opts->{'species'}};
    my $annotation=${$opts->{'annotation'}};
    my $genome=${$opts->{'genome'}};
    my $template=${$opts->{'template'}};
    my $read_length=${$opts->{'readlength'}}; 
    my $preprocess_trim_length=${$opts->{'preprocess_trim_length'}};
    my $mismatches=${$opts->{'mismatches'}};

    # If the preprocess trim length is constant set the new read length
    if ($preprocess_trim_length=~/^=?(\d)+$/o) {
	my $trim=$preprocess_trim_length;
	$trim=~s/^=//;
	print STDERR "Adjusting read length by $trim\n";
	print STDERR "$preprocess_trim_length\n";
	$read_length -= $trim;
    }

    my $log_fh=shift;

    my $dbh=MySQL_DB_Connect($database);

    my $species_id=gset_field_id($dbh,
				 $species,
				 'species_info',
				 'species_id',
				 'species',
				 '',
				 $log_fh);
    my $genome_id=gset_field_id($dbh,
				$genome,
				'genome_files',
				'genome_id',
				'genome',
				$species_id,
				$log_fh);

    my $annotation_id=gset_field_id($dbh,
				    $annotation,
				    'annotation_files',
				    'annotation_id',
				    'annotation',
				    $species_id,
				    $log_fh);

    print $log_fh "Checking for the presence of $exp_id in the database...\n";
    my $table='experiments';
    my ($query,$sth,$count);

    $query ='SELECT experiment_id ';
    $query.="FROM $table ";
    $query.='WHERE experiment_id = ? AND project_id = ?';
    $sth=$dbh->prepare($query);
    $count=$sth->execute($exp_id,$project_id);

    my $update=0;
    if ($count == 1) {
	print $log_fh "$exp_id is already present in the database for $project_id\n";
	print STDERR "$exp_id is already present as part of $project_id. Are you
sure you want to overrite the previous information?(y/N)";
	my $yes=<STDIN>;
	chomp($yes);
	if ($yes=~/^y/i){
	    $update=1;
	    print $log_fh "Updating entry\n";
	    $query ="UPDATE $table ";
	    $query.='SET species_id = ?,';
	    $query.='genome_id = ?, annotation_id = ?, template_file = ?,';
	    $query.='read_length = ?, mismatches = ? ';
	    $query.='WHERE experiment_id = ? AND project_id = ?';
	}
    } elsif ($count > 1) {
	die "Project ID with experiment ID combination is present many times. This should not happen\n";
    } else {
	print $log_fh "$exp_id is not present in the database for $project_id\n";
	print $log_fh "Adding a new entry\n";
	$update=1;
	$query ="INSERT INTO $table ";
	$query.='SET species_id = ?,';
	$query.='genome_id = ?, annotation_id = ?, template_file = ?,';
	$query.='read_length = ?, mismatches = ?, experiment_id = ?,';
	$query.='project_id = ?';
    }

    if ($update) {
	# Insert the info into the database
	print $log_fh "Executing: $query\n";
	$sth=$dbh->prepare($query);
	$sth->execute($species_id,
		      $genome_id,$annotation_id,$template,
		      $read_length,$mismatches,$exp_id,$project_id);
    }
}

sub add_proj_info {
    my $opts=shift;
    my $database=${$opts->{'commondb'}};
    my $proj_id=${$opts->{'project'}};
    my $project_desc=${$opts->{'projdesc'}};

    my $log_fh=shift;

    my $dbh=MySQL_DB_Connect($database);

    if ($project_desc) {
	print $log_fh "Inserting description $project_desc of $proj_id into the database...\n";
    } else {
	print $log_fh "No project description supplied\n";
	return();
    }

    my $table='projects';
    my ($query,$sth,$count);

    $query ='SELECT proj_description ';
    $query.="FROM $table ";
    $query.='WHERE project_id = ? AND proj_description IS NOT NULL';
    $sth=$dbh->prepare($query);
    $count=$sth->execute($proj_id);

    my $overwrite=0;
    if ($count == 1) {
	print STDERR "$proj_id already has a description. Do you want to overwrite it?(y/n)\n";
	my $reply=<STDIN>;
	chomp($reply);
	if ($reply=~/^y/i) {
	    $overwrite=1;
	    $query ="UPDATE $table ";
	    $query.='SET proj_description = ? ';
	    $query.='WHERE project_id = ?';
	}
    } elsif ($count > 1) {
	die "Project ID with experiment ID combination is present many times. This should not happen\n";
    } else {
	$overwrite=1;
	$query ="INSERT INTO $table ";
	$query.='SET proj_description = ?,';
	$query.='project_id = ?';
    }

    if ($overwrite) {
	print $log_fh "Adding project description\n";

	# Insert the info into the database
	print $log_fh "Executing: $query\n";
	$sth=$dbh->prepare($query);
	$sth->execute($project_desc,$proj_id);
    }
}

sub add_exp_info {
    my $opts=shift;
    my $log_fh=shift;

    my $database=${$opts->{'commondb'}};
    my $proj_id=${$opts->{'project'}};
    my $exp_id=${$opts->{'experiment'}};
    my $preprocess_trim_length=${$opts->{'preprocess_trim_length'}};
    
    my %vals;
    $vals{'CellType'}=${$opts->{'cellline'}};
    $vals{'Compartment'}=${$opts->{'compartment'}};
    $vals{'exp_description'}=${$opts->{'run_description'}};
    $vals{'RNAType'}=${$opts->{'rnafrac'}};
    $vals{'Bioreplicate'}=${$opts->{'bioreplicate'}};

    # If the trimm length is not constant add the info to the description
    if ($preprocess_trim_length) {
	$vals{'exp_description'}.=" Original reads were trimmed by $preprocess_trim_length nucleotides";
    }
    
    my $dbh=MySQL_DB_Connect($database);

    my $table='experiments';
    my ($query,$sth,$count);
    
    my $overwrite_all=0;
    foreach my $key (keys %vals) {
	my $value=$vals{$key};
	unless ($value) {
	    print $log_fh "No $key supplied form $exp_id\n";
	    next;
	}
	
	$query ="SELECT $key ";
	$query.="FROM $table ";
	$query.='WHERE project_id = ? AND experiment_id = ? ';
	$sth=$dbh->prepare($query);
	$count=$sth->execute($proj_id,$exp_id);
	
	my $overwrite=$overwrite_all;
	if ($count == 1) {
	    my ($existing)=$sth->fetchrow_array() || '';
	    # If the keys are already the same or if the key is not provided
	    # skip
	    unless ($value) {next;}
	    if ($existing eq $value) {next;}

	    # If the key already has a value check if it should be altered
	    if ($existing) {
		print STDERR "$exp_id from $proj_id may already have a $key set to $existing. Do you want to overwrite $existing with $value?(y/n/all)\n";
		my $reply=<STDIN>;
		chomp($reply);
		if ($reply=~/^y/i) {
		    $overwrite=1;
		} elsif ($reply=~/^all$/i) {
		    $overwrite=1;
		    $overwrite_all=1;
		} else {
		    next;
		}
	    } else {
		$overwrite=1
	    }
	    $query ="UPDATE $table ";
	    $query.="SET $key = ? ";
	    $query.='WHERE project_id = ? and experiment_id = ?';
	} elsif ($count > 1) {
	    die "Project $proj_id with experiment $exp_id combination is present many times. This should not happen\n";
	} else {
	    $overwrite=1;
	    $query ="INSERT INTO $table ";
	    $query.="SET $key = ?,";
	    $query.='project_id = ? and experiment_id = ?';
	}

	if ($overwrite) {
	    print $log_fh "Setting $key for experiment $exp_id from $proj_id to $value \n";
	    # Insert the info into the database
	    print $log_fh "Executing: $query\n";
	    $sth=$dbh->prepare($query);
	    $sth->execute($value,$proj_id,$exp_id);
	}
    }
}

###
# Some subroutines to clean a directory and start from scratch
sub clear_tables {
    my $opts=shift;
    my $proj_id=${$opts->{'project'}};
    my $exp_id=${$opts->{'experiment'}};
    my $commondb=${$opts->{'database'}};
    my $prefix=$proj_id.'_'.$exp_id;

    my %tables=%{get_tables_hash($prefix)};

    foreach my $table (keys %tables) {
	my $value=$tables{$table};
	my $command="mysql $commondb -e 'DROP TABLE IF EXISTS $value'";
	run_system_command($command);
    }
}

# These are subs that will clear the directories in the project directory
### TO DO make an additionla level of clean so the genome directory and maybe
# some other ones are not completely wiped
sub clear_dirs {
    my $options=shift;

    # First check the directory structure and create any directories that are
    # not there
    print STDERR "Checking directory structure...\n";

    $options->{'localdir'}=~s/\/$//;
    my $localdir=${$options->{'localdir'}};
    print STDERR "Setting localdir as: $localdir\n";
    my @directories=('mysql',
#		     'GEMIndices',
		     'genome',
		     'transcriptome',
		     'junctions',
		     'splitmapping',
		     'recursivemap',
		     'exons',
		     'trimmedReads',
		     'graphs',
		     'logs',
		     "$localdir",
		     'clusters',
		     'exclusion',
		     'SAM',
		     'work');
    # The -L option is necessary to print the same path for all nodes
    my $project_dir=getcwd();
    print STDERR "Using dir: ".$project_dir,"\n";

    foreach my $dir (@directories) {
	my $command="rm -r $dir";
	if (-r $dir) {
	    run_system_command($command);
	}
    }
}

# Delete some files that may be left in the directory
sub clear_files {
    my $opts=shift;
    my $proj_id=${$opts->{'project'}};
    my $exp_id=${$opts->{'experiment'}};

    my $prefix=$proj_id.'_'.$exp_id;
    my $command;

    # Remove any table files that are left over
    $command="rm $prefix*";
    run_system_command($command);

    # Remove leftover log files
    $command='rm *.log';
    run_system_command($command);

    # remove the pipeline file
    $command='rm RNAseq_pipeline.txt';
    run_system_command($command);

    # remove possible temporary files left by the build sequences part
    $command='rm *.tmp.sequences.txt';
    run_system_command($command);
}

# Delete the information from the projects and experiments tables
sub clear_common_tables {
    my $opts=shift;

    my $proj_id=${$opts->{'project'}};
    my $exp_id=${$opts->{'experiment'}};
    my $commondb=${$opts->{'commondb'}};
    my $proj_tab='projects';
    my $exp_tab='experiments';

    my $dbh=MySQL_DB_Connect($commondb);

    my ($query,$sth,$count);

    # Delete the entry from the experiments table
    $query ="DELETE FROM $exp_tab ";
    $query.='WHERE project_id = ? AND experiment_id = ?';
    $sth=$dbh->prepare($query);
    print STDERR "Executing: $query\n";
    $sth->execute($proj_id,$exp_id);

    # Check if there are any entries left in the experiments table
    # corresponding to the project and if not remove it also from the
    # projects table
    $query ="SELECT experiment_id FROM $exp_tab ";
    $query.='WHERE project_id = ?';
    $sth=$dbh->prepare($query);
    print STDERR "Executing: $query\n";
    $count=$sth->execute($proj_id);

    if ($count > 0) {
	print STDERR $count,"\tExperiments left in $exp_tab corresponding to $proj_id\n";
    } else {
	print STDERR "No experiments left in $exp_tab corresponding to $proj_id\n";
	# Delete the entry from the projects table
	$query ="DELETE FROM $proj_tab ";
	$query.='WHERE project_id = ?';
	$sth=$dbh->prepare($query);
	print STDERR "Executing: $query\n";
	$sth->execute($proj_id);
    }
}

1;
