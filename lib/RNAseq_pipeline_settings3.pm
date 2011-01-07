package RNAseq_pipeline_settings3;

# Export to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('get_species_names_sub',
	    'read_config_file','read_file_list','get_dbh','check_db',
	    'get_numbers_from_annotation','get_saturation_curve',
	    'get_gff_chrom_freq','get_junction_no','get_unique_exons',
	    'get_unique_transcripts','get_common_dbh',
	    'get_gene_from_trans_sub','get_gene_from_exon_sub',
	    'is_paired','send2cluster',
	    'get_junction_type_sub','get_gene_from_junc_sub',
	    'get_gene_junctions_sub',
	    'get_gene_from_short_junc_sub',
	    'get_coords_from_exon_sub',
	    'get_coords_from_junc_id_sub',
	    'get_feature_overlap_sub','get_feature_overlap_split1',
	    'get_pair_id','get_lane_id','get_dataset_id','get_mapping_fh');

use strict;
use warnings;

# Load other modules required
use POSIX qw(uname);

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# This file contains subroutines used for building the tables in the pipeline as
# well as extracting data from them

# Load subs from other modules
use RNAseq_pipeline3 ('MySQL_DB_Connect','get_fh','run_system_command',
		      'parse_gff_line');
use Cwd;

# This sub will read the configuration settings from the file project_config.txt
### OK
sub read_config_file {
    my %configuration;
    my $config_file=shift;

    # Set a default for the configuration file if not provided
    unless ($config_file) {
	$config_file='project_config.txt';
    }

    # Read the configuration file and extract the values for the different
    # variables
    my $read_fh=get_fh($config_file);
    while (my $line=<$read_fh>) {
	chomp($line);
	my @line=split("\t",$line);
	if (@line != 2 ) {
	    die "Each line should be a tab separated key-value pair.\nThis is not the case for line:\t $line.\nPlease check the file $config_file";
	}
	$configuration{$line[0]}=$line[1];
    }
    close($config_file);

    return(\%configuration);
}

# This sub will read the read.list.txt file that contains a list of the
# files to be processed
### OK
sub read_file_list {
    my %files;

    my %options=%{read_config_file()};
    my $file_list=$options{'FILELIST'};

    # Read the information form the read.list.txt file
    my $read_fh=get_fh($file_list);
    while (my $line=<$read_fh>) {
	chomp($line);
	my @line=split("\t",$line);
	# This check is really not necessary, but it will give a better error
	# message than the one provided by the _dataset rule (which will check
	# for incorrecly formed read.list.txt files
	if (@line < 3 ) {
	    die "Each line in $file_list should be a tab separated list containing file name, pair ID, lane ID and optionally a group ID.\nThis is not the case for line:\t $line.\nPlease check the file $file_list";
	}
	my $file=shift(@line);
	$files{$file}=[@line];
    }
    close($file_list);

    return(\%files);
}

# Get the information on the pair id and the lane id
### TO DO
# These subs should be merged into one that uses the datasets table
sub get_pair_id {
    my $files=shift;
    my %pairs;

    foreach my $file (keys %{$files}) {
	my $filebase=$file;
	$filebase=~s/.fa(stq)?$//;
	push @{$pairs{$files->{$file}->[0]}},$filebase;
    }

    return(\%pairs);
}

sub get_lane_id {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	my $filebase=$file;
	$filebase=~s/.fa(stq)?$//;
	push @{$lanes{$files->{$file}->[1]}},$filebase;
    }

    return(\%lanes);
}

sub get_dataset_id {
    my $dbh=shift;
    my $table=shift;

    my %lanes;
    my ($query,$sth);

    $query ='SELECT pair_id, lane_idx ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($lane_id,$lane_index)=$sth->fetchrow_array()) {
	$lanes{$lane_id}=$lane_index;
    }

    return(\%lanes);
}

sub is_paired {
    my $pairids=shift;
    my $pair=shift;

    my $paired;
    if (@{$pairids->{$pair}} == 1) {
	$paired='single';
    } elsif (@{$pairids->{$pair}} == 2) {
	$paired='paired';
    } else {
	die "Incorrect number of pairs\n";
    }

    return($paired);
}


# This subroutine should build a sub that will return the gene to which an
# transcript belongs
sub get_gene_from_trans_sub {
    my %options=%{read_config_file()};
    my $dbh=get_dbh(1);
    my $table=$options{'EXONSCLASSTABLE'};

    # For saving time
    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT gene_id ';
    $query.="FROM $table ";
    $query.='WHERE transcript_id = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $trans=shift;

	unless ($cache{$trans}) {
	    $count=$sth->execute($trans);

	    if ($count != 1) {
		die "No gene in $table corresponds to $trans\n";
	    } else {
		my ($gene)=$sth->fetchrow_array();
		$cache{$trans}=$gene;
	    }
	}
	return($cache{$trans});
    };
    return($subroutine);
}

# This should get the gene id from an exon id
sub get_gene_from_exon_sub {
    my %options=%{read_config_file()};
    my $dbh=shift;
    my $table=$options{'EXONSCLASSTABLE'};

    # For saving time
    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT gene_id ';
    $query.="FROM $table ";
    $query.='WHERE exon_id = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $exon=shift;

	unless ($cache{$exon}) {
	    $count=$sth->execute($exon);

	    if ($count == 0) {
		die "No single gene in $table corresponds to $exon\n";
	    } else {
		my ($gene)=$sth->fetchrow_array();
		push @{$cache{$exon}}, $gene;
	    }
	}
	return($cache{$exon});
    };
    return($subroutine);
}

# This subroutine should build a sub that will return the gene to which an
# exon junction belongs
sub get_gene_from_junc_sub {
    my %options=%{read_config_file()};
    my $dbh=shift;
    my $table=$options{'JUNCTIONSTABLE'};

    # For saving time, as the junctions table is huge
    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT gene_id ';
    $query.="FROM $table ";
    $query.='WHERE junction_id = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $junc=shift;

	unless ($cache{$junc}) {
	    $count=$sth->execute($junc);

	    if ($count != 1) {
		die "No gene in $table corresponds to $junc\n";
	    } else {
		my ($gene)=$sth->fetchrow_array();
		$cache{$junc}=$gene;
	    }
	}
	return($cache{$junc});
    };
    return($subroutine);
}

# This sub should take a short junction_id and extract the gene it belongs to
# exon junction belongs
sub get_gene_from_short_junc_sub {
    my %options=%{read_config_file()};
    my $dbh=shift;
    my $table=$options{'JUNCTIONSTABLE'};

    # For saving time, as the junctions table is huge
    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT DISTINCT gene_id ';
    $query.="FROM $table ";
    $query.='WHERE chr = ? AND start = ? AND end = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $junc=shift;
	### TO DO fix for problematic chromosomes
	my ($chr,$start,$splice,$end)=split('_',$junc);

	unless ($cache{$junc}) {
	    $count=$sth->execute($chr,$start,$end);

	    if ($count == 0) {
		die "No gene in $table corresponds to $junc\n";
	    } else {
		while (my ($gene)=$sth->fetchrow_array()) {
		    push @{$cache{$junc}}, $gene;
		}
	    }
	}
	return($cache{$junc});
    };
    return($subroutine);
}

# This subroutine should build a sub that will return the junctions belonging to
# a certain gene
sub get_gene_junctions_sub {
    my %options=%{read_config_file()};
    my $dbh=shift;
    my $table=$options{'JUNCTIONSTABLE'};

    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT junction_id ';
    $query.="FROM $table ";
    $query.='WHERE gene_id = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $gene=shift;

	unless ($cache{$gene}) {
	    $count=$sth->execute($gene);

	    if ($count < 1) {
		die "No junctions in $table correspond to $gene\n";
	    } else {
		my @junctions;
		while (my ($junc)=$sth->fetchrow_array()) {
		    push @junctions,$junc;
		}
		$cache{$gene}=[@junctions];
	    }
	}
	return($cache{$gene});
    };
    return($subroutine);
}

# This subroutine should build a sub that will return the coordinates of an exon
# from the id of the exon
sub get_coords_from_exon_sub {
    my %cache;

    my $subroutine=sub {
	my $exon=shift;

	unless ($cache{$exon}) {
	    my @exon=split('_',$exon);
	    my $strand=pop(@exon);
	    my $end=pop(@exon);
	    my $start=pop(@exon);
	    my $chr=join('_',@exon);
	    if ($chr && $start && $end && $strand) {
		$cache{$exon}=[$chr,$start,$end,$strand];
	    } else {
		die "I can't identify $exon\n";
	    }
	}
	return($cache{$exon});
    };
    return($subroutine);
}

# This subroutine should build a sub that will return the coordinates of an
# junction from the id of the junction
sub get_coords_from_junc_id_sub {
    my %cache;

    my $subroutine=sub {
	my $junction=shift;

	unless ($cache{$junction}) {
	    my @junction=split('_',$junction);
	    my $end=pop(@junction);
	    my $start=pop(@junction);
	    my $chr=join('_',@junction);
	    if ($chr && $start && $end) {
		$cache{$junction}=[$chr,$start,$end];
	    } else {
		die "Unable to extract coordinates from $junction\n";
	    }
	}
	return($cache{$junction});
    };
    return($subroutine)
}

sub split_id {
    my $id=shift;
    $id=~s/splice//;
    $id=~s/_+/_/g;
    my @id=split('_',$id);
    my $end=pop(@id);
    my $start=pop(@id);
    my $chr=join('_',@id);
    if ($chr && $start && $end) {
	return([$chr,$start,$end])
    } else {
	die "Unable to extract coordinates from $id\n";
    }
}

# Get junction type from the database
sub get_junction_type_sub {
    my $table=shift;
    my $debug=shift || 0;
    my $dbh=get_dbh(1);
    my ($query,$sth,$count);
    my %types;

    # TO DO Select the table automatically

    $query ='SELECT DISTINCT type ';
    $query.="FROM $table ";
    $query.='WHERE chr = ? AND start = ? AND end = ?';
    $sth=$dbh->prepare($query);

    my $junction_type=sub {
	my $junction=shift;

	# This should avoid probles with junctions that have or do not have the
	# splice linker
	my $junction_id=$junction;
	$junction_id=~s/_splice//;
	# This should avoid problems with chromosomes containing underscores
	my @junction=split('_',$junction_id);
	my $end=pop(@junction);
	my $start=pop(@junction);
	my $chr=join('_',@junction);
	unless (exists $types{$junction}) {
	    $count=$sth->execute($chr,$start,$end);
	    unless ($count) {
		die "No junction type for $junction\n";
	    }
	    my %junctypes;
	    while (my ($junctype)=$sth->fetchrow_array()) {
		$junctypes{$junctype}++;
	    }

	    my @types=keys %junctypes;
	    my $type='';
	    if (@types > 1) {
		if ($junctypes{'known'}) {
		    $type='known';
		} else {
		    $type='unknown';
		}
		if ($debug) {
		    warn "More than one type for $junction. Classified as $type\n";
		}
	    } else {
		$type=shift(@types);
	    }
	    $types{$junction}=$type;
	}
	return($types{$junction});
    };

    return($junction_type);
}

# This subroutine will get the database associated to this specific project
sub get_dbh {
    my $common=shift;

    my $database;
    my %options=%{read_config_file()};
    my $host=$options{'HOST'};

    if ($common) {
	$database=$options{'COMMONDB'};
    } else {
	$database=$options{'DB'};
    }

    unless ($database && $host) {
	die "Both database and host are required in the config file\n";
    }

    my $dbh=MySQL_DB_Connect($database,
			     $host);
    if ($dbh) {
	return($dbh);
    } else {
	die "failed to connect to $database\n";
    }
}
    
# This subroutine will get the database containing the common information
### TO DO remove it and change every place where it appears to get_dbh(1)
sub get_common_dbh {
#    my %options=%{read_config_file()};
#    my $database=$options{'COMMONDB'};
#    my $host=$options{'HOST'};

#    unless ($database && $host) {
#	die "Both database and host are required in the config file\n";
#    }

#    my $dbh=MySQL_DB_Connect($database,
#			     $host);

    warn "This subroutine is deprecated. Use get_dbh(1)\n";
    my $dbh=get_dbh(1);
    return($dbh);
}

# This sub checks for the existence of a table
### OK
sub check_db {
    my $table=shift;
    my $outtable=shift;

    my %options=%{read_config_file()};
    if ($outtable) {
#	$outtable=$options{'PREFIX'}.'_junctions.txt';
	print STDERR "Setting outtable to $outtable\n";
    } else {
	die "No out table provided\n";
    }

    my $dbh=get_dbh(1);
    my $present=0;

    my ($query,$sth);
    $query ='SELECT count(*) ';
    $query.="FROM $table";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    my $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    
    my $tablefh=get_fh($outtable,1);
    if (@{$results}) {
	$present=1;
	# Print the table location for the junctions of this experiment
	print $tablefh join("\t",
			    $table,
			    "Present"),"\n";
    } else {
	# Continue, as the table must be created
	print STDERR $table,"\tIs not present\n";
	print $tablefh join("\t",
			    $table,
			    "Generated"),"\n";
    }
    close ($tablefh);

    return($present);
}

# This sub will return the name of the species by default,
# and if it is given an argument it will return the abbreviated name
sub get_species_names_sub {
    my %species=('Homo sapiens' => 'Hs',
		 'Pan troglodytes' => 'Pt',
		 'Mus musculus' => 'Mum',
		 'Rattus norvegicus' => 'Rn',
		 'Loxodonta africana' => 'La',
		 'Echinops telfairi' => 'Et',
		 'Macaca mulatta' => 'Mam',
		 'Oryctolagus cuniculus' =>'Oc',
		 'Danio rerio' => 'Dr',
		 'Pongo pygmaeus' => 'Pp',
		 'Microcebus murinus' => 'Mim',
		 'Otolemur garnettii' => 'Og',
		 'Takifugu rubripes' => 'Tr',
		 'Gallus gallus' => 'Gg',
		 'Monodelphis domestica' => 'Md',
		 'Tetraodon nigroviridis' => 'Tn',
		 'Xenopus tropicalis' => 'Xt',
		 'Gasterosteus aculeatus' => 'Ga',
		 'Drosophila melanogaster' => 'Dm',
		 'Acyrthosiphon pisum' => 'Ap');
    my %synonyms=('mouse' => 'Mus musculus',
		  'fly' => 'Drosophila melanogaster',
		  'human' => 'Homo sapiens',
		  'chimp' => 'Pan troglodytes',
		  'rat' => 'Rattus norvegicus',
		  'elephant' => 'Loxodonta africana',
		  'tenrec' => 'Echinops telfairi',
		  'rabbit' => 'Oryctolagus cuniculus',
		  'macaque' => 'Macaca mulatta',
		  'zebrafish' => 'Danio rerio',
		  'orangutan' => 'Pongo pygmaeus',
		  'chicken' => 'Gallus gallus',
		  'bushbaby' => 'Otolemur garnettii',
		  'opossum' => 'Monodelphis domestica',
		  'aphid' => 'Acyrthosiphon pisum');
    my $names_sub= sub {
	my $species=shift;
	my $abbrev=shift;
	my $ret_val;
	my $sp_low_case=lc($species);
	if (exists $species{$species}) {
	    $ret_val=$species;
	} elsif (exists $synonyms{$sp_low_case}) {
	    $ret_val=$synonyms{$sp_low_case};
	} else {
	    $ret_val=0;
	}
	if ($abbrev) {
	    $ret_val=$species{$ret_val};
	}
	return($ret_val);
    };
    return($names_sub);
}

# This should give us some info on the features in the annotation
sub get_numbers_from_annotation {
    my $type=shift;
    my %config=%{read_config_file()};
    my $dbh=get_dbh();
    my $table=$config{'PREFIX'}.'_extract_annotation';

    my ($query,$sth);
    $query ='SELECT features ';
    $query.="FROM $table ";
    $query.='WHERE filename = ?';
    $sth=$dbh->prepare($query);
    my $count=$sth->execute($type);
    unless ($count &&
	    ($count == 1)) {
	die "Incorrect number of records ($count) extracted from $table\n";
    }
    my ($number)=$sth->fetchrow_array();
    return($number);
}

sub get_unique_exons {
    my %config=%{read_config_file()};
    my $table=$config{'EXONSCLASSTABLE'};
    my $dbh=get_dbh(1);

    my ($query,$sth);
    $query ='SELECT count(DISTINCT exon_id) ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    my ($number)=$sth->fetchrow_array();

    return($number);
}

sub get_unique_transcripts {
    my %config=%{read_config_file()};
    my $table=$config{'EXONSCLASSTABLE'};
    my $dbh=get_dbh(1);

    my ($query,$sth);
    $query ='SELECT count(DISTINCT transcript_id) ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    my ($number)=$sth->fetchrow_array();

    return($number);
}

sub get_junction_no {
    my $type=shift;
    my %config=%{read_config_file()};
    my $table=$config{'JUNCTIONSTABLE'};

    my $dbh=get_dbh(1);

    my ($query,$sth);

    print STDERR "Extracting junctions data from $table\n";

    $query ='SELECT count(junction_id) ';
    $query.="FROM $table";
    if ($type) {
	$query.=' WHERE type = ?';
    }
    $sth=$dbh->prepare($query);
    if ($type) {
	$sth->execute($type);
    } else {
	$sth->execute();
    }

    my ($count)=$sth->fetchrow_array();
    return($count);
}

# Get saturation curve from a hash of hashes where the first key is the file
# an the second the feature contained in it
sub get_saturation_curve {
    my $detected=shift;
    my $cumulative={};
    my $left=$detected;
    my @curve;

    while (keys %{$left}) {
	($cumulative,$left)=get_next_saturation_point($cumulative,$left,\@curve);
    }

    return(\@curve);
}

sub get_next_saturation_point {
    my $cumulative=shift;
    my $left=shift;
    my $curve=shift;

    my %totals;
    my %new_cumulative;

    foreach my $file (keys %{$left}) {
	# Add the current cumulative
	foreach my $feature (keys %{$cumulative}) {
	    $totals{$file}{$feature}=1;
	}

	# Add the files genes
	foreach my $feature (keys %{$left->{$file}}) {
	    $totals{$file}{$feature}=1;
	}
    }

    # Check which file adds more
    my $old_feat_no=keys %{$cumulative};
    my $features=0;
    my $extra_features=0;
    my $next_file='';
    foreach my $file (keys %totals) {
	my $new_features=keys %{$totals{$file}};
	if ($new_features > $features) {
	    $features=$new_features;
	    %new_cumulative=%{$totals{$file}};
	    $next_file=$file;
	    $extra_features=$new_features - $old_feat_no;
	}
    }

    push @{$curve},[$next_file,
		    $features,
		    $extra_features];

    delete $left->{$next_file};

    return(\%new_cumulative,$left);
}

# Some subs to check the junctions
sub get_read_length {
    my %options=%{read_config_file()};
    my $dbh=get_dbh();
    my ($query,$sth,$count);
    my $table=$options{'PREFIX'}.'_read_stats';

    $query ='SELECT ReadLength ';
    $query.="FROM $table";

    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless ($count &&
	    ($count > 0)) {
	die "No read length obtained from $table\n;";
    }
    my $min_length;
    
    while (my ($length)=$sth->fetchrow_array()) {
	if ($min_length) {
	    if ($length != $min_length) {
		warn "Reads with different lengths found, choosing shortest and continuing\n";
	    }
	    if ($length < $min_length) {
		$min_length=$length;
	    }
	} else {
	    $min_length=$length;
	}
    }
    $sth->finish();
    $dbh->disconnect();
    return($min_length);
}

sub get_gff_chrom_freq {
    my $readlength=get_read_length();
    my $subroutine= sub {
	my $infile=shift;

	my $infh=get_fh($infile);
	my %chromosomes;
	my $short=0;

	while (my $line=<$infh>) {
	    chomp($line);
	    my @line=split("\t",$line);
	    my $chr=$line[0];
	    my $length=$line[4] - $line[3] + 1;
	    if ($length == $readlength) {		
		$chromosomes{$chr}++;
	    } elsif ($length > $readlength) {
		warn "Mapped read is longer then $readlength for $chr\n";
		$chromosomes{$chr}++;
	    } else {
		$short++;
	    }
	}
	close($infh);

	if ($short) {
	    print STDERR $short,"\tRead maps were shorter than $readlength\n";
	}

	return(\%chromosomes);
    };
    return($subroutine)
}

sub get_mapping_fh {
    my $pair=shift;
    my $type=shift;

    my %options=%{read_config_file()};
    my $genomedir=$options{'GENOMEDIR'};
    my $junctionsdir=$options{'JUNCTIONSDIR'};
    my $splitdir=$options{'SPLITMAPDIR'};
    my $recursivedir=$options{'RECMAPDIR'};

    my %endings=('genome' => '.gem.map',
		 'junctions' => '.gem.map',
		 'junctionsgenomic' => '.gem.map.gen.coords',
		 'split' => '.unmapped.gem.split-map',
		 'recursive' => '.unmapped.gem.split-map.gz.recursive.map');

    my %directories=('genome' => $genomedir,
		     'junctions' => $junctionsdir,
		     'junctionsgenomic' => $junctionsdir,
		     'split' => $splitdir,
		     'recursive' => $recursivedir);

    my $filename=$directories{$type}.'/'.$pair.$endings{$type};

    # Check if the file (or a compressed version of it) exists
    my $fh;
    if (-r $filename) {
	print STDERR "Opening $filename\n";
	$fh=get_fh($filename);
    } elsif (-r $filename.'.gz') {
	print STDERR $filename,"\tIs gzipped, changing to $filename.gz\n";
	$filename.='.gz';
	$fh=get_fh($filename);
    } elsif($type eq 'recursive') {
	print STDERR "Recursive mapping is not present merging the rest\n";
	$fh=get_fh($filename);
#	return(0);
    } else {
	die "I can't find $filename or $filename.gz\n";
    }
    
    return($fh);
}

### Sub to send stuff to the cluster
sub send2cluster {
    my $file=shift;
    my $queue=shift;
    my $jobname=shift;
    my %options=%{read_config_file()};

    unless ($queue) {
	$queue=$options{'CLUSTER'};
    }

    my $logs=$options{'LOGS'};

    # Submit to the cluster
    my $command="qsub -q $queue $file";
    sleep(1);
    my $execute=`$command`;

    # monitor the execution
    my ($job_id)=(split(/\s+/,$execute))[2];
    $job_id=~s/\..+$//;
    print STDERR "Running job $job_id on $queue\n";
    
    $command="qstat -j $job_id 2>&1";
    print STDERR "Waiting for job $job_id";
    while (1) {
	my @running=`$command`;
	my $success=$?;
	my $finished=1;

	while (my $line=shift(@running)) {
	    if ($line=~/^=+$/) {
		next;
	    } elsif ($line=~/^job_number/o) {
		# The job is queued or running
		my @line=split(/\s+/,$line);
		$finished=0;
		print STDERR '.';
		last;
	    } elsif ($line=~/^Following\sjobs\sdo\snot\sexist/o) {
		# The job has finished
		last;
	    } else {
		# There is a problem
		print STDERR $line,"\n";
		print STDERR "\nProblem\n",@running,"\n";
		die;
	    }
	}
	
	if ($finished) {
	    print STDERR "done\n";
	    last;
	} else {
	    sleep(10);
	}
    }

    # collect the output of the cluster node into the log file for the
    # step
    # First indicate which job the output beloiongs to:
    $command="echo $job_id >> $logs/$jobname.$queue.log";
    run_system_command($command);

    # Add the command output
    $command="cat $jobname.[eo]$job_id* >> $logs/$jobname.$queue.log";
    run_system_command($command);

    # Remove the files
    $command="rm $jobname.[eo]$job_id*";
    run_system_command($command);

    return($job_id);
}

###
# Here are a series of subroutines used to run sarahs overlap program on files
# regardless of their size
### TO DO
# We need to make sure the correct queue is chosen
sub get_feature_overlap_sub {
    my $parallel=shift;
    my $paralleltmp=shift;
    my $bindir=shift;
    my $m_value=shift;
    my $size=1000000;
    my %subs;

    my $flags='-v -ucsc';
    if ($m_value) {
	$flags.=" $m_value";
    } else {
	$flags.=' -m -10';
    }
    
    $subs{'default'}= sub {
	my $features=shift;
	my $annotation=shift;
	my $stranded=shift;
	my $outfile=shift;
	my $tmpdir=shift;

	if ($stranded) {
#	print STDERR "WARNING:Overlapping as Unstranded due to -m -st problem in overlap\n"
	    $flags.=' -st 1';
	}

	# Get the otput file name
	my $overlapfn=$outfile;
	$overlapfn=~s/.gz//;
	$overlapfn.='.overlap.gz';
	
	# Check if the overlap file exists
	if (-r $overlapfn) {
	    print STDERR $overlapfn,"\tAlready exists. Skipping...\n";
	    return();
	}

	# run overlap
	print STDERR "Running overlap on $annotation and $features\n";
	
	# Split the file
	if (-r $features) {
	    print STDERR "Splitting $features into $size lines fragments\n";
	} else {
	    die "Can't read $features\n";
	}
	my @files=split_file($features,
			     $size,
			     $tmpdir);
	my @overlap_files;
	
	# run overlap for each of them
	foreach my $file (@files) {
	    my $outfn=$file.'.overlap';
	    run_overlap($annotation,
			$file,
			$flags,
			$outfn);
	    push @overlap_files, $outfn;
	}
	
	# Combine the overlap files
	my @over_files=combine_overlap(\@overlap_files,
				       $overlapfn);
	
	# clean up
	clean_up(@files,
		 @over_files);
    };


    $subs{'parallel'}=sub {
	my $features=shift;
	my $annotation=shift;
	my $stranded=shift;
	my $outfile=shift;
	my $tmpdir=$paralleltmp;

	# Get the otput file name
	my $overlapfn=$outfile;
	$overlapfn=~s/.gz//;
	$overlapfn.='.overlap.gz';

	# Check if the overlap file exists
	if (-r $overlapfn) {
	    print STDERR $overlapfn,"\tAlready exists. Skipping...\n";
	    return();
	}

	# Check which responding nodes nodes have a threshold lower than 4
	my @available_clients=get_client_loads();

	# Check if we are running in a local directory, and if so restrict
	# clients to those that are local.
	my $localparalleldir=0;
	if ($tmpdir!~/^\/users/) {
	    my $runninghost=determine_host();
	    $localparalleldir=1;
	    print STDERR "WARNING: The tmp dir suplied for parallel jobs is local\n";
	    print STDERR "I will restrict available nodes to those from $runninghost\n";
	    my @localclients;
	    foreach my $client (@available_clients) {
		if ($runninghost=~/$client/) {
		    push @localclients, $client;
		}
	    }
	    @available_clients=@localclients;

	    my $nodes=@available_clients;
	    print STDERR $nodes,"\tnodes available locally on $runninghost\n";
	}

	my $client_number=0;
	$client_number=@available_clients;

	if ($client_number==0) {
	    die "There are no free nodes\n";
	}

	if ($stranded) {
#	print STDERR "WARNING:Overlapping as Unstranded due to -m -st problem in overlap\n"
	    $flags.=' -st 1';
	}

	# run overlap
	print STDERR "Running overlap on $annotation and $features\n";
	
	# Split the file
	if (-r $features) {
	    print STDERR "Splitting $features into $size lines fragments\n";
	} else {
	    die "Can't read $features\n";
	}
	my @files=split_file($features,
			     $size,
			     $tmpdir);
	my @overlap_files;
	
	# run overlap for each of them
	my @left=@files;
	my @pids;
	my %nodes;
	while(@left ||
	      (@available_clients < $client_number)) {
	    # While will not exit until all jobs are done and all clients are
	    # free
	    if (@available_clients && @left) {
		my $file=shift(@left);
		my $node=shift(@available_clients);
		my $child_pid=fork();
		$nodes{$child_pid}=$node;
		my $outfn=$file.'.overlap';
		# Temporary files seem named by overlap using the time, so 
		# if we don't sleep we can have problems with this
		sleep(1);
		if ($child_pid) {
		    push @pids, $child_pid;
		    push @overlap_files, $outfn;
		} else {  # child exec date		
		    print STDERR "Procesing $file in $node\n";
		    run_overlap($annotation,
				$file,
				$flags,
				$outfn,
				$node);
		}
	    } else {
		# Put in a while loop in case a child dies while another corpse
		# is being collected, to avoid leaving zombies
		my $ended=wait();
		sleep(1);
		if ($ended > 0) {
		    push @available_clients, $nodes{$ended};
		    print STDERR "Job\t",$ended,"\tFinished in $nodes{$ended}\n";
		}
	    }
	    sleep(1);
	}
	
	sleep(1);
	print STDERR "All partial jobs finished for $overlapfn\n";
	
	# Combine the overlap files
	my @over_files=combine_overlap(\@overlap_files,
				       $overlapfn);
	
	# clean up
	clean_up(@files,
		 @over_files);
    };

   $subs{'cluster'}=sub {
	my $features=shift;
	my $annotation=shift;
	my $stranded=shift;
	my $outfile=shift;
	my $tmpdir=$paralleltmp;

	# Get the otput file name
	my $overlapfn=$outfile;
	$overlapfn=~s/.gz//;
	$overlapfn.='.overlap.gz';

	# Check if the overlap file exists
	if (-r $overlapfn) {
	    print STDERR $overlapfn,"\tAlready exists. Skipping...\n";
	    return();
	}

	if ($stranded) {
	    $flags.=' -st 1';
	}

	# run overlap
	print STDERR "Running overlap on $annotation and $features\n";
	
	# Split the file
	if (-r $features) {
	    print STDERR "Splitting $features into $size lines fragments\n";
	} else {
	    die "Can't read $features\n";
	}
	my @files=split_file($features,
			     $size,
			     $tmpdir);

	# Create the file pairs fro overlap
	my @file_pairs;
	my @overlap_files;
	foreach my $file (@files) {
	    my $outfn=$file.'.overlap';
	    push @file_pairs,[$file,$outfn];
	    push @overlap_files, $outfn;
	}


	# build the submission file
	my $jobname='RNAseqOverlap';
	my $subfile=build_run_overlap_submission(\@file_pairs,
						 $bindir,
						 $flags,
						 $annotation,
						 $jobname);

	# submit to the cluster
	my $queue;
	my $job_id=send2cluster($subfile,
				$queue,
				$jobname);

	# clean up
	my $command="rm $subfile";
	print STDERR "Executing: $command\n";
	system($command);
	
	# Combine the overlap files
	my @over_files=combine_overlap(\@overlap_files,
				       $overlapfn);
	
	# clean up
	clean_up(@files,
		 @over_files);

	# Remove the log files from the cluster
	# This should no longer be necessary with the new send2cluster
	$command="rm RNAseqOverlap.?".$job_id."*";
	print STDERR "Executing: $command\n";
#	system($command);
    };

    if ($subs{$parallel}) {
	return($subs{$parallel});
    } else{
	print STDERR "Unknown type $parallel. Setting to default (serial)\n";
	return($subs{'default'});
    }
}

sub split_file {
    my $file=shift;
    my $size=shift;
    my $tmpdir=shift;
    my @splitfiles;


    my $infh=get_fh($file);
    my $lineno=0;
    my $fileno=0;
    my $outfn=$file.'.split.'.$fileno;
    if ($tmpdir) {
	$outfn=~s/.*\///;
	$tmpdir=~s/\/$//;
	$outfn=$tmpdir.'/'.$outfn;
    }
    my $outfh=get_fh($outfn,1);
    push @splitfiles,$outfn;
    while (my $line=<$infh>) {
	$lineno++;
	unless ($lineno % $size) {
	    close($outfh);
	    $fileno++;
	    $outfn=$file.'.split.'.$fileno;
	    if ($tmpdir) {
		$outfn=~s/.*\///;
		$tmpdir=~s/\/$//;
		$outfn=$tmpdir.'/'.$outfn;
	    }
	    push @splitfiles,$outfn;
	    $outfh=get_fh($outfn,1);
	}
	print $outfh $line;
    }
	    
    @splitfiles=sort @splitfiles;

    return(@splitfiles);
}

sub run_overlap {
    my $annotation=shift;
    my $features=shift;
    my $flags=shift;
    my $outfile=shift;
    my $host=shift;

    my $command;
    if ($host) {
	$command="ssh $host nice overlap $annotation $features $flags > $outfile";
	print STDERR "Executing: $command\n";
	exec($command);
    } else {
	$command="overlap $annotation $features $flags > $outfile";
	print STDERR "Executing: $command\n";
	system($command);
    }
}

sub combine_overlap {
    my $files=shift;
    my $overlapfn=shift;
    
    # Because we are running overlap on the same file there will always be
    # the same number of lines
    my @lines;
    my %feats;
    my @overlap_files;
    print STDERR "Combining overlap files...\n";
    foreach my $file (@{$files}) {
	my $overlapfile=$file;
	if (-r $overlapfile) {
	    print STDERR "Processing $overlapfile\n";
	} else {
	    warn "Unable to find readable file $overlapfile\n";
	}
	push @overlap_files,$overlapfile;

	my $overfh=get_fh($overlapfile);
	my $line_no=0;
	while (my $line=<$overfh>) {
	    chomp($line);
	    my %line=%{parse_gff_line($line)};

	    # Add the number of overlaps
	    my $hits;

	    # Add the locations if they are present
	    $hits=$line{'feature'}{'list_feat2:'};
	    if ($hits ne '.') {
		$feats{$line_no}[1].=$hits;
	    }
	    $feats{$line_no}[0]+=$line{'feature'}{'nb_ov_feat2:'};
	    unless ($lines[$line_no]) {
		$lines[$line_no]=$line;
	    }

	    $line_no++;
	}
	close($overfh);
    }

    print STDERR "Merging...\n";
    my $overlapfh=get_fh($overlapfn,1);
    for (my $i=0;$i<@lines;$i++) {
	my $line=$lines[$i];
	my $feats=$feats{$i}[0];
	my $hits=$feats{$i}[1];

	# Add a '.' for the cases with no hits in any file
	unless ($hits) {
	    $hits='.';
	}

	if ($hits) {
	    unless ($line=~s/(nb_ov_([^\s])+: )"\d+";( list_([^\s])+: )"([^\s])+";$/$1"$feats";$3"$hits";/) {
		warn "Problem with overlap line\n";
	    }
	} else {
	    warn "Problem no hits field in: $line\n";
	}

	print $overlapfh $line,"\n";
    }
    print STDERR "done\n";
    return (@overlap_files);
}

sub clean_up {
    my @files=@_;
    print STDERR "Cleaning up...\n";
    my $initial=@files;
    my $result;
    print STDERR $initial,"\tFiles to be removed\n";
    $result=unlink @files;
    print STDERR 'unlink ',@files,"\n";
    print STDERR $result,"\tFiles removed\n";
    if ($initial == $result) {
	print STDERR "successfully removed temporary files\n";
    } else {
	warn "Problem removing files\n";
    }
}

# For parallel running
sub get_client_loads {
    my @out = ();
    my $TIMEOUT = 10;
    my $threshold = 7;
    my %thresholds=('corb'=> 4,
		    'icarus'=> 4,
		    'cel'=> 4,
		    'foc' => 4,
		    'sun' => 4,
		    'cuc' => 4,
		    'tmatik' => 4, 
		    'palm' => 4,
		    'hulk' => 4,
		    'tofo' => 4
	);
    my @clients = ('cel',
#		   'corb',
		   'foc',
#		   'palm',
		   'tofo',
		   'hulk',
		   'sun',
		   'icarus', 
		   'cuc', 
		   'tmatik');

    # check load on all clients...
    print STDERR "Checking client load for the last 1/5/15 minute(s)\n";
    foreach my $client (@clients) {
	my $load = '';
	eval {
	    local $SIG{ALRM} = sub { die 'alarm time out' };
	    alarm $TIMEOUT;
	    $load = `ssh $client uptime`;
	    alarm 0;
	    1; # return value from eval on normalcy
	} or warn "load check on $client timed out after $TIMEOUT seconds.\n";
	unless ($load) {
	    next;
	}
	$load =~ /load average: (.*?), (.*?), (.*?)\s+/;
	my ($min_1,$min_5,$min_15) = ($1,$2,$3);
	my $max = 0;
	$max = $min_1 if ($min_1 > $max);
	$max = $min_5 if ($min_5 > $max);
	$max = $min_15 if ($min_15 > $max);
	while ($max < $thresholds{$client}) {
	    push @out, $client;
	    $thresholds{$client}--;
	}
    }

    # # Print out the available resources
    my $client_number=@out;
    print STDERR $client_number,"\tClients available\n";
    my %free;
    foreach my $client (@out){
	$free{$client}++;
    }
    foreach my $client (keys %free) {
	print STDERR join("\t",
			  $client,
			  $free{$client}),"\n";
    }

    return (@out)
}

sub determine_host {
    my ($host)=(uname)[1];
    return($host);
}

# For cluster running
sub build_run_overlap_submission {
    my $pairs=shift;
    my $bindir=shift;
    my $flags=shift;
    my $annotation=shift;
    my $jobname=shift;

    unless ($jobname) {
	die "No jobname provided\n";
    }

    print STDERR 'Building submission file...';
    my $filenum=@{$pairs};
     
    unless(@{$pairs}) {
	die "No input supplied\n";
    }
    
    # Get the input and output files
    my @infiles;
    my @outfiles;
    foreach my $pair (@{$pairs}) {
	push @infiles,$pair->[0];
	push @outfiles,$pair->[1];
    }

    # Print the submission file
    my $subfile="subfile.$$.job";
    my $outfh=get_fh($subfile,1);
    
    print $outfh <<FORMEND;
# Get the job name
#\$ -N $jobname
    
# Set the array jobs
#\$ -t 1-$filenum

# Request 8 cpus this cannot be done, but we can request memmory
#\$ -l h_vmem=16G

# Write in to the current working directory
#\$ -cwd 
export PATH=\$PATH:/soft/bin
infiles=(@infiles)
outfiles=(@outfiles)

export infile=\${infiles[\$SGE_TASK_ID-1]}
export outfile=\${outfiles[\$SGE_TASK_ID-1]}

echo \$HOSTNAME >&2
$bindir/overlap $annotation \$infile $flags > \$outfile
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile);
}

# This has to be rewritten for the cluster
sub get_feature_overlap_split1 {
    my $file1=shift;
    my $file2=shift;
    my $stranded=shift;
    my $outfile=shift;
    my $tmpdir=shift;
    my $size=1000000;
    my $flags='-v -m 0 -ucsc';

    if ($stranded) {
#	print STDERR "WARNING:Overlapping as Unstranded due to -m -st problem in overlap\n"
	$flags.=' -st 1';
    }

    # run overlap
    print STDERR "Running overlap on $file1 and $file2\n";

    # Split the file
    if (-r $file1) {
	print STDERR "Splitting $file1 into $size lines fragments\n";
    } else {
	die "Can't read $file1\n";
    }
    my @files=split_file($file1,
			 $size,
			 $tmpdir);
    my @overlap_files;
    
    # run overlap for each of them
    foreach my $file (@files) {
	my $outfn=$file.'.overlap';
	run_overlap($file,
		    $file2,
		    $flags,
		    $outfn);
	push @overlap_files, $outfn;
    }

    # Combine the overlap files
    my $overlapfn=$outfile;
    $overlapfn=~s/.gz//;
    $overlapfn.='.overlap.gz';
    my $files=join(' ',@overlap_files);
    my $command="cat $files | gzip -9 -c > $overlapfn";
    system($command);

    # clean up
    clean_up(@files,
	     @overlap_files);

}

1;
