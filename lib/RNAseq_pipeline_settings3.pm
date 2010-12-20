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
	    'is_paired',
	    'get_junction_type_sub','get_gene_from_junc_sub',
	    'get_gene_junctions_sub',
	    'get_gene_from_short_junc_sub',
	    'get_coords_from_exon_sub',
	    'get_coords_from_junc_id_sub',
	    'get_pair_id','get_lane_id','get_dataset_id','get_mapping_fh');

use strict;
use warnings;

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
use RNAseq_pipeline3 qw(MySQL_DB_Connect get_fh);
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
1;
