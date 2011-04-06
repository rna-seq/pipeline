package RNAseq_pipeline3;
# DGK 2008-2009 CRG

# Export subroutines to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');
push @EXPORT_OK,('MySQL_DB_Connect','get_fh','get_log_fh');
push @EXPORT_OK,('check_table_existence','check_file_existence');
push @EXPORT_OK,('check_field_existence');
push @EXPORT_OK,('check_gff_file','check_fasta_file','get_md5sum');
push @EXPORT_OK,('run_system_command');
push @EXPORT_OK,('print_gff','get_sorted_gff_fh','parse_gff_line');
push @EXPORT_OK,('get_files_from_table_sub');
push @EXPORT_OK,('get_annotation_from_gtf','get_exon_list_from_gtf',
		 'build_annotated_junctions');
push @EXPORT_OK,('cluster_gff','get_gff_from_junc_id','get_sorted_gff_fh',
		 'gff_line_check','format_as_gff');
push @EXPORT_OK,('get_chr_subseq','get_feature_coverage');
push @EXPORT_OK,('get_feature_overlap');
push @EXPORT_OK,('get_Average',
		 'get_exon_mappable_length_sub','get_exon_sequences',
		 'get_exon_coverage_1000nt');
push @EXPORT_OK,('get_distinct_col_value_from_tab','plot_heatmap');
push @EXPORT_OK,('get_feat_mappable_length_sub');

use strict;
use warnings;

# Load other modules required
# POSIX subs
use Cwd;

# Database modules
use DBI;
use DBD::mysql;

# Bio::Perl modules
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;

# Connect to a mysql database
# It determines the host under which it is running and connects to the
# mysql server using the name and password contained in .my.cnf
# Modyfy for the correct host determination if the DB is not running on rusc

# Modified 27.01.2009 to connect to rusc without adding localhost
sub MySQL_DB_Connect {
    my ($database, $host) = @_;

    # This piece of code is only requred if the host is different
    unless ($host) {
	$host='pou';
	print STDERR "WARNING: No host was found, so I've set it to $host\n";
    }

    my $datasource = "DBI:mysql:$database;host=$host";
    my $dbh;
    my $cnf_file='.my.cnf';

    # Check for .my.cnf
    if (-e "$ENV{HOME}/$cnf_file") {
	$datasource .= ";mysql_read_default_file=$ENV{HOME}/$cnf_file";
	$dbh = DBI->connect($datasource, undef, undef, {RaiseError => 1});
    } else {
	print STDERR "Unable to find $ENV{HOME}/.my.cnf\n";
    }
    return $dbh;
}

# Check if a table exists in the database
### 4.0 ok
sub check_table_existence {
    my $dbh=shift;
    my $table=shift;

    print STDERR "Checking database for $table...";

    my ($query,$sth);
    $query ='SELECT count(*) ';
    $query.="FROM $table";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    my $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    my $present=0;
    
    if (@{$results}) {
	# Print the table location for the junctions of this experiment
	print STDERR "Present\n";
	$present=1;
    } else {
	print STDERR "Absent\n";
    }
    return($present);
}

sub check_field_existence {
    my $dbh=shift;
    my $table=shift;
    my $column=shift;

    unless ($column) {
	die "No column name provided\n";
    }

    print STDERR "Checking database for $column in $table...";

    my ($query,$sth);
    $query ='SELECT count(*) ';
    $query.="FROM $table";

    $sth = $dbh->column_info(undef,undef,$table,$column);

    my $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    my $present=0;
    
    if (@{$results}) {
	# Print the table location for the junctions of this experiment
	print STDERR "Present\n";
	$present=1;
    } else {
	print STDERR "Absent\n";
    }
    return($present);
}

# Get a small subroutine to run system commands
sub run_system_command {
    my $command=shift;
    my $log_fh=shift;

    if ($log_fh) {
	print $log_fh "Executing: $command\n";
	print $log_fh `$command`,"\n";
    } else {
	$log_fh=*STDERR;
	print $log_fh "Executing: $command\n";
	system($command);
    }
    sleep(1);
}

# Get the md5sum of a file from the absolute path and return the file name
# without the initial path and the md5sum
sub get_md5sum {
    my $filepath=shift;
    my $filename=$filepath;
    $filename=~s/.+\///;

    my $sum=`md5sum $filepath`;
    my ($md5sum)=split(' ',$sum);
    
    return($filename,$md5sum);
}

# Read or write from a file (gzipped or not)
### OK
sub get_fh {
    my $filename=shift;
    my $write=shift;

    chomp($filename);
    my $openstring;
    my $fh;
    
    if ($write) {
	if ($filename=~/.gz$/) {
	    warn "Piping $filename through gzip\n";
	    # Use gzip -7 as the increase in compression between 7 an 9 is quite
	    # small and the extra time invested quite large
	    $openstring="| gzip -7 -c > $filename";
	} else {
	    $openstring=">".$filename;
	}
    } else {
	if ($filename=~/.gz$/) {
	    warn "Opening gzipped file $filename\n";
	    $openstring="gunzip -c $filename |";
	} else {
	    $openstring=$filename;
	}
    }
    
    open($fh,$openstring) ||
	die "Unable to open $filename: $!,$?\n";

    return($fh);
}

# Get a file handle for logging purposes
### OK
sub get_log_fh {
    my $logfn=shift;
    my $debug=shift;

    my $log_fh;
    if ($debug) {
	# Redirect the log to a filehandle
	print STDERR "WARNING: Debugging mode\n";
	$log_fh=*STDERR;
	# Create a bogus log file in order to avoid problems downstream if the
	# file does not exist
	my $command="touch $logfn";
	print STDERR "Executing: $command\n";
	system($command);
    } else {
	$log_fh=get_fh($logfn,1);
    }
    return($log_fh);
}

# Check if a file exists and is readable by the user
sub check_file_existence {
    my $file=shift;
    my $logfh=shift;

    if (-r $file) {
	if ($logfh) {
	    print $logfh "$file is present and readable\n";
	} else {
	    print STDERR "$file is present and readable\n";
	}
    } else {
	die "$file is not readable\n";
    }
}

# Check if the fasta file identifiers contain any spaces and warn if thise is
# the case
sub check_fasta_file {
    my $file=shift;
    my $fileok=1;

    print STDERR "Checking the genome file...";
    my $in  = Bio::SeqIO->new(-file => $file ,
			      -format => 'Fasta');

    while ( my $seq = $in->next_seq() ) {
	my $seqid=$seq->display_id();
	my $desc=$seq->desc();
	if ($desc=~/[=~:;]/o) {
	    print STDERR "Presence of special characters in the header may cause parsing problems after mapping\n";
	    print STDERR "Check: $desc\n";
	    $fileok=0;
	} elsif ($seqid=~/^\s*\w+\s*$/o) {
	    next;
	} else {
	    print STDERR "WARNING: Spacing in the fasta identifiers may cause problems: $seqid\n";
	}
    }
    $in->close();
    if ($fileok) {
	print STDERR "Fasta file seems fine\n";
    }
    return($fileok);
}

# Check if a file conformas to the UCSC specifications of GTF
sub check_gff_file {
    my $file=shift;

    print STDERR "Checking annotation file...";
    my %features;
    my $infh=get_fh($file);
    while (my $line=<$infh>) {
	if ($line=~/^#/o) {
	    next;
	}
	my %line=%{parse_gff_line($line)};
	$features{$line{'type'}}++;
    }

    # Issue some warnings if the file looks strange
    if ($features{'exon'} < 1) {
	die "No exons found in $infh\n";
    } elsif ($features{'exon'} &&
	     $features{'exon'} < 10000) {
	print STDERR "WARNING: Only $features{'exon'} exons found in $file, maybe it is incomplete\n";
    }
    if ($features{'gene'} &&
	$features{'gene'} < 1) {
	print STDERR "WARNING: No genes found in $file\n";
    }
    if ($features{'transcript'} &&
	$features{'transcript'} < 1) {
	print STDERR "WARNING: No transcripts found in $file\n";
    }
    print STDERR "Annotation seems fine\n";

    close($infh);
}

# This should get all the distinct values present in a certain column of a
# MySQL table. The database handle as well as the column and table names must
# be provided
sub get_distinct_col_value_from_tab {
    my $dbh=shift;
    my $column=shift;
    my $table=shift;

    my @values;
    my ($query,$sth);
    # Get the samples

    $query ="SELECT DISTINCT $column ";
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute(); 

    while (my ($value)=$sth->fetchrow_array()) {
	push @values, $value;
    }

    return(@values);
}

# get a reference to a list of reow one to a list of colimns and one to  data
# frame organized as a hash of hasses containing some values and buid a
# heatmap using R
sub plot_heatmap {
    my $columns=shift;
    my $rows=shift;
    my $data=shift;

    
    # Print columns
    my $columnsfn="$$.columns.txt";
    my $colfh=get_fh($columnsfn,1);
    foreach my $col (@{$columns}) {
	my $col_print=$col;
	$col_print=~s/_(gene|exon)_RPKM//;
	$col_print=~s/Stranded/St/;
	$col_print=~s/[Tt]rimmed/Tr/;
	$col_print=~s/Ginge?r?a?s?/Ging/;
	$col_print=~s/AAXX//;
	print $colfh $col_print,"\n";
    }
    close($colfh);

    # Print data and rows
    my $datafn="$$.data.txt";
    my $rowsfn="$$.rows.txt";
    my $datafh=get_fh($datafn,1);
    my $rowfh=get_fh($rowsfn,1);
    foreach my $gene (sort keys %{$rows}) {
	my @values;
	my $non_zero=0;
	foreach my $sample (@{$columns}) {
	    if (exists $data->{$gene} &&
		(exists $data->{$gene}{$sample})) {
		my $datum=0;
		if ($data->{$gene}{$sample}) {
		    $datum=$data->{$gene}{$sample};
		    $non_zero=1;
		}
		push @values, $datum; 
	    } else {
		push @values, 0;
	    }
	}
	if ($non_zero) {
	    print $datafh join("\t",
			       @values),"\n";
	    print $rowfh $gene,"\n";
	}
    }
    close($datafh);
    close($rowfh);

    # Build R command file
    my $r_file="$$.commands.r";
    my $rfh=get_fh($r_file,1);
    my $r_string;
    $r_string ="Data<-read.table(\"$$.data.txt\")\n";
    $r_string.="Rows<-read.table(\"$$.rows.txt\")\n";
    $r_string.="Cols<-read.table(\"$$.columns.txt\")\n";
    $r_string.="library(\"gplots\")\n";
    $r_string.="Data.mat<-as.matrix(Data)\n";
    $r_string.="postscript(\"heatmap.ps\")\n";	
    $r_string.="hm<-heatmap.2(Data.mat,labRow=Rows\$V1,labCol=Cols\$V1,symm=F,scale='none',density.info='none',col=topo.colors(100),trace='none',margins=c(15,10),dendrogram='column')\n";
    $r_string.="dev.off()\n";

    # Write the order of the rows from bottom to top
    $r_string.="write.table(hm\$rowInd,\"$$.row.order.txt\",col.names=F,row.names=F)\n";
    $r_string.="write.table(hm\$colInd,\"$$.col.order.txt\",col.names=F,row.names=F)\n";
    $r_string.="q()\n";

    print $rfh $r_string;
    close($rfh);

    # Execute R command
    my $command="R --vanilla --slave < $r_file";
    system($command);

    # print the order in which the genes should be examined.
    my $orderfn="$$.row.order.txt";
    my $orderfh=get_fh($orderfn);
    my @roworder;
    while (my $row=<$orderfh>) {
	chomp($row);
	$row--;
	push @roworder,$row;
    }

    # get the order in which the samples are shown in the figure
    my $colorderfn="$$.col.order.txt";
    my $colorderfh=get_fh($colorderfn);
    my @colorder;
    while (my $col=<$colorderfh>) {
	chomp($col);
	$col--;
	push @colorder,$columns->[$col];
    }
    @{$columns}=@colorder;

    # clean up
    $command="rm $r_file $columnsfn $rowsfn $datafn $orderfn $colorderfn";
    system($command);

    return(@roworder);
}



### Process the annotation to extract the different parts of information we are
# interested in
# This will get all the exons from the annotation file and also take a gene
# file where it will add the exon for each of the genes
sub get_exon_list_from_gtf {
    my $exonfile=shift;
    my $genes=shift;
    my $repeated_exons='rep.exons.txt';

    my %exons;
    my %remove;

    %exons=%{get_annotation_from_gtf($exonfile,
				     '',
				     'exons')};

    foreach my $exon (keys %exons) {
	    my $gene_id=$exons{$exon};
	    $genes->{$gene_id}->{$exon}=1;
    }

    print STDERR "done\n";

    my $count=keys %exons;
    print STDERR $count,"\tExons obtained\n";
    $count=keys %{$genes};
    print STDERR $count,"\tGenes\n";

    return(\%exons);
}

# Extract from a gtf annotation all the gene objects using the parse_gff sub
### OK
sub get_annotation_from_gtf {
    my $file=shift;
    my $log_fh=shift;
    my $subset=shift || '';

    my $fh=get_fh($file);
    my %genes;
    my %trans;
    my %exons;
    my %remove_exons;
    my %excluded;
    my $count=0;
    my %trans_count;

    # If we have no $log_fh redirect to STDERR
    unless ($log_fh) {
	$log_fh=*STDERR;
    }

    # Check the subset
    if ($subset eq 'exons') {
	print $log_fh "extracting only exons from $file\n";
    } else {
	print $log_fh "extracting genes, transcripts and exons from $file\n";
    }

    print $log_fh "Reading $file\n";
    print $log_fh "WARNING: Skipping entries located in chr random, chr U(nknown), haplotypes and EnsEMBL assembly exceptions\n";

    while (my $line=<$fh>) {
	$count++;
	chomp($line);

	# Skip possible comments
	if ($line=~/^#/o) {
	    next;
	}

	# use the parse_gff  subroutine to parse the lines
	my ($chr,$type,$start,$end,$strand,$frame,$info,$source);
	my $line=parse_gff_line($line);
	$chr=$line->{'chr'};
	$type=$line->{'type'};
	$start=$line->{'start'};
	$end=$line->{'end'};
	$frame=$line->{'frame'};

	# Skip entries we are not interested in
	# Skip non-exon entries
	unless ($type=~/^(exon|transcript|gene)$/o) {
	    next;
	}

	# Skip random and haplotype chromosomes
	# We will not skip the random ones, as there are some species with
	# annotated and expressed genes in these chromosomes, particularly
	# those where the genome assembly is not very good
#	if ($chr=~/random/io) {
#	    next;
#	}
	if ($chr=~/hap/o) {
	    next;
	} elsif ($chr=~/^chrU/o) {
	    next;
	} elsif ($chr=~/^Un\./o) {
	    # This is for EnsEMBL cow
	    next;
	} elsif ($chr=~/^(chr)?HSCHR/o) {
	    next;
	} elsif ($chr=~/^AAFC03011182/o) {
	    # EnsEMBL cow
	    next;
	}

	### TO DO Fix some naming issues that may occurr
	# If the chromosomes are not named as chr in the file name them so this
	# may cause some problems if we are looking at contigs etc... but
	# it should only activate if the chromosomens are named as the humans
	# but with no chr
#	if ($chr!~/^(chr|contig|scaffold|supercontig)/) {
#	    $chr=~s/^/chr/;
#	}
	# To prevent problems with the naming of the chromosomes we will change
	# the chrMT to chrM
#	$chr=~s/chrMT/chrM/;

	# Check the strand
	if ($line->{'strand'} eq '+') {
	    $strand=1;
	} elsif ($line->{'strand'} eq '-') {
	    $strand=-1;
	} else {
	    warn "Unknown strand $strand\n";
	}

	# Get the gene_id and transcript_id info
	my ($gene_id,$trans_id);
	$gene_id=$line->{'feature'}{'gene_id'};
	unless ($gene_id) {
	    $gene_id=$line->{'feature'}{'ID'};
	}
	$trans_id=$line->{'feature'}{'transcript_id'};

	# Complain if there is anything missing
	unless ($gene_id) {
	    print $log_fh "No gene_id found in:\n",$line,"\n";
	    next;
	}
	unless ($trans_id) {
	    if ($type ne 'gene') {
		print $log_fh "No transcript_id for $gene_id\n";
		next;
	    }
	}
	unless ($start && $end && $strand) {
	    warn "Problem with $info. Excluding $gene_id\n";
	    $excluded{$gene_id}=1;
	}

	# If we have the required information proceed to build the gene,
	# transcript and exon objects
	if ($type eq 'gene') {
	    # Get the gene feature
	    if (exists $genes{$gene_id}) {
		# Complain if two genes with the same ID are annotated on
		# different chromosomes
		unless ($genes{$gene_id}{'chr'} eq $chr) {
		    warn $gene_id ," is present both in ",$genes{$gene_id}{'chr'}," and in $chr so it will be excluded as ambiguous\n";
		    $excluded{$gene_id}=1;
		}
		# Also complain if the script enters this loop, as in any case
		# it means the gene entry is repeated
		# This shouldn't happen, but does in the Drosophila annotation
		# at least. I think the genes affected are some for which
		# structure is not completely resolved
		warn "Gene $gene_id may be a repeated entry (or perhaps the gtf file is unsorted and exons appear before their corresponding gene)\n";
	    } else {
		my $gene_feat = Bio::SeqFeature::Gene::GeneStructure->new(
		    -start        => $start,
		    -end          => $end,
		    -strand       => $strand,
		    -primary      => 'gene', # -primary_tag is a synonym
		    -display_name => $gene_id);
		$genes{$gene_id}{'gene'}=$gene_feat;
		$genes{$gene_id}{'chr'}=$chr;
		$genes{$gene_id}{'start'}=$start;
		$genes{$gene_id}{'end'}=$end;
		$genes{$gene_id}{'strand'}=$line->{'strand'};
		$genes{$gene_id}{'type'}=$type;
		$genes{$gene_id}{'feature'}=$line->{'feature'};
	    }
	} elsif ($type =~ /^transcript$/o) {
	    # Get the transcript feature
	    # If the gene object does not exist build it on the fly. This could
	    # happen with unordered files or files with only the exons
	    unless ($genes{$gene_id}{'gene'}) {
		my $gene_feat = Bio::SeqFeature::Gene::GeneStructure->new(
		    -start        => $start,
		    -end          => $end,
		    -strand       => $strand,
		    -primary      => 'gene', # -primary_tag is a synonym
		    -display_name => $gene_id);
		$genes{$gene_id}{'gene'}=$gene_feat;
		$genes{$gene_id}{'chr'}=$chr;
	    }
	    unless (exists ($genes{$gene_id}{'transcripts'}{$trans_id})) {
		my $trans_feat = Bio::SeqFeature::Gene::Transcript->new(
		    -start        => $start,
		    -end          => $end,
		    -strand       => $strand,
		    -primary      => 'transcript', # -primary_tag is a synonym
		    -display_name => $trans_id);
		$genes{$gene_id}{'transcripts'}{$trans_id}=$trans_feat;
	    }
	    # Add the transcript information if required
	    if ($subset eq 'transcript') {
		$trans{$trans_id}{'gene_id'}=$gene_id;
		$trans{$trans_id}{'transcript_id'}=$trans_id;
		$trans{$trans_id}{'chr'}=$chr;
		$trans{$trans_id}{'type'}=$type;
		$trans{$trans_id}{'feature'}=$line->{'feature'};
	    }
	} elsif ($type =~ /exon/o) {
	    my $exon_id=join('_',
			     $chr,
			     $start,
			     $end,
			     $strand);
	    # It may be useful to have the gene info in this hash, but we will
	    # remove those that are assigned to multiple genes
	    if ($exons{$exon_id} && 
		($exons{$exon_id} ne $gene_id)) {
		$remove_exons{$exon_id}=1;
	    } else {
		$exons{$exon_id}=$gene_id;
	    }

	    # If we only need the basic exon list skip transcript and gene
	    # building
	    if ($subset eq 'exons') {
		next;
	    }

	    # Get the exon feature
	    # If the gene object does not exist build it on the fly
	    unless ($genes{$gene_id}{'gene'}) {
		my $gene_feat = Bio::SeqFeature::Gene::GeneStructure->new(
		    -start        => $start,
		    -end          => $end,
		    -strand       => $strand,
		    -primary      => 'gene', # -primary_tag is a synonym
		    -display_name => $gene_id);
		$genes{$gene_id}{'gene'}=$gene_feat;
		$genes{$gene_id}{'chr'}=$chr;
		$genes{$gene_id}{'type'}='gene';
	    }
	    # If the transcript object does not exist build it on the fly
	    unless ($genes{$gene_id}{'transcripts'}{$trans_id}) {
		my $trans_feat = Bio::SeqFeature::Gene::Transcript->new(
		    -start        => $start,
		    -end          => $end,
		    -strand       => $strand,
		    -primary      => 'transcript', # -primary_tag is a synonym
		    -display_name => $trans_id);
		$genes{$gene_id}{'transcripts'}{$trans_id}=$trans_feat;
		$trans{$trans_id}{'gene_id'}=$gene_id;
		$trans{$trans_id}{'transcript_id'}=$trans_id;
		$trans{$trans_id}{'chr'}=$chr;
		$trans{$trans_id}{'type'}='transcript';
	    }
	    my $exon_feat = Bio::SeqFeature::Gene::Exon->new(
		-start        => $start,
		-end          => $end,
		-strand       => $strand,
		-primary      => 'exon', # -primary_tag is a synonym
		-display_name => $exon_id);
	    if (($frame ne '.') &&
		($frame >= 0)) {
		$exon_feat->frame($frame);
		$exon_feat->is_coding(1);
	    }
	    $genes{$gene_id}{'transcripts'}{$trans_id}->add_exon($exon_feat);
	}
    }
    close($fh);

    # If we only need the exons return this
    if ($subset eq 'exons') {
	$count=keys %exons;
	print STDERR $count, "\tExon entries obtained\n";
	# Remove those exons that map to multiple genes and print to a file
	# if the file does not exist
	my $repeated_exons=$file;
	$repeated_exons=~s/.gtf$/.rep.exons.gtf/;
	$repeated_exons=~s/.*\///;
	if (-r $repeated_exons) {
	    print STDERR $repeated_exons," is already present. Skipping...\n";
	} else {
	    my $repeatfh=get_fh($repeated_exons,1);
	    foreach my $exon (keys %remove_exons) {
		print $repeatfh join("\t",
				     $exons{$exon},
				     $exon),"\n";
	    }
	    close($repeatfh);
	}
	# Do not remove the exons
	foreach my $exon (keys %remove_exons) {
#	    delete $exons{$exon};
	}
	$count=keys %remove_exons;
	print STDERR $count,"\tExons mapping to multiple genes removed\n";
	return(\%exons);
    } elsif ($subset eq 'gene') {
	$count=keys %genes;
	print STDERR $count, "\tGene entries obtained\n";
	return(\%genes);
    } elsif ($subset eq 'transcript') {
	$count=keys %trans;
	print STDERR $count, "\tTranscript entries obtained\n";
	return(\%trans);
    } elsif ($subset) {
	print STDERR "UNknown subset $subset requested. Ignoring...\n";
    }

    # If we want the full annotation continue
    # Add all the transcripts to their corresponding gene feature
    foreach my $gene_id (keys %genes) {
	foreach my $trans_id (keys %{$genes{$gene_id}{'transcripts'}}) {
	    $genes{$gene_id}{'gene'}->add_transcript($genes{$gene_id}{'transcripts'}{$trans_id});
	}
	delete $genes{$gene_id}{'feature'}{'transcript_id'};
    }
    print $log_fh $count,"\tLines read\n";

    # Remove those genes that have been flagged as ambiguous by any of the
    # previous checks
    print $log_fh "Removing ambiguous sequences\n";
    foreach my $id (keys %excluded) {
	delete $genes{$id};
    }
    $count=keys %excluded;
    print $log_fh $count,"\tEntries excluded because of ambiguous genomic location\n";

    $count=keys %genes;
    print $log_fh $count, "\tGene entries obtained\n";

    return(\%genes);
}

# Build all annotated junctions from a hash with gene objects
### OK
sub build_annotated_junctions {
    my $exons=shift;
    my $log_fh=shift;
    my %annotated_junctions;

    # If we have no $log_fh redirect to STDERR
    unless ($log_fh) {
	$log_fh=*STDERR;
    }

    print $log_fh "Extracting all annotated junctions\n";
    foreach my $gene (keys %{$exons}) {
	foreach my $trans ($exons->{$gene}->{'gene'}->transcripts()) {
	    my $strand=$trans->strand();
	    my @exons=$trans->exons_ordered();
	    for(my $i=0;$i<@exons - 1; $i++) {

		my $exon1=$exons[$i];
		my $exon2=$exons[$i + 1];

		my $exon1_id=$exon1->display_name();
		my $exon2_id=$exon2->display_name();

		my $junction_id=$exon1_id.'_splice_'.$exon2_id;

		if ($strand == -1) {
		    $junction_id=$exon2_id.'_splice_'.$exon1_id;
		}
		$annotated_junctions{$junction_id}=$gene;
	    }
	}
    }
    my $count=keys %annotated_junctions;
    print $log_fh $count,"\tAnnotated junctions retrieved\n";
    return(\%annotated_junctions);
}

# This sub will take a dbh & a table (one of the mapping results tables) and it
# will return a sub that will take as input a lane name and return the
# corresponding mapping file
sub get_files_from_table_sub {
    my $dbh=shift;
    my $table=shift;

    # First check if the table exists and if not die graciously
    my $present=check_table_existence($dbh,
				      $table);

    unless ($present) {
	die "Required table $table does not seem to be present in the database...\nI'm quitting\n";
    }

    my %cache;

    my ($query,$sth,$count);

    $query ='SELECT filename ';
    $query.="FROM $table ";
    $query.='WHERE LaneName = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $lane=shift;

	unless ($cache{$lane}) {
	    $count=$sth->execute($lane);

	    if ($count != 1) {
		die "No file in $table corresponds to $lane\n";
	    } else {
		my ($filename)=$sth->fetchrow_array();
		$cache{$lane}=$filename;
	    }
	}
	return($cache{$lane});
    };
    return($subroutine);
}

sub print_gff {
    my $outfh=shift;
    my $readfile=shift;
    my $read_id=shift;
    my $read=shift;
    my $type=shift || '.';

    # Remove preceding directories
    $readfile=~s/.*\///o;

    my $string= join(' ',
		     'read_id','"'.$read_id.'";',
		     'mismatches','"'.$read->[4].'";',
		     'qualities','"'.$read->[5].'";',
		     'matches','"'.$read->[6].'";');


    print $outfh join("\t",
		      $read->[0],
		      $readfile,
		      $type,
		      $read->[1],
		      $read->[2],
		      '.',
		      $read->[3],
		      '.',
		      $string),"\n";
}

# This will take the line returned by the parse_gff_line subroutine
sub format_as_gff {
    my $line=shift;

    my $chr=$line->{'chr'};
    my $start=$line->{'start'};
    my $end=$line->{'end'};
    my $strand=$line->{'strand'};
    my $frame=$line->{'frame'};

    my $source=$line->{'source'};
    my $type=$line->{'type'};
    my $score=$line->{'score'};

    my @features;
    foreach my $feature (keys %{$line->{'feature'}}) {
	push @features, $feature.' "'.$line->{'feature'}{$feature}.'"';
    }
    my $features=join('; ',@features);

    return([$chr,
	    $source,$type,
	    $start,$end,
	    $score,
	    $strand,$frame,$features])
}
    
sub cluster_gff {
    my $infh=shift;
    my $clustername=shift;
    my $threshold=shift || 1;

    my $tmpfn=$clustername;
    my $tmpfh=get_fh($tmpfn,1);

    my $chr;
    my $start;
    my $end;
    my $incl=0;
    my $strand;
    my $type='cl';


    while (my $line=<$infh>) {
	my %line=%{parse_gff_line($line)};
	if ($chr) {
	    if ($chr eq $line{'chr'}) {
		if (feat_overlap($start,$end,
				 $line{'start'},$line{'end'})) {
		    if ($line{'end'} > $end) {
			$end=$line{'end'};
		    }
		    unless ($strand eq $line{'strand'}) {
			$strand='.';
		    }
		    $incl++;
		} else {
		    if ($incl >= $threshold) {
			print $tmpfh join("\t",
					  $chr,
					  'cluster',
					  $type,
					  $start,
					  $end,
					  '.',
					  $strand,
					  '.',
					  "included \"$incl\";"),"\n";
		    }
		    # Start new feature on same chr
		    $start=$line{'start'};
		    $end=$line{'end'};
		    $strand=$line{'strand'};
		    $incl=1;
		}
	    } else {
		# print line
		if ($incl > $threshold) {
		    print $tmpfh join("\t",
				      $chr,
				      'cluster',
				      $type,
				      $start,
				      $end,
				      '.',
				      $strand,
				      '.',
				      "included \"$incl\";"),"\n";
		}
		# Start new feature
		$chr=$line{'chr'};
		$start=$line{'start'};
		$end=$line{'end'};
		$strand=$line{'strand'};
		$incl=1;	
	    }
	} else {
	    #Start first feature
	    $chr=$line{'chr'};
	    $start=$line{'start'};
	    $end=$line{'end'};
	    $strand=$line{'strand'};
	    $incl=1;
	}
    }
    # print last line
    if ($incl > $threshold) {
	print $tmpfh join("\t",
			  $chr,
			  'cluster',
			  $type,
			  $start,
			  $end,
			  '.',
			  $strand,
			  '.',
			  "included \"$incl\";"),"\n";
    }
    close($tmpfh);
}

sub feat_overlap {
    my $begin1=shift;
    my $end1=shift;
    my $begin2=shift;
    my $end2=shift;

    if (($end1 >= $begin2) &&
	($begin1 <= $end2)) {
	return(1);
    } else {
	return(0);
    }
}

# This subroutine will take the gtf file from the junction mappuings and it will
# transform these mappings into genomic coordinates
sub get_gff_from_junc_id {
    my $infiles=shift;
    my $outfiles=shift;

    foreach my $file (@{$infiles}) {
	my $infh=get_fh($file);
	my $outfile=$file;
	$outfile=~s/gz$/parsed.gz/;
	my $outfh=get_fh($outfile,1);
	while (my $line=<$infh>) {
	    my %coords=%{parse_gff_line($line)};
	    my ($chr1,$start1,$end1,$strand1,$splice,
		$chr2,$start2,$end2,$strand2)=split('_',$coords{'chr'});
	    my $mappedlength=$coords{'length'};

	    unless($strand1 eq $strand2) {
		warn "Conflicting strands: ".$coords{'chr'}."\n";
	    }

	    unless($chr1 eq $chr2) {
		die "Conflicting Chromosomes: ".$coords{'chr'}."\n";
	    }

	    # Get the length of each junction half. This should be the end minus
	    # the start plus 1, and also the total length
	    my $half_length1=$end1 - $start1 + 1;
	    my $half_length2=$end2 - $start2 + 1;
	    my $total_length=$half_length1 + $half_length2;

	    # Using the length of each fragment calculate where the map starts
	    # and ends
	    # Adjust start1
	    my ($rstart1,$rend1,$rstart2,$rend2,$roffset)=(0,0,0,0,0);
	    if ($strand1 == 1) {
		if ($coords{'end'} <= $half_length1) {
		    # The read is completely included in the first fragment
		    $rstart1=$start1 + $coords{'start'} - 1;
		    $rend1=$start1 + $coords{'end'} - 1;
		} elsif ($coords{'start'} > $half_length1) {
		    # The read is completely included in the second fragment
		    $rstart1=$start2 + $coords{'start'} - $half_length1 - 1;
		    $rend1=$start2 + $coords{'end'} - $half_length1 - 1;
		} elsif (($coords{'start'} <= $half_length1) &&
			 ($coords{'end'} > $half_length1)) {
		    # The read spans the junction
		    $rstart1=$start1 + $coords{'start'} - 1;
		    $rend1=$end1;
		    $rstart2=$start2;
		    $rend2=$start2 + $coords{'end'} - $half_length1 - 1;
		    $roffset=$rend1 - $rstart1 + 1;
		} else {
		    die "Read problem\n";
		}
	    } elsif ($strand1 == -1) {
		# Here we basically want to mirror the coordinates, so we need
		# the total length of the junction and we will substract each
		# coordinate from this. That way we can actually position the
		# reads the same way as we would in the plus strand
		my $alt_start=$total_length - $coords{'end'} + 1;
		my $alt_end=$total_length - $coords{'start'} + 1;
		if ($alt_end <= $half_length1) {
		    # The read is completely included in the first fragment (which
		    # is the lower coordinate one)
		    $rstart1=$start1 + $alt_start - 1;
		    $rend1=$start1 + $alt_end - 1;
		} elsif ($alt_start > $half_length1) {
		    # The read is completely included in the second fragment
		    $rstart1=$start2 + $alt_start - $half_length1 - 1;
		    $rend1=$start2 + $alt_end - $half_length1 - 1;
		} elsif (($alt_start <= $half_length1) &&
			 ($alt_end > $half_length1)) {
		    # The read spans the junction
		    $rstart1=$start1 + $alt_start - 1;
		    $rend1=$end1;
		    $rstart2=$start2;
		    $rend2=$start2 + $alt_end - $half_length1 - 1;
		    $roffset=$rend1 - $rstart1 + 1;
		} else {
		    die "Read problem\n";
		}
	    } else {
		warn "Unknown strand\n";
	    }

	    # Calculate the strand to output, as this should be with regard to
	    # the genome and not to the junction
	    my $outstrand;
	    if ($strand1 == 1) {
		if ($coords{'strand'} eq '+') {
		    $outstrand='+';
		} elsif ($coords{'strand'} eq '-') {
		    $outstrand='-';
		}
	    } elsif ($strand1 == -1) {
		if ($coords{'strand'} eq '+') {
	    $outstrand='-';
		} elsif ($coords{'strand'} eq '-') {
		    $outstrand='+';
		}
	    } else {
		die "Strand problem\n";
	    }

	    if ($start1 > $start2) {
		warn "coordinate_problem\n";
		print STDERR join("\t",
				  $chr1,
				  $start1 - 1,
				  $end2,
				  $coords{'id'},
				  $coords{'strand'}),"\n";;
	    }
    
	    if (($chr1 eq $chr2) &&
		($strand1 eq $strand2)) {
		my $strand=$strand1;

		print $outfh join("\t",
				  $chr1,
				  'junc',
				  'split',
				  $start1,
				  $end1,
				  '.',
				  $strand,
				  '.',
				  "junction_id \"$coords{'chr'}\";"),"\n";
		print $outfh join("\t",
				  $chr2,
				  'junc',
				  'split',
				  $start2,
				  $end2,
				  '.',
				  $strand,
				  '.',
				  "junction_id \"$coords{'chr'}\";"),"\n";
	    } else {
		print STDERR join("\t",
				  $chr1,$start1,$end1,$strand1,
				  $chr2,$start2,$end2,$strand2),"\n";
		return();
	    }
	}
	close($infh);
	close($outfh);
	push @{$outfiles},$outfile;
    }
}

# Read or write from a file (gzipped or not)
sub get_sorted_gff_fh {
    my $filenames=shift;
    my $stagger=shift;
    my $tmpdir=shift;
    my $openstring;
    my $fh;
    
    $openstring='sort -k1,1 -k4,4n -k5,5n ';
    if ($tmpdir) {
	$openstring.="-T $tmpdir ";
    }
    if ($stagger) {
	$openstring.='| uniq ';
    }
    if ($filenames->[0]=~/.gz$/) {
	warn "Opening gzipped files\n";
	$openstring=join(' ',
			 'zcat',
			 @{$filenames},
			 '|').$openstring;
    } else {
	$openstring=join(' ',
			 'cat',
			 @{$filenames},
			 '|').$openstring;
    }
    $openstring.='|';
    print STDERR "executing: $openstring\n";
    
    open($fh,$openstring) ||
	die "Unable to open files: $!,$?\n";

    return($fh);
}

# Parse strictly and return and die if any problems are encountered 
sub parse_gff_line {
    my $line=shift;
    my $log_fh=shift;

    # If no log_fh is provided redirect to STDERR
    unless ($log_fh) {
	$log_fh=*STDERR;
    }

    chomp($line);
    my %line;

    my @line=split("\t",$line);

    unless (@line >=8) {
	die "Insufficient fields in $line\n";
    }

    $line{'chr'}=$line[0];
    $line{'start'}=$line[3];
    $line{'end'}=$line[4];
    $line{'strand'}=$line[6];
    $line{'frame'}=$line[7];

    $line{'source'}=$line[1];
    $line{'type'}=$line[2];
    $line{'score'}=$line[5];
    
    # Parse the info fields
    if ($line[8]) {
	$line[8]=~s/^\s*//o;
	my @info;

	if ($line[8]=~/;/o) {
	    $line[8]=~s/;$//o;
	    @info=split(/; /,$line[8]);
	} elsif ($line[8]=~/^gene_id/o) {
	    # Check is there is only one key-value pair as in this case there
	    # is not necessarily a ';' If this is the case the attribute should
	    # be gene_id
		@info=($line[8]);
	} else {
	    die "Info field does not seem correct in:\n\t$line\n";
	}


	# If the info field is present print it
	if (@info) {
	    for (my $i=0;$i<@info;$i+=1) {
		my ($key,$value)=split(' ',$info[$i],2);

		if ($value) {
		    $value=~s/"//og;
		    # Add something in order to be compatible with the old
		    # version of overlap

		    # Allow for multiple values of the same tag
		    if ((exists $line{'feature'})&&
			(exists $line{'feature'}{$key})) {
			if ($key eq 'gene_id') {
			    die "gene_id is presence more than once in $line\n";
			} elsif ($key eq 'transcript_id') {
			    die "transcript_id is presence more than once in $line\n";
			} else {
			    $line{'feature'}{$key}.=",$value";
			}
		    } else {
			$line{'feature'}{$key}=$value;
		    }
		} else {
		    if ($line{'type'}!~/(exon|transcript|gene)/o) {
			print $log_fh "Unknown key $key in $line{'type'} line $line\n";
		    } else {
			die "No value found for key $key in $line{'type'} line $line\n";
		    }
		}
	    }
	} else {
	    die "Info field is missing from gtf in:\n$line\n";
	}
    }

    return(\%line);
}

sub get_chr_subseq {
    my $genomefile=shift;
    my $tmpdir=shift;
    my $bordering=shift;
    my $reverse=shift;

    my %chromosomes;
    my %cache;

    # Generate the chromosomes
    print STDERR "Generating $genomefile index...";
    # Here we reindex because the script seems not to find the required
    # index if not when the index is in a different dir
    my $db=Bio::DB::Fasta->new($genomefile,
			       -reindex => 1);
    print STDERR "done\n";

    # Build the subroutine for the search
    my $get_seq_from_chrom= sub {
	my $line=shift;

	if (@{$line} == 9) {
	    my $chr=$line->[0];
	    my $start=$line->[3];
	    my $end=$line->[4];
	    if ($bordering) {
		$start-=$bordering;
		$end+=$bordering;
	    }
	    my $strand=$line->[6];
	    my $seq;
	    
	    my $key=join('_',
			 $chr,
			 $start,
			 $end);
	    
	    unless ($cache{$key}) {
		# Get the slice
		my $seq;
		if ($strand eq '-') {
		    $seq=$db->subseq($chr,$end,$start);
		} elsif ($strand eq '+') {
		    $seq=$db->subseq($chr,$start,$end);
		} elsif ($strand eq '.') {
		    if ($reverse) { 
			$seq=$db->subseq($chr,$end,$start);
		    } else {
			$seq=$db->subseq($chr,$start,$end);
		    }
		} else {
		    warn "Problem with strand\n";
		}

		if ($bordering && $bordering > 0) {
		    my $seq1=substr($seq,0,$bordering);
		    my $seq2=substr($seq,-$bordering);
		    $seq=join("\t",
			      $seq1,
			      $seq2);
		}
		
		$cache{$key}=$seq;
		chomp($cache{$key});
	    }
	    return($cache{$key});
	} else {
	    warn "Problem with the arguments\n";
	}
    };
    return($get_seq_from_chrom);
}

### TO DO rewrite to use the parse_gff_line subroutine
sub gff_line_check {
    my $line=shift;
    chomp($line);
    my @parsed;
    my $parsed;

    my @line=split("\t",$line);
    my $type=$line[2];

    unless (@line >=8) {
	return({});
    }

    @parsed=@line[0..7];
    
    # Parse the info fileds
    if ($line[8]) {
	$line[8]=~s/\s+$//;
	$line[8]=~s/;$//;
	my @info1=split('; ',$line[8]);
	my @info;
	foreach my $frag (@info1) {
	    $frag=~s/"//g;
	    my ($key,$val)=split(' ',$frag,2);
	    $val=~s/\s+/_/g;
	    push @info,$key,$val;
	}
	my @features;
	
	for (my $i=0;$i<@info;$i+=2) {
	    unless (defined ($info[$i + 1])) {
		print STDERR "Error at:\t$line\n";
		die "Incorrect number of fields\n";
	    }
	    my $feature='"'.$info[$i + 1].'";';
	    push @features, $info[$i],$feature;
	}
	my $features=join(' ',@features);
	push @parsed,$features;
    }
    $parsed=join("\t",@parsed);

    return($parsed,$type);
}

###
# Here are a series of subroutines used to run sarahs overlap program on files
# regardless of their size
### TO DO
sub get_feature_overlap {
    warn "THIS SHOULD NEVER HAPPEN\n";
    my $features=shift;
    my $annotation=shift;
    my $stranded=shift;
    my $outfile=shift;
    my $tmpdir=shift;
    my $size=1000000;
    my $flags='-v -m -10 -ucsc';

    print STDERR "WARNING: THIS subroutine is slower than the parallel versio\n";

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

}

# This subroutine will split the first instead of the second file into fragments
# before overlapping
### TO DO

# Some statistics subroutines
# Get the average of an array of numbers
sub get_Average {
    my $values=shift;
    my $decimals=shift;
    $decimals=_set_decimals($decimals);
    my $count=@{$values};
    my $total=0;

    foreach my $value (@{$values}) {
	$total+=$value;
    }

    return(sprintf("%.${decimals}f",$total/$count));
}

sub _set_decimals {
    my $decimals=shift;
    # Set the number of decimal places unless it is specified in the sub call
    unless ($decimals) {
	$decimals=2;
    }
    return($decimals);
}

sub get_feat_mappable_length_sub {
    my $dbh=shift;
    my $read_length=shift;
    my $table=shift;;

    my ($query,$sth,$count);
    # Determine ith the necessary table is present
    print STDERR "Checking for the presence of $table...\n";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    
    if (@{$results}) {
	# Print the table location for the junctions of this experiment
	print STDERR join("\t",
			  $table,
			  "Present"),"\n";
    } else {
	# Exit with an error
	warn $table,"\tIs not present...\nI'm dying...\n";
	sleep(rand(10));
	warn "Still dying, I guess...\n";
	sleep(rand(10));
	warn "Just about there, I think, never been here before...\n";
	sleep(rand(10));
	warn "I can see the light, must be close...\n";
	sleep(rand(10));
	die "Oh come on!..You must be kidding!\n";
    }

    # Determine the closest read length in the database
    $query ='SELECT DISTINCT read_length ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();
    
    my $adjusted_length;
    while (my ($length)=$sth->fetchrow_array()) {
	if ($adjusted_length) {
	    if (($read_length - $adjusted_length) > abs($read_length - $length)) {
		$adjusted_length=$length;
	    }
	} else {
	    $adjusted_length=$length;
	}
    }

    unless ($count) {
	die "No entries in $table. Please generate required mappability for $read_length and rerun\n";
    }

    my $diff=abs($adjusted_length - $read_length);
    if ($diff == 0) {
	print STDERR $read_length, "\tPresent in $table\n";
    } elsif ($diff <= 10) {
	print STDERR $read_length, "\tNot present in the table. Using closest length: $adjusted_length\n";
    } else {
	die $read_length, "\tNot present in $table. No similar length present. Please generate mappability for this length and rerun\n";
    }

    # Generate the subroutine
    my %cache;

    $query ='SELECT mappable_length ';
    $query.="FROM $table ";
    $query.='WHERE feat_id = ? AND read_length = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $exon_id=shift;

	unless ($cache{$exon_id}) {
	    $count=$sth->execute($exon_id,$adjusted_length);

	    if ($count != 1) {
		die "No exon in $table corresponds to $exon_id\n";
	    } else {
		($cache{$exon_id})=$sth->fetchrow_array();
	    }
	}
	return($cache{$exon_id});
    };
    return($subroutine);
}

sub get_exon_mappable_length_sub {
    my $dbh=shift;
    my $read_length=shift;
    my $table='exon_mappable_lengths';

    my ($query,$sth,$count);
    # Determine ith the necessary table is present
    print STDERR "Checking for the presence of $table...\n";

    $sth = $dbh->table_info(undef,undef,$table,"TABLE");

    $count=$sth->execute();
    my $results=$sth->fetchall_arrayref();
    
    if (@{$results}) {
	# Print the table location for the junctions of this experiment
	print STDERR join("\t",
			  $table,
			  "Present"),"\n";
    } else {
	# Exit with an error
	warn $table,"\tIs not present...\nI'm dying...\n";
	sleep(rand(10));
	warn "Still dying, I guess...\n";
	sleep(rand(10));
	warn "Just about there, I think, never been here before...\n";
	sleep(rand(10));
	warn "I can see the light, must be close...\n";
	sleep(rand(10));
	die "Oh come on!..You must be kidding!\n";
    }

    # Determine the closest read length in the database
    $query ='SELECT DISTINCT read_length ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();
    
    my $adjusted_length;
    while (my ($length)=$sth->fetchrow_array()) {
	if ($adjusted_length) {
	    if (($read_length - $adjusted_length) > abs($read_length - $length)) {
		$adjusted_length=$length;
	    }
	} else {
	    $adjusted_length=$length;
	}
    }

    unless ($count) {
	die "No entries in $table. Please generate required mappability for $read_length and rerun\n";
    }

    my $diff=abs($adjusted_length - $read_length);
    if ($diff == 0) {
	print STDERR $read_length, "\tPresent in $table\n";
    } elsif ($diff <= 10) {
	print STDERR $read_length, "\tNot present in the table. Using closest length: $adjusted_length\n";
    } else {
	die $read_length, "\tNot present in $table. No similar length present. Please generate mappability for this length and rerun\n";
    }

    # Generate the subroutine
    my %cache;

    $query ='SELECT mappable_length ';
    $query.="FROM $table ";
    $query.='WHERE exon_id = ? AND read_length = ?';
    $sth=$dbh->prepare($query);

    my $subroutine=sub {
	my $exon_id=shift;

	unless ($cache{$exon_id}) {
	    $count=$sth->execute($exon_id,$adjusted_length);

	    if ($count != 1) {
		die "No exon in $table corresponds to $exon_id\n";
	    } else {
		($cache{$exon_id})=$sth->fetchrow_array();
	    }
	}
	return($cache{$exon_id});
    };
    return($subroutine);
}

sub get_exon_sequences {
    my $file=shift;
    my %exons;

    print STDERR "Extracting exon info from $file\n";
    my $infh=Bio::SeqIO->new(-file => $file,
			     -format => 'fasta');

    while (my $seq=$infh->next_seq()) {
	my $sequence=$seq->seq();
	my $seq_id=$seq->display_id();
	if($sequence) {
	    $exons{$seq_id}=$sequence;
	} else {
	    die "No sequence for $seq_id\n";
	}
    }
    $infh->close();
    print STDERR "done\n";
    my $count=keys %exons;
    print STDERR $count,"\tExon sequences obtained\n";

    return(\%exons);
}

sub get_exon_coverage_1000nt {
    my $infn=shift;
    my $features=shift;
    my $readlength=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage,$length,$type)=split("\t",$line);
	$length=$length - $readlength + 1;

	unless ($length > 0) {
	    $length=0;
	    warn "No length for $feature setting to 1 (could be insertions)\n";
	    $length=1;
	}

	# Normalize by length as reads per 1000 nt
	$coverage=$coverage * (1000 / $length);
	if ($coverage < 0) {
	    die "Coverage is below zero ????\n";
	}
	$features->{$feature}+=$coverage;
    }

    print STDERR "done\n";
}

1;
