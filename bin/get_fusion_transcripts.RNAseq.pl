#!/soft/bin/perl
# DGK

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
# This script should take as input the uniquely mapped reads that map to the
# transcriptome and determine if there are any reads that map between different
# transcript. If so it will classify these into intragenic and intergenic

use RNAseq_pipeline3 qw(get_fh get_log_fh parse_gff_line);
use RNAseq_pipeline_settings3 ('get_gene_from_trans_sub','get_dbh',
			      'read_file_list','read_config_file');
use RNAseq_GEM3 ('parse_gem_line');

# Declare some variables
my $prefix;
my $transdir;
my $gendir;
my $debug;
my %stats;

# Get a log file
my $log_fh=get_log_fh('get_fusion_transcripts.RNAseq.log',
		      $debug);

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$transdir=$options{'TRANSDIR'};
$gendir=$options{'GENOMEDIR'};

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};
my %lane2file=%{get_lane2file(\%files)};

# Get the require subroutines.
# If the dbh is opened here and the file handles after there will be an erro
# that is coaused by the forks in the filehandles closing before the dbh
# connection is closed
*trans2gene=get_gene_from_trans_sub();
my $rpkmtable=$prefix.'_gene_RPKM';
*geneRPKM=get_gene_RPKM_sub($rpkmtable);

# Read and process the transcript raw mapping files
foreach my $lane (keys %lanes) {
    print $log_fh "Getting info for lane: $lane\n";

    # Determine if the reads are paired or not
    my $type;
    if (keys %{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (keys %{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    # Get a file for the output
    my $outfile=$transdir.'/'.$lane.'.'.$type.'.trans.fusions.txt.gz';
    my $outfh=get_fh($outfile,1);

    # If the reads are paired process them using the 4 files
    print $log_fh "Reads identified as paired\n";
    if ($type eq 'paired') {
	# First get the mapping files
	my ($lane1,$lane2)=sort keys %{$lanes{$lane}};
	$lane1=$lane2file{$lane1};
	$lane2=$lane2file{$lane2};
	print STDERR $lane1,"\n";

	my $transfile1=$transdir.'/'.$lane1.'.gem.map';
	if (-r $transfile1) {
	    print STDERR $transfile1,"\tPresent\n";
	} elsif (-r $transfile1.'.gz') {
	    $transfile1.='.gz';
	    print STDERR $transfile1,"\tPresent\n";
	}

	my $genfile1=$gendir.'/'.$lane1.'.gem.map';
	if (-r $genfile1) {
	    print STDERR $genfile1,"\tPresent\n";
	} elsif (-r $genfile1.'.gz') {
	    $genfile1.='.gz';
	    print STDERR $genfile1,"\tPresent\n";
	}

	my $transfile2=$transdir.'/'.$lane2.'.gem.map';
	if (-r $transfile2) {
	    print STDERR $transfile2,"\tPresent\n";
	} elsif (-r $transfile2.'.gz') {
	    $transfile2.='.gz';
	    print STDERR $transfile2,"\tPresent\n";
	}
	my $genfile2=$gendir.'/'.$lane2.'.gem.map';
	if (-r $genfile2) {
	    print STDERR $genfile2,"\tPresent\n";
	} elsif (-r $genfile2.'.gz') {
	    $genfile2.='.gz';
	    print STDERR $genfile2,"\tPresent\n";
	}
	
	if ((-r $transfile1)  && 
	    (-r $genfile1) &&
	    (-r $transfile2)  && 
	    (-r $genfile2)) {
	    print $log_fh "Processing transcript files $transfile1 and $transfile2\n";
	    print $log_fh "Processing genome files $genfile1 and $genfile2\n";
	} else {
#	    warn $lane,"\tFiles missing. Skipping\n";
#	    next;
	    die "Can't find all necessary files\n";
	}

	# Open all the files
	my $transfh1=get_fh($transfile1);
	my $transfh2=get_fh($transfile2);
	my $genfh1=get_fh($genfile1);
	my $genfh2=get_fh($genfile2);

	my $lineNo=0;
	# Read through all the files to determine what is going on
	while (my $trans1line=<$transfh1>) {
	    my $trans2line=<$transfh2>;
	    my $gen1line=<$genfh1>;
	    my $gen2line=<$genfh2>;
	    
	    $lineNo++;
	    # Parse the transcript line
	    my %trans1hit=%{parse_gem_line($trans1line)};
	    my %trans2hit=%{parse_gem_line($trans2line)};

	    # Get the read ID
	    my $readId;
	    my $readid1=$trans1hit{'id'};
	    my $readid2=$trans2hit{'id'};
	    $readid1=~s/\/?p?[12]$//;
	    $readid2=~s/\/?p?[12]$//;

	    # Check the Ids are the same
	    if ($readid1 eq $readid2) {
		$readId=$readid1;
	    } else {
		print $log_fh $trans1line;
		print $log_fh $trans2line;
		print $log_fh $lineNo,"\n";
		die "Read Ids do not match: $readid1 - $readid2\n";
	    }

	    # Get the lenght of the read
	    my $length=$trans1hit{'length'};

	    # Check that the mappings to the transcript are uique in both cases
	    if (($trans1hit{'matches'}=~/^(0:)*1(:0)*$/) &&
		($trans2hit{'matches'}=~/^(0:)*1(:0)*$/)){
	    } else {
		$stats{'not_unique'}++;
		next;
	    }

	    # Check if the reads map to the same transcript
	    # First get the mismatch string
	    my ($transid1,$mismatches1)=split(':',$trans1hit{'hits'});
	    my ($transid2,$mismatches2)=split(':',$trans2hit{'hits'});
	    
	    if ($transid1 eq $transid2) {
		$stats{'eq_trans'}++;
		next;
	    }

	    # If this is not the case we have a candidate. In this case we will
	    # check the string of mismatches, and see if it is the same, if it 
	    # is we will select the hit, if not we will report is as having a
	    # different location.

	    # We will first check the number of mismatches. However this is
	    # not sufficient to discriminate the problematic cases so even
	    # if it is the same we will after check the actuoal substitutions,
	    # however if here we see differences we can save the other check

	    # First get the genomic hit
	    my %gen1hit=%{parse_gem_line($gen1line)};
	    my %gen2hit=%{parse_gem_line($gen2line)};

	    # Get the chromosome and location of the hits
	    my ($chr1,$genmismatches1,$genestart1);
	    my ($chr2,$genmismatches2,$genestart2);
	    my ($transstart1,$transstart2);

	    my $genomeequal=0;
	    ### TO DO
	    # Add the junctions to the comparison, we should check if a read
	    # that does not find a genomic hit finds a junctions hit instead
	    if (($trans1hit{'matches'} eq $gen1hit{'matches'}) &&
		($trans2hit{'matches'} eq $gen2hit{'matches'})) {
		# This is a promissing candidates
		$stats{'equalmatchesboth'}++;

		# First check the mismatch string and get as a side effect the
		# startting position and strand in the transcript
		$mismatches1=~s/^([RF]\d+)//;
		$transstart1=$1;
		$mismatches2=~s/^([RF]\d+)//;
		$transstart2=$1;

		# Get the chromosome and location from the genomic hits
		($chr1,$genmismatches1)=split(':',$gen1hit{'hits'});
		($chr2,$genmismatches2)=split(':',$gen2hit{'hits'});
		$genmismatches1=~s/^([RF]\d+)//;
		$genestart1=$1;
		$genmismatches2=~s/^([RF]\d+)//;
		$genestart2=$1;

		# Compare the mismatch strings
		# Here we are excluding those cases where the last nucleotide 
		# of the the read is a mismatch in bot genome and transcriptome,
		# so probaly the mappings are the same, positions even if not
		# the same sequences (one spliced and one not possibly)
		if (($mismatches1 eq $genmismatches1) &&
		    ($mismatches2 eq $genmismatches2)) {
		    $stats{'equalmismatchesboth'}++;
		    $genomeequal=1;
		} elsif ($mismatches1 eq $genmismatches1) {
		    $stats{'equalTranGensmismatches2'}++;
		    print $log_fh "Different mismatch strings\n";
		    print $log_fh $trans2line;
		    print $log_fh $gen2line,"\n";
		} elsif ($mismatches2 eq $genmismatches2) {
		    $stats{'equalTransGenmismatches1'}++;
		    print $log_fh "Different mismatch strings\n";
		    print $log_fh $trans1line;
		    print $log_fh $gen1line,"\n";
		} else {
		    $stats{'equalTransGenmismatchesBoth'}++;
		    print $log_fh "Different mismatch strings\n";
		    print $log_fh $trans2line;
		    print $log_fh $gen2line;
		    print $log_fh $trans1line;
		    print $log_fh $gen1line,"\n";
		}
	    } elsif ($trans1hit{'matches'} eq $gen1hit{'matches'}) {
		$stats{'diffmismatches2'}++;
	    } elsif ($trans2hit{'matches'} eq $gen2hit{'matches'}) {
		$stats{'diffmismatches1'}++;
	    } else {
		$stats{'noequalmatches'}++;
	    }
	    unless ($genomeequal) {
		next;
	    }

	    # Get the gene Ids for each of the transcripts
	    my $gene1id=trans2gene($transid1);
	    my $gene2id=trans2gene($transid2);

	    # Get the RPKM value for each of the genes
	    my $gene1rpkm=geneRPKM($gene1id,$lane);
	    my $gene2rpkm=geneRPKM($gene2id,$lane);

	    # Classify the candidate into Variant or fusion
	    my $type='Variant';
	    if ($gene1id ne $gene2id) {
		$type='Fusion';
	    }

	    # Once we get here we have a candidate, so we will take it and
	    # Save the pair in a hash indexed by the read
	    my $member1=[$transstart1,
			 $transid1,
			 $chr1,
			 $genestart1,
			 $gene1id,
			 $gene1rpkm];
	    my $member2=[$transstart2,
			 $transid2,
			 $chr2,
			 $genestart2,
			 $gene2id,
			 $gene2rpkm];
	    print $outfh join("\t",
			      $readId,
			      $length,
			      $type,
			      @{$member1},
			      @{$member2},
			      $lane),"\n";
	}
	# Close all the files
	close($transfile1);
	close($transfile2);
	close($genfile1);
	close($genfile2);
    } else {
	print STDERR "Finding fusions from pairs cannot be done for single ends unless I find my cristal ball.\n";
	print STDERR "Sorry for any inconveniences and have a nice day :)\n";
    }
    close($outfh);
}

# print out the stats
foreach my $key (keys %stats) {
    print $log_fh join("\t",
		      $key,
		      $stats{$key}),"\n";
}

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}

sub get_lane2file {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[1]}=$file;
	$lanes{$files->{$file}->[1]}=~s/.fa(stq)*$//;
    }

    return(\%lanes);
}

sub get_gene_RPKM_sub {
    my $table=shift;
    my $dbh=get_dbh();

    my %cache;

    my ($query,$sth);
    $query ='SELECT RPKM ';
    $query.="FROM $table ";
    $query.='WHERE gene_id = ? AND LaneName = ?';
    $sth=$dbh->prepare($query);

    my $subroutine= sub {
	my $gene=shift;
	my $lane=shift;

	unless ($cache{$gene}) {
	    my $count=$sth->execute($gene,$lane);
	    if ($count > 1) {
		die "Too many hits recovered";
	    } else {
		my ($rpkm)=$sth->fetchrow_array();
		$cache{$gene}=$rpkm || 0;
	    }
	}
	return($cache{$gene});
    };
    return($subroutine);
}
