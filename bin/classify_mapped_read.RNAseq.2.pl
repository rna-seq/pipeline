#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script will take a file with the exon annotation, another with the gene
# annotation and a set of read files.
# It will classify each read into exonic, intronic or intergenic in an exclusive
# manner If a read overlaps an exon it is exonic, if it doesn't but it overlaps
# a gene annotation it is intronic, and if it does not overlap a gene it is
# intergenic.

###
### To do
### Those reads that map to a junction should be included in the exonic
### Parallelize
### recalculate from the merged bam files

use RNAseq_pipeline3 ('get_fh','parse_gff_line','get_feature_overlap_sub');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list');
use POSIX qw(getcwd);

my $prefix;
my $exondir;
my $genomedir;
my $exonsfile;
my $genesfile;
my $tmpdir;
my $stranded;
my $parallel='cluster';
my $paralleltmp;
my $bindir;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$exondir=$options{'EXONDIR'};
$genomedir=$options{'GENOMEDIR'};
$exonsfile=$exondir.'/'.$prefix.'.exon.gtf';
$genesfile=$genomedir.'/'.$prefix.'.gene.gtf';
$tmpdir=$options{'LOCALDIR'};
$stranded=$options{'STRANDED'};
$paralleltmp=$options{'PROJECT'}.'/'.$options{'LOCALPARALLEL'};
$bindir=$options{'BIN'};

# Get the subroutines
*run_overlap=get_feature_overlap_sub($parallel,
				     $paralleltmp,
				     $bindir,
				     '-m 0');

# Check input
unless (-r $exonsfile) {
    die $exonsfile,"\tIs not readable\n";
}
unless(-r $genesfile) {
    die $genesfile,"\tIs not readable\n";
}

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

foreach my $lane (keys %lanes) {
    my $type;
    if (keys %{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (keys %{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    my %read_types;

    my $uniquemaps=$genomedir.'/'.$lane.'.'.$type.'.unique.gtf.gz';
    my $file=$uniquemaps;
    my $temporary;
    print STDERR "Processing $uniquemaps...\n";

    # Decompress the files
    print STDERR "Decompressing $file...";
    unless ($file=~s/.gz//) {
	warn "$file is not compressed\n";
    }
    $file=~s/.*\///;
    $file=$tmpdir.'/'.$file;
    my $command="gunzip $uniquemaps -c > $file";
    system($command);
    $temporary=$file;
    print STDERR "done\n";

    # Run overlap on the exons
    my $flags=' -v -m 0';
    my $exonsout=$file.'.exons';
    $exonsout=~s/.*\///;
    $exonsout=$tmpdir.'/'.$exonsout;
    get_feature_overlap_split1($file,
			       $exonsfile,
			       $stranded,
			       $exonsout,
			       $tmpdir);

    # Run overlap on the genes
    my $genesout=$file.'.genes';
    $genesout=~s/.*\///;
    $genesout=$tmpdir.'/'.$genesout;
    get_feature_overlap_split1($file,
			       $genesfile,
			       $stranded,
			       $genesout,
			       $tmpdir);

    # Classify the reads and count them
    classify_reads($exonsout,
		   $genesout,
		   \%read_types);

    # Keep the overlap files
    $command= "mv $exonsout.overlap.gz $exondir";
    print STDERR "Executing: $command\n";
    system($command);
    $command= "mv $genesout.overlap.gz $genomedir";;
    print STDERR "Executing: $command\n";
    system($command);

    if ($temporary) {
	clean_up2($temporary);
    }

    # Print out the results
    my @results;
    foreach my $type ('exonic','intronic','intergenic','total') {
	$read_types{$type}->[1]=sprintf "%.2f",($read_types{$type}->[0] / 
						$read_types{'total'}->[0]) * 100;
	push @results,@{$read_types{$type}};
    }
    print join("\t",
	       $lane,
	       @results),"\n";
}

exit;

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

sub clean_up2 {
    while (my $file=shift) {
	my $command="rm $file";
	system($command);
    }
}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
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

sub classify_reads {
    my $exonsover=shift;
    my $genesover=shift;
    my $types=shift;

    my $exonsfh=get_fh($exonsover.'.overlap.gz');
    my $genesfh=get_fh($genesover.'.overlap.gz');

    print STDERR "Classifying reads...";

    while (my $line1=<$exonsfh>) {
	my $line2=<$genesfh>;
	if ($line1 && $line2) {
	    $types->{'total'}->[0]++;
	} else {
	    warn "One line is undefined??\n";
	}
	my %line1=%{parse_gff_line($line1)};
	my %line2=%{parse_gff_line($line2)};

	my $exonoverlap=$line1{'feature'}{'ov_feat2:'};
	my $geneoverlap=$line2{'feature'}{'ov_feat2:'};

	# decide on the type of read it is
	if ($exonoverlap > 0) {
	    $types->{'exonic'}->[0]++;
	} elsif ($geneoverlap > 0) {
	    $types->{'intronic'}->[0]++;
	} elsif (($exonoverlap == 0) &&
		 ($geneoverlap == 0)) {
	    $types->{'intergenic'}->[0]++;
	} else {
	    warn "We have a problem...\n";
	}
    }
    close($exonsfh);
    close($genesfh);
    print STDERR "done\n";
}
