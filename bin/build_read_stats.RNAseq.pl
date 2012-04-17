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
# This script will take as an input a list  of fasta or fastq files and it will
# output a summary of the read statistics in the file

# Input: Read files
# Output three stats files with inforamtion regarding the read qualities

### TO DO
# Add the median to the qualities per position
# Add the avertage to the N's per position
# Add a sub to guess the quality type and rerun the fastq parsing with the
# correct settings if necessary

# Load some modules
use RNAseq_pipeline3 qw(get_fh get_log_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);
use Bio::SeqIO;
use Bio::DB::Sam;
use Tools::Bam qw(process_bam_file);

# Get some options from the configuration file
my $species;
my $project;
my $tmpdir;
my $file_list;
my $prefix;
my $qualities;
my $debug=0;
my $genomefn;

my %options=%{read_config_file()};
$species=$options{'SPECIES'};
$project=$options{'PROJECT'};
$tmpdir=$options{'LOCALDIR'};
$file_list=$options{'FILELIST'};
$prefix=$options{'PREFIX'};
$qualities=$options{'QUALITIES'};
$genomefn=$options{'GENOMESEQ'};

# Get the files we are going to process
my %files=%{read_file_list()};

my %stats;
my %quality_stats;
my %ntpos; # This will contain the number of each of the different nucleotides
           # at a given position

# Get a log file
my $log_fh=get_log_fh('build_read_stats.RNAseq.log',
		      $debug);

foreach my $infile (keys %files) {
    print $log_fh "Processing $infile\n";

    my $lanename=$files{$infile}->[1];

    # Identify the file type
    my $filetype;
    if ($infile=~/.fa(sta)?$/) {
	$filetype='fasta';
    } elsif ($infile=~/.fastq$/) {
	$filetype='fastq';
    } elsif ($infile=~/.bam$/) {
	$filetype='bam';
    } else {
	warn "Unknown filetype for $infile.\nFiletype is guessed from the file ending (fa, fasta or fastq), if the file namings are different please fix them.\n";
	next;
    }
    print $log_fh "$infile identified as $filetype\n";

    my ($good,$bad,$total,$length,$unique);

    # Process tha file
    if ($filetype eq 'fasta') {
	($good,$bad,$total,$length,$unique)=process_fasta_file($infile,
							       \%stats,
							       \%quality_stats,
							       $tmpdir,
							       $log_fh,
							       $lanename);
    } elsif ($filetype eq 'fastq') {
	($good,$bad,$total,$length,$unique)=process_fastq_file($infile,
							       \%stats,
							       \%quality_stats,
							       $tmpdir,
							       $log_fh,
							       $lanename,
							       \%ntpos,
							       $qualities);
    } elsif ($filetype eq 'bam') {
	($good,$bad,$total,$length,$unique)=process_bam_file($infile,
							     \%stats,
							     \%quality_stats,
							     $tmpdir,
							     $log_fh,
							     $lanename,
							     \%ntpos,
							     $qualities,
							     $genomefn);
    } else {
	die "How did I get here???\n";
    }

    $stats{$infile}=[$species,$project,
		     $length, $total,
		     $good,$bad,$unique,$lanename];

}


# Print out Nucleotide distribution
print $log_fh 'Building the nt distribution data...';
my $ntposfh=get_fh("${prefix}_ambiguous.txt",1);
#my $max=0;
foreach my $lanename (sort keys %ntpos) {
    foreach my $pos (sort {$a <=> $b} keys %{$ntpos{$lanename}}) {
#	if ($ntpos{$lanename}{$pos}[0] &&
#	    $ntpos{$lanename}{$pos}[0] > $max) {
#	    $max = $ntpos{$lanename}{$pos}[0];
#	} 

	print $ntposfh join("\t",
			    $lanename,
			    $pos,
			    @{$ntpos{$lanename}{$pos}}),"\n";
    }
}
close ($ntposfh);
print $log_fh "done\n";

# Print out quality pos data
print $log_fh 'Building the quality pos data...';
my $qualityposfh=get_fh("${prefix}_qualitiespos.txt",1);
foreach my $lanename (sort keys %quality_stats) {
    foreach my $pos (sort {$a <=> $b} keys %{$quality_stats{$lanename}}) {
	print $qualityposfh join("\t",
				 $lanename,
				 $pos,
				 $quality_stats{$lanename}{$pos}),"\n";
    }
}
close ($qualityposfh);
print $log_fh "done\n";

# Adjust the max a bit so there is some space under the legend
#$max*=1.2;

# Print out the read summaries
print $log_fh 'Building the summary file...';
my $summaryfh=get_fh("${prefix}_read_stats.txt",1);
# Sort on the laneid in order to be consistent with the other sortings
foreach my $file (sort {$stats{$a}->[7] cmp $stats{$b}->[7]} keys %stats) {
    my $filename=$file;
    $filename=~s/.*\///;
    print $summaryfh join("\t",
			  $filename,
			  @{$stats{$file}}),"\n";
}
close($summaryfh);
print $log_fh "done\n";

# Skip the plotting of the graphs in R
#print $log_fh 'Plotting the pictures...';
#plot_graphs_R(\%stats,
#	      $max,
#	      $prefix);
#print $log_fh "done\n";

close($log_fh);

exit;

### TO DO
# Fix this script so the lane names are obtained in a better way and there is no
# chance of them being incorrectly ordered. I think this is done, although
# it could be more elegant. This is because all the files are sorted by lanename
# so the order should always be the same
sub get_lanes {
    my $file=shift;
    my $infh=get_fh("$file.txt");
    my @lanes;
    while (my $line=<$infh>) {
	chomp($line);
	my @line=split("\t",$line);
	my $name=$line[8];
	push @lanes,"'".$name."'";
    }
    my $name_str=join(',',
		      @lanes);
    return($name_str);
}

sub plot_graphs_R {
    my $stats=shift;
    my $max=shift;
    my $prefix=shift;
    my $length=0;
    my $names=get_lanes("${prefix}_read_stats");
    my @names=split(',',$names);
    my $lanes=@names;
    my $name_size=1 - (0.05 * $lanes);
    if ($name_size < 0.2) {
	$name_size=0.2;
    }
    $lanes*=5;
    my %figures=('ps' => "postscript(\"${prefix}_read_stats.ps\")\n");

    # Get the maximum read length
    foreach my $file (keys %{$stats}) {
	if ($length < $stats->{$file}->[2]) {
	    $length=$stats->{$file}->[2];
	}
    }

    # Build the R command file
    my $execution_file="execution.$$.r";
    my $exec_fh=get_fh($execution_file,1);
    my $r_string;

    # Read the data
    $r_string.="qualities<-read.table(\"${prefix}_qualitiespos.txt\",sep=\"\t\")\n";
    $r_string.="ambiguous<-read.table(\"${prefix}_ambiguous.txt\",sep=\"\t\")\n";
    $r_string.="summary<-read.table(\"${prefix}_read_stats.txt\",sep=\"\t\")\n";

    # Calculate the qualities part
    $r_string.="qualsfile<-split(qualities,qualities\$V1)\n";

    # Calculate the barplots
    $r_string.="matrix<-as.matrix(summary[,5:8] / 1000000)\n";
    $r_string.="matrixfrac<-as.matrix(summary[,6:7] / summary[,5])\n";

    # Calculate ambiguous
    $r_string.="ambfile<-split(ambiguous,ambiguous\$V1)\n";

    # Setup the figures
    for my $type (keys %figures) {
	$r_string.=$figures{$type};
	$r_string.="par(oma=c(2,0,2,0))\n";
	$r_string.="layout(matrix(1:6,nrow=2,byrow=T),widths=rep(c(3,1,3),2))\n";
	$r_string.="cols=rainbow($lanes)\n";
	$r_string.="syms=rep(c(0:4),$lanes)\n";
	
	# Plot barplots
	$r_string.="barplot(t(matrix),beside=T,main=\"1:Read distribution\",col=c('orange','green','red','yellow'),ylab='Million reads',names.arg=c($names),cex.names=0.8,las=2)\n";
	$r_string.="axis(1,labels=F)\n";
	$r_string.="par(mfg=c(1,2))\n";
	$r_string.="old.mar<-par(\"mar\")\n";
	$r_string.="par(mar=c(0,0,0,0))\n";
	$r_string.="legend('center',legend=c('Total','Unambiguous','Ambiguous','Unique'),fill=c('orange','green','red','yellow'))\n";
	$r_string.="par(mar=old.mar)\n";
	$r_string.="barplot(t(matrixfrac),main=\"2: Fraction of reads with ambiguous bases\",col=c('green','red'),ylab='Fraction',names.arg=c($names),cex.names=$name_size,las=2)\n";
	$r_string.="axis(1,labels=F)\n";
	
	# Plot qualities
	$r_string.="plot(qualsfile[[1]]\$V2,qualsfile[[1]]\$V3,xlim=c(0,$length),xlab=\"Position\",ylab=\"Quality\",type=\"l\",col=cols[1],pch=0:4,main=\"3: Quality score by position\")\n";
	if (keys %{$stats} > 1) {
	    $r_string.="for (n in 2:length(qualsfile)) {lines(qualsfile[[n]]\$V2,qualsfile[[n]]\$V3,col=cols[[n]])}\n";
	}
	
	$r_string.="legend(\"topright\",c($names),col=cols,lty=1,pch=syms,cex=0.6,ncol=2)\n";
	
	# Plot ambiguous
	$r_string.="par(mfg=c(2,3))\n";
	$r_string.="plot(ambfile[[1]]\$V2,ambfile[[1]]\$V3,xlim=c(0,$length),ylim=c(0,$max),xlab=\"Position\",ylab=\"Ambiguous\",type=\"l\",col=cols[1],pch=0:4,main=\"Number of ambiguous bases (N) per position\")\n";
	if (keys %{$stats} > 1) {
	    $r_string.="for (n in 2:length(ambfile)) {lines(ambfile[[n]]\$V2,ambfile[[n]]\$V3,col=cols[[n]])}\n";
	}
	$r_string.="legend(\"topright\",c($names),col=cols,lty=1,pch=syms,cex=0.6,ncol=2)\n";
	$r_string.="title(main=\"$species $project Read summary\",outer=T)\n";
	$r_string.="dev.off()\n";
    }
	
    print $exec_fh $r_string;
    close($exec_fh);

    # execute the R file
    my $command="R --vanilla --quiet --slave < $execution_file";
    system($command);
    $command="rm $execution_file";
    system($command);
}

# This sub will process the fastq files a bit faster than the old one.
# There is however a problem, the qualities are assumed to be solexa illumina
# If this is not the case the script should be rerun with the correct settings.
### TO DO
# Guess the settings and after rerun using the correct quality type if we 
# guessed wrong
sub process_fastq_file_old {
    my $infn=shift;
    my $stats=shift;
    my $quals_pos=shift;
    my $tmpdir=shift;
    my $log_fh=shift;
    my $laneid=shift;
    my $ntpos=shift;
    my $qualities=shift;

    my $unique=0;
    my $read_length=0;
    my $ambiguous_reads=0;
    my $good_reads=0;
    my $total_reads=0;
    my $line_type=0;

    # determine the quality variant
    my $variant='illumina';
    if ($qualities eq 'phred') {
	$variant='sanger';
    }


    # Open the fasta file
    my $infh=Bio::SeqIO->new(-file => $tmpdir.'/'.$infn,
			     -format => 'fastq',
			     -variant => $variant);

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    my %ntindex=('A' => 1,
		 'C' => 2,
		 'G' => 3,
		 'T' => 4);

    # Process the fastq file
    while (my $seq_obj=$infh->next_seq()) {
	$total_reads++;

	my $seq=$seq_obj->seq();
	my $quals=$seq_obj->qual();

	# Count the number of unique reads
	print $tmpfh $seq,"\n";

	# Get the minimum read length
	my $seq_length=length($seq);
	if (!$read_length ||
	    ($read_length > $seq_length)) {
	    $read_length=$seq_length;
	}

	# Count the number of sequences with N's (ambiguous bases)
	if ($seq=~/N/) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}

	# Build the distribution of the N's and qualities per position
	my @nucleotides=split('',$seq);
	my @qualities=@{$quals};

	unless (scalar(@nucleotides) == scalar(@qualities)) {
	    die "Quality and nucleotide lengths differ\n";
	}

	# Using a single loop to go through both nt and qualities should save
	# processing time
	for (my $i=0;$i<@nucleotides;$i++) {
	    my $pos=$i + 1;

	    # Set the qualities
	    my $qual=$qualities[$i];
	    $quals_pos->{$laneid}->{$pos}+=$qual;

	    # Set the nucleotides
	    if ($nucleotides[$i]=~/[^ACTG]/o) {
		$ntpos->{$laneid}->{$pos}->[0]++;
	    } else {
		$ntpos->{$laneid}->{$pos}->[$ntindex{$nucleotides[$i]}]++;
		$ntpos->{$laneid}->{$pos}->[0]+=0;
	    }
	}
    }
    close($tmpfh);
    $infh->close();

    ### TO DO
    # Calculate the average qualities.
    # This is only useful if we know the type
    # of qualities We should guess the quality type here and after that print
    # the actual meaning of the qualities
    foreach my $pos (keys %{$quals_pos->{$laneid}}) {
	$quals_pos->{$laneid}->{$pos}=sprintf "%.2f",($quals_pos->{$laneid}->{$pos} / $total_reads);
    }

    # Get the unique reads without going out of the roof using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}

sub process_fasta_file {
    my $infn=shift;
    my $stats=shift;
    my $quals_pos=shift;
    my $tmpdir=shift;
    my $log_fh=shift;
    my $laneid=shift;
    my $ntpos=shift;

    my $unique=0;
    my $read_length=0;
    my $ambiguous_reads=0;
    my $good_reads=0;
    my $total_reads=0;

    # Open the fasta file
    my $infh=Bio::SeqIO->new(-file => $tmpdir.'/'.$infn,
			     -format => 'fasta');

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    my %ntindex=('A' => 1,
		 'C' => 2,
		 'G' => 3,
		 'T' => 4);

    # Process the fasta file
    while (my $seq_obj=$infh->next_seq()) {
	$total_reads++;

	my $seq=$seq_obj->seq;
	# Count the number of unique reads
	print $tmpfh $seq,"\n";

	# Get the minimum read length
	my $seq_length=length($seq);
	if (!$read_length ||
	    ($read_length > $seq_length)) {
	    $read_length=$seq_length;
	}

	# Count the number of Ns in the sequence (ambiguous bases)
	if ($seq=~/N/) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}

	# Build the distribution of the N's Actually get the non standard bases
	my @nucleotides=split('',$seq);
	for (my $i=0;$i<@nucleotides;$i++) {
	    my $pos=$i + 1;

	    # Set the qualities
	    $quals_pos->{$laneid}->{$pos}=50;

	    # Set the nucleotides
	    if ($nucleotides[$i]=~/[^ACTG]/o) {
		$ntpos->{$laneid}->{$pos}->[0]++;
	    } else {
		$ntpos->{$laneid}->{$pos}->[$ntindex{$nucleotides[$i]}]++;
		$ntpos->{$laneid}->{$pos}->[0]+=0;
	    }
	}
    }
    close($tmpfh);
    $infh->close();

   # Get the unique reads without going out of the rood using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}

sub process_fastq_file {
    my $infn=shift;
    my $stats=shift;
    my $quals_pos=shift;
    my $tmpdir=shift;
    my $log_fh=shift;
    my $laneid=shift;
    my $ntpos=shift;
    my $qualities=shift;

    my $unique=0;
    my $read_length=0;
    my $ambiguous_reads=0;
    my $good_reads=0;
    my $total_reads=0;
    my %sequence;
    my $line_type=0;

    # determine the quality variant
    my $variant='illumina';
    if ($qualities eq 'phred') {
	$variant='sanger';
    }

    # Open the fastq file. Here we will use our own sub, because some of the
    # files received do not comply whith the quality values they should actually
    # have and this makes them unreadable by the Bio::SeqIO::fastq module which
    # will always crash
    my $infh;
    $infh=get_fh($tmpdir.'/'.$infn);

    my $tmpfn=$$.'.tmp.sequences.txt';
    my $tmpfh=get_fh($tmpfn,1);

    my %ntindex=('A' => 1,
		 'C' => 2,
		 'G' => 3,
		 'T' => 4);

    # Process the fastq file
    while (my $line=<$infh>) {
	chomp($line);

	# Decide which line we are in
	if (($line_type == 0) &&
	    ($line=~s/^@//o)) {
	    # Initialize the sequence hash
	    %sequence=('id' => $line);
	    $line_type=1;
	    $total_reads++;
	    $line=<$infh>;
	    chomp($line);
	} elsif (($line_type==1) &&
		 ($line=~/^\+/o)) {
	    $line_type=2;
	    $line=<$infh>;
	    chomp($line);
	}

	if ($line_type == 1) {
	    $sequence{'seq'}.=$line;
	    next;
	} elsif ($line_type == 2) {
	    $sequence{'qual'}.=$line;
	    $line_type=0;
	    if (length($sequence{'qual'}) < length($sequence{'seq'})) {
		warn "Quality string shorter than sequence,skipping???\n";
		next;
	    } elsif (length($sequence{'qual'}) > length($sequence{'seq'})) {
		warn "Quality string longer than sequence,skipping???\n";
		next;
	    }
	} else {
	    warn "Problem with $infn. Shouldn't be here\n";
	}

	my $seq=$sequence{'seq'};
	# Count the number of unique reads
	print $tmpfh $seq,"\n";

	# Get the minimum read length
	my $seq_length=length($seq);
	if (!$read_length ||
	    ($read_length > $seq_length)) {
	    $read_length=$seq_length;
	}

	# Count the number of sequences with N's (ambiguous bases)
	if ($seq=~/N/o) {
	    $ambiguous_reads++;
	} elsif ($seq=~/\./o) {
	    $ambiguous_reads++;
	} else {
	    $good_reads++;
	}

	# Build the distribution of the N's and qualities per position
	my @nucleotides=split('',$seq);
	my @qualities=split('',$sequence{'qual'});
	unless (scalar(@nucleotides) == scalar(@qualities)) {
	    die "Quality and nucleotide lengths differe\n";
	}

	# Using a single loop to go through both nt and qualities should save
	# processing time
	for (my $i=0;$i<@nucleotides;$i++) {
	    my $pos=$i + 1;

	    # Set the qualities
	    my $qual=ord($qualities[$i]);
	    if ($variant eq 'sanger') {
		$qual-=33;
	    } else {
		$qual-=64;
	    }
	    $quals_pos->{$laneid}->{$pos}+=$qual;

	    # Set the nucleotides
	    if ($nucleotides[$i]=~/[^ACTG]/o) {
		$ntpos->{$laneid}->{$pos}->[0]++;
	    } else {
		$ntpos->{$laneid}->{$pos}->[$ntindex{$nucleotides[$i]}]++;
	    }
	}
    }
    close($tmpfh);
    $infh->close();

    ### TO DO
    # Calculate the average qualities.
    # This is only useful if we know the type
    # of qualities We should guess the quality type here and after that print
    # the actual meaning of the qualities
    foreach my $pos (keys %{$quals_pos->{$laneid}}) {
	$quals_pos->{$laneid}->{$pos}=sprintf "%.2f",($quals_pos->{$laneid}->{$pos} / $total_reads);
    }

    # Get the unique reads without going out of the roof using RAM
    my $command="sort -T $tmpdir $tmpfn | uniq| wc -l";
    $unique=`$command`;
    $command="rm $tmpfn";
    system($command);

    $unique=~s/[^\d]//g;
    return($good_reads,$ambiguous_reads,$total_reads,$read_length,$unique);
}
