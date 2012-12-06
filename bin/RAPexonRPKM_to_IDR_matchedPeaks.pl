#!/usr/bin/perl
#RAPexonRPKM_to_IDR_matchedPeaks.pl

###################################################################################################################################################
### This script parses next RAP output for exon and gene RPKMs two bioreplicates into 'match peak format" for IDR analysis (version 31Aug 2010) ###

# RAP output: exon_id	RPKM	datasets(lanes)_id
# chr20_3209741_3210098_-1	0.052	All
# chr7_154774929_154775454_-1	5.625	All
# chr7_135047751_135047938_-1	6.301	All

#matched peakformat (used for IDR code from 06-24-10 by Quanha Li)
# # chr20_3209741_3210098_-1	0.052	chr20_3209741_3210098_-1	0.052	#peak present in both reps
# # chr7_154774929_154775454_-1	5.625	-1	-1				#peak absent in second file
# # -1	-1	chr7_135047751_135047938_-1	6.301				#peak absent in first file

#usage: ~/bin/RAPexonRPKM_to_IDR_matchedPeaks.pl --input1 file1 --input2 file2 --output somemeaningfulname

#e.g.:  ~/bin/RAPexonRPKM_to_IDR_matchedPeaks.pl --input1 /users/rg/projects/NGS/Projects/ENCODE/hg19main/GingerasCarrie/001WC/exons/All.exon.rpkm.pooled.txt.gz --input2 /users/rg/projects/NGS/Projects/ENCODE/hg19main/GingerasCarrie/002WC/exons/All.exon.rpkm.pooled.txt.gz --output test
##################################################################################################################################################

use strict;
use warnings;

use Getopt::Long;

########################################################
## Input/ output &checks
########################################################

#define some variables
my $input1;    #bioreplicate 1
my $input2;    #bioreplicate 2
my $output;

GetOptions(
	'input1|i1=s' => \$input1,
	'input2|i2=s' => \$input2,
	'output|o=s' => \$output);

open( my $infile1, "zcat $input1 | " ) or die "couldn't open input file";
open( my $infile2, "zcat $input2 | " ) or die "couldn't open input file";

die "I have no input files\n" unless $infile1 or $infile2;

open( my $outfile, "> $output.matchedPeaks.txt");

##########################################
## create 'matched peak IDR file'
###########################################


my %hash;                           			#create hash

print "creating hash...\n";

while (my $line =  <$infile1>) {                	#read in the first file
	chomp($line);
	my @tmp1 = split( "\t", $line );
	my $ID = $tmp1[0];   			      	#extract coordinates and
	my $value1 = $tmp1[1];                        	#signal value

	$hash{$ID} = [ $value1, 0 ];                 	#fill the hash

}


## print number of elements predicted in Replicate1
print scalar( keys %hash ) . " elements present in " . $input1 . "\n";

my $count;

while (my $line = <$infile2>) {                          #read in the second file
	chomp($line);
	my @tmp2 = split( "\t", $line );
	my $ID = $tmp2[0];    				 #extract coordinates and
	my $value2 = $tmp2[1];                       	 #signal value

	if ( $hash{$ID} ) {                          	 #keep filling the hash

		# print $hash{$ID} . "exists adding 2nd value\n";

		$hash{$ID}->[1] = $value2; 					#1. if the elementID already exists as key, add the second signal value (bioreplicate2) only

	} else { 										#2. if the elementID does not yet exists as key, create it and add no signal value for bioreplicate1 and signal value

	  # print "key doesnt exist, adding key and value\n";					#for bioreplicate2

		$hash{$ID} = [ 0, $value2 ];	
	}

  $count++;

}


print $count . " elements present in" . $input2 . "\n";

foreach my $ID ( keys %hash ) {    					# fill the dummy array with the elementIDs and signal value...
		
		my @values = ( -1, -1, -1, -1 );    		# make a dummy array
		
		if ( $hash{$ID}->[0] ) {    				# bioreplicate1 -> column1 = elementID; column2 = signal value
	
			$values[0] = $ID;
			$values[1] = $hash{$ID}->[0];			# use reference to look up the signal value

		}

		if ( $hash{$ID}->[1] ) {				    # bioreplicate2 -> column3 = elementID; column4 = signal value

			$values[2] = $ID;
			$values[3] = $hash{$ID}->[1];			#use reference to look up the signal value

		}
		
	print $outfile join("\t", @values) . "\n";
}


close $outfile;
