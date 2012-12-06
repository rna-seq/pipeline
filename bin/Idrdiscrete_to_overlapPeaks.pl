#!/usr/bin/perl
# idrdiscrete_to_overlapPeaks.pl

##############################################################################################
# Objective: This is to parse the idr output (discrete code version 2) "*-categories.txt" into the "*-overlapped-peaks.txt" fromat from the former version of the idr code
################################################################################################
#usage: 

#infile: e.g. Ging.001WCMinus.vs.002WCMinus.sj.-categories.txt
# "score1" "score2" "counts" "idr" "IDR"
# -1 1 44925 0.999999305744421 0.564235553965311
# -1 2 14186 0.999999999773239 0.620293514313659
# -1 3 5453 0.99999999999999 0.638605454439926

# ref: e.g. Ging.001WC.vs.002WC.sj.matchedPeaks.txt
# chr22_39714597_39715732	3	chr22_39714597_39715732	7
# chr3_52012098_52021325	1	-1	-1
# -1	-1	chr14_73420888_73420706	1
# chr1_36383945_36383619	1	-1	-1
# -1	-1	chr19_9473760_9486992	1

#outfile
# #head Ging.001WC.vs.002WC.sj.-overlapped-peaks.txt
# "ID1" "score1" "ID2" "score2" "idr"
# "1" "chr22_39714597_39715732" 3 "chr22_39714597_39715732" 7 0.879666118541438
# "2" "chr1_24304142_24305203" 2 "chr1_24304142_24305203" 3 0.855663346513335
# "3" "chr13_25825913_25825991" 61 "chr13_25825913_25825991" 27 0.000604801404084587
# "4" "chr19_3750182_3750580" 2 "chr19_3750182_3750580" 4 0.84651367958781

################################################################################################

use strict ; 
use Data::Dumper ; 
use Getopt::Long;

#define some variables
my $filename1;    				# *-categories.txt
my $filename2;      				# *matchedPeaks.txt
my $filename3;					# output file prefix


#get options

GetOptions(
	'infile|i=s' => \$filename1,
	'ref|r=s' => \$filename2,
	'outfile|o=s' => \$filename3);


# create hash
my %hash ; 


open(my $file_1 , "< $filename1") or die $! ; 		#read *-categories.txt

<$file_1>;						#skip the first line
while( my $line = <$file_1>){

  my @parts = split(/\s+/ , $line) ; 			#assign scores
  my $score1 = $parts[0] ; 
  my $score2 = $parts[1] ; 

  $hash{ $score1 . '_' . $score2 } = $parts[3] ;  	#fill hash (use unique identifier)

}

#print Dumper(\%hash) ; 


open (my $file2 , "< $filename2") or die $! ; 				#read *matchedPeaks.txt
open( my $outfile, "> $filename3.dc2.-overlapped-peaks.txt" ) or die $! ; 	#open output file

print $outfile " \t\"ID1\"\t\"score1\"\t\"ID2\"\t\"score2\"\t\"idr\"\n" ;  


while(my $line_2 = <$file2>){				

    my @parts2 = split(/\s+/ , $line_2) ; 			

    my $junction1 = $parts2[0] ; 					#assign junction ids and scores
    
    my $f2_score1 = $parts2[1] ; 
    my $junction2 = $parts2[2] ; 
    my $f2_score2 = $parts2[3] ; 

    my $junction_new ; 							# alternative assignment if id only present in one bioreplicate (1,2 col or 3,4 col)
    if ($f2_score1 == -1 ){
	$junction_new = $junction2 ;
    }else{
	$junction_new = $junction1
    }

   
    my $idr = $hash{$f2_score1 . '_' . $f2_score2} ; 				# pull out the (idr) value from the hash; (equivalent to something like 											$part[0], which pulls out the first element from the array @part)
    
  
  print $outfile $. . "\t$junction_new\t$f2_score1\t$junction_new\t$f2_score2\t$idr\n" ;  

}

