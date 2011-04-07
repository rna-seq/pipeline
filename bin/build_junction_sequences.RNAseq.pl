#!/soft/bin/perl

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
# This script should use the database to collect the junction combinations
# and the exons fasta file in order to build a set of all junction sequences

use RNAseq_pipeline3 qw(get_fh get_exon_sequences);
use RNAseq_pipeline_settings3 qw(read_config_file get_dbh);
use Bio::SeqIO;

my $exonsfile;
my $junctionstable;
my $commondb;
my $junctionsfile;

my %options=%{read_config_file()};
$exonsfile=$options{'EXONSFASTA'};
$junctionstable=$options{'JUNCTIONSTABLE'};
$commondb=$options{'COMMONDB'};
$junctionsfile=$options{'JUNCTIONSFASTA'};

print STDERR "Extracting all junction sequences\n";

unless ($exonsfile &&
	$junctionstable) {
    die "ERROR:An exons file and a junctions table are required\n";
}

# Check if the file exists already
my $present=check_file($junctionsfile);

if ($present) {
    print STDERR "$junctionsfile is present. Skipping\n";
    exit;
}

# Connect to the databases
my $dbhcommon=get_dbh(1);

print STDERR "Processing entries from $junctionstable\n";

# Get the sequence of the exons
my %exons=%{get_exon_sequences($exonsfile)};

# Extract all the exon combinations for each of the junctions
my %junctions=%{get_junction_sequences($dbhcommon,
				       $junctionstable,
				       \%exons)};

# Close the database connection when finished
$dbhcommon->disconnect();

# Print out all junction sequences
print STDERR "Printing junction sequences in $junctionsfile\n";
my $juncfh=get_fh($junctionsfile,1);
foreach my $junc (keys %junctions) {
    my $seq=$junctions{$junc};
    print $juncfh ">$junc\n";
    for (my $pos=0;$pos < length($seq) ; $pos += 60) {
	print $juncfh substr($seq,$pos,60),"\n";
    }
}
close($junctionsfile);
print STDERR "done\n";

exit;

sub get_junction_sequences {
    my $dbh=shift;
    my $table=shift;
    my $exons=shift;

    my %junctions;

    my ($query,$sth);

    $query ='SELECT junction_id, exon1_id, exon2_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($junction_id,$exon1_id,$exon2_id) =$sth->fetchrow_array()) {
	my $seq1=$exons->{$exon1_id};
	my $seq2=$exons->{$exon2_id};

	# We need only one strand, as the other is the same and should have been
	# checked by a previous step
	my @name1=split('_',$exon1_id);
	my ($chr1,$start1,$end1,$strand1);
	$strand1=pop(@name1);
	$end1=pop(@name1);
	$start1=pop(@name1);
	$chr1=join('_',@name1);

	unless ($seq1 && $seq2) {
#	    print STDERR "exon1: $seq1\n";
#	    print STDERR "exon2: $seq2\n";
	    warn "Exon sequence missing $exon1_id $exon2_id. Skipping...\n";
	    next;
	}

	# Depending on the strand we need to put one or the other exon first
	# as in the junction ID always the 5' exon is first regrdless of the
	# strand
	# We no longer need information on mismatches, etc... as the junction
	# mappngs will be checked when mapped
	if ($strand1 == 1) {
	    $junctions{$junction_id}=$seq1.$seq2;
	} elsif ($strand1 == -1) {
	    $junctions{$junction_id}=$seq2.$seq1;
	} else {
	    die "Wrong strand $strand1 in $exon1_id\n";
	}
    }
    return(\%junctions);
}

sub check_file {
    my $file=shift;

    my $present=0;

    my $tablefh=get_fh($options{'PREFIX'}.'_junction_seqs.txt',1);
    if (-r $file) {
	$present=1;
	# Print the table location for the junctions of this experiment
	print $tablefh join("\t",
			    $file,
			    "Present"),"\n";
    } else {
	# Continue, as the table must be created
	print STDERR $file,"\tIs not present\n";
	print $tablefh join("\t",
			    $file,
			    "Generated"),"\n";
    }
    close ($tablefh);

    return($present);
}
