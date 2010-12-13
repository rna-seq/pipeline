#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take the necessary information from the table containing
# the junctions identified and the split_mapping_breakdown and it will build
# a table with the numbers and percentages (when appliable)

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file','read_file_list',
			       'get_junction_type_sub');

# Declare some variables
my $prefix;
my $junctiondir;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$junctiondir=$options{'JUNCTIONSDIR'};

# Get the table and file names we need
my $splittable=$prefix.'_split_mapping_breakdown';
my $junctable=$prefix.'_junctions';

# Get the dbh
my $dbh=get_dbh();

# Get the junctions annotation table
my $juncannot=get_junct_annot_table($junctable,
				    $dbh);

# Get some useful subs
*junct_type=get_junction_type_sub($juncannot);

# get the file names
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# Collect the information from the junctions
my %totals=%{get_total_junctions($juncannot)};
my %junctions;
foreach my $lane (keys %lanes) {
    my $type;
    if (keys %{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (keys %{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    # Get all the junctions files
    my @junctionfns;
    foreach my $track (keys %{$lanes{$lane}}) {
	my $juncfilename=$junctiondir.'/'.$track.'.single.unique.overlap.total';
	push @junctionfns, $juncfilename;
    }

    # Process each of the junction files
    foreach my $juncfilename (@junctionfns) {
	if (-r $juncfilename) {
	    print STDERR "Processing $juncfilename\n";
	} else {
	    die "Can't read $juncfilename\n";
	}
	my $juncfh=get_fh($juncfilename);

	while (my $line=<$juncfh>) {
	    my ($juncid,$reads)=split("\t",$line);

	    # Allow compatibility with older files
	    if ($juncid=~/splice_/) {
		$juncid=~s/splice_//;
	    }
	    my $type=junct_type($juncid);
#	    my $type=('known','novel')[rand(2)];
	    $junctions{$type}{$juncid}=1;
	}
    }
}

# collect the information from the split_mapping_breakdown table
my $splitjunc=get_split_junctions($splittable,
				  $dbh);

# Print out the summary
my $detectedknown=keys %{$junctions{'known'}};
my $detectednovel=keys %{$junctions{'novel'}};
print join("\t",
	   'Known',
	   $detectedknown,
	   $totals{'known'}),"\n";
print join("\t",
	   'Novel',
	   $detectednovel,
	   $totals{'novel'}),"\n";
print join("\t",
	   'Unannotated',
	   $splitjunc,
	   '\N'),"\n";

exit;

sub get_split_junctions {
    my $table=shift;
    my $dbh=shift;

    my ($query,$sth);

    $query ='SELECT count(DISTINCT location) ';
    $query.="FROM $table ";
    $query.='WHERE type != ?';
    $sth=$dbh->prepare($query);
    $sth->execute('close');

    my ($locations)=$sth->fetchrow_array();

    return($locations);
}

# This should get a list of types from the junctions table as well as the
# total number of junctions of that type in the database
# There is a slight redundancy due to alternative splicing. This should be
# checked
sub get_total_junctions {
    my $table=shift;
    my $dbh=get_dbh(1);

    my %totals;

    my ($query,$sth);

    $query ='SELECT type, count(*) ';
    $query.="FROM $table ";
    $query.='GROUP BY type';

    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($type,$number)=$sth->fetchrow_array()) {
	$totals{$type}=$number;
    }

    $dbh->disconnect();

    return(\%totals);
}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}

sub get_junct_annot_table {
    my $table=shift;
    my $dbh=shift;

    my ($query,$sth);

    $query ='SELECT table_id ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    my $count=$sth->execute();

    unless ($count == 1) {
	die "Problem getting the table info from $table\n";
    }

    my ($annottable)=$sth->fetchrow_array();

    return($annottable);
}
