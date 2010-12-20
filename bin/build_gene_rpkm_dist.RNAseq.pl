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
# This scritp will take information from the tables containing the RNAseq data
# and the RNAseq pooled data and it will build the distribution of teh RPKMs

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file');

# Declare some variables
my $prefix;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};

# Get the table names we need
my $rpkm=$prefix.'_gene_RPKM';
my $rpkm_pooled=$prefix.'_gene_RPKM_pooled';

# Get the information for the histogram
my %hist;
my $dbh=get_dbh();
get_info_from_table($dbh,
		    $rpkm,
		    %hist);
get_info_from_table($dbh,
		    $rpkm_pooled,
		    %hist);

# Print out the result
foreach my $value (keys %hist) {
    foreach my $set (keys %{$hist{$value}}) {
	print join("\t",
		   $hist{$value}{$set},
		   $value,
		   $set),"\n";
    }
}

exit;

sub get_info_from_table {
    my $dbh=shift;
    my $table=shift;
    my $hist=shift;

    my ($query,$sth,$count);
    $query ='SELECT * ';
    $query.="FROM $table ";
    $sth=$dbh->prepare($query);
    $count=$sth->execute();

    unless($count && 
	   ($count > 0)) {
	die "No entries found in $table\n";
    }

    while (my ($gene,$rpkm,$set)=$sth->fetchrow_array()) {
	my $value=int($rpkm + 0.5);
	if ($value) {
	    $hist{$value}{$set}++;
	}
    }
}
