#!/soft/bin/perl
# DGK 2009

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
# This script should take a g[tf]f file and plot the number of mappings in known
# and in novel junctions in each of the files supplied

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_dbh get_coords_from_junc_id_sub get_junction_type_sub);

my $species;
my $project;
my $outfile;
my $table;
my $junctionsdir;

my %options=%{read_config_file()};
$table=$options{'JUNCTIONSTABLE'};
$species=$options{'SPECIES'};
$project=$options{'PROJECTID'};
$junctionsdir.=$options{'JUNCTIONSDIR'};

# Get subroutines
*junction_type=get_junction_type_sub($table);
*junction_coords=get_coords_from_junc_id_sub();

# Set the outfile
$outfile=$options{'PREFIX'}.'_unique_maps_junctions';


unless ($outfile) {
    die "No output file name\n";
}

my %files=%{read_file_list()};

my %lanes=%{get_lanes(\%files)};

foreach my $lane (keys %lanes) {
    my $type='single';
    my %distribution;
    my $infilename=$junctionsdir.'/'.$lane.'.'.$type.'.unique.gtf.gz';
    print STDERR "Processing $infilename\n";

    unless (-r $infilename) {
	die "$infilename is not readable\n";
    }

    build_distribution($infilename,
		       \%distribution);

    foreach my $junction (keys %distribution) {
	my $type=$distribution{$junction}[1];
	my @coords=@{junction_coords($junction)};
	if ($type eq 'novel') {
	    print join("\t",
		       $junction,
		       @coords,
		       $type,
		       $distribution{$junction}[0],
		       $lane),"\n";
	}
    }


    print STDERR "done\n";
}

exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[1]}++;
    }

    return(\%lanes);
}

sub build_distribution {
    my $file=shift;
    my $dist=shift;
    my $infh=get_fh($file);

    while (my $line1=<$infh>) {
	chomp($line1);
	my @line1=split("\t",$line1);

	my $line2=<$infh>;
	chomp($line2);

	my @line2=split("\t",$line2);
	if ($line1[8] ne $line2[8]) {
	    die "Wrong junctions:\n$line1\n$line2\n";
	}

	# get the splice coordinates
	my $chr=$line1[0];
	my $start=$line1[4];
	my $end=$line2[3];
	
	my $junc_id=join('_',$chr,$start,$end);
	my $type=junction_type($junc_id);

	$dist->{$junc_id}[0]++;
	$dist->{$junc_id}[1]="$type";
    }

    close($infh);
}
