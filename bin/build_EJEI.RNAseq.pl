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

use Pod::Usage;
use RNAseq_pipeline3 qw(get_fh parse_gff_line);
use RNAseq_pipeline_settings3 ('read_config_file','get_dbh','read_file_list',
			      'get_gene_from_short_junc_sub'
    );

# Declare some variables
my $help;
my $man;
my $prefix;
my $junctionsdir;
my $tmpdir;

# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Read the configuration file
my %options=%{read_config_file()};
$prefix=$options{'PREFIX'};
$junctionsdir=$options{'JUNCTIONSDIR'};
$tmpdir=$options{'LOCALDIR'};

# Connect to the database
my $dbh=get_dbh();
my $commondbh=get_dbh(1);

# get some useful subroutines
*junc2gene=get_gene_from_short_junc_sub($commondbh);

# Get the lane names;
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# Read and process the junction overlap.total files
foreach my $lane (keys %lanes) {
    my $type;
    if (keys %{$lanes{$lane}} == 1) {
	$type='single';
    } elsif (keys %{$lanes{$lane}} == 2) {
	$type='paired';
    } else {
	die "Unknown type\n";
    }

    # Get the junction files
    my @junctionfns;
    foreach my $track (keys %{$lanes{$lane}}) {
	my $juncfilename=$junctionsdir.'/'.$track.'.single.unique.overlap.total';
	push @junctionfns, $juncfilename;
    }

    # Read in the junction coverage
    my %junction_reads;
    my %gene_reads;
    foreach my $juncfilename (@junctionfns) {
	if (-r $juncfilename) {
	    print STDERR "Processing $juncfilename\n";
	} else {
	    die "Can't read $juncfilename\n";
	}
	get_feature_coverage_junctions($juncfilename,
				       \%junction_reads,
				       \%gene_reads);
    }

    # Print out the information
    my $outfile=$junctionsdir.'/'.$lane.'.'.$type.'.EJEI.txt.gz';
    process_features(\%junction_reads,
		     \%gene_reads,
		     $outfile,
		     $lane);
}



exit;

sub process_features {
    my $junclist=shift;
    my $gene_reads=shift;
    my $output=shift;
    my $lane=shift;

    print STDERR "Collecting events...\n";

    my $outfh=get_fh($output,1);

    foreach my $junc (keys %{$junclist}) {
	my $gene=junc2gene($junc);
	my $total=$gene_reads->{$gene};
	my $ejei=sprintf "%.3f",($junclist->{$junc} / $total);

	# Print results only if they are positive
	print $outfh join("\t",
			  $gene,
			  $junc,
			  $ejei,
			  $total,
			  $lane),"\n";
    }
    close($outfh);

    print STDERR "\rDone\n";
}

sub get_feature_coverage_junctions {
    my $infn=shift;
    my $features=shift;
    my $gene_reads=shift;

    print STDERR "Getting mappings from $infn...";

    my $infh=get_fh($infn);

    while (my $line=<$infh>) {
	chomp($line);
	my ($feature,$coverage)=split("\t",$line);

	my $gene_id=junc2gene($feature);
	$gene_reads->{$gene_id}+=$coverage;

	$features->{$feature}+=$coverage;
    }
    print STDERR "done\n";
}

sub get_lanes {
    my $files=shift;
    my %lanes;
    
    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[0]}{$files->{$file}->[1]}=1;
    }

    return(\%lanes);
}

__END__

=head1 NAME
    
    get_EJEI.pl
    
=head1 SYNOPSIS
    
get_EJEI.pl [options] 
    
  Options:
    -help            brief help message
    -man             full documentation
        
=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.
    
=back
    
=head1 DESCRIPTION
    
    get_EJEI.pl should calculate for each exon in the study the EJEI (which is
    the Exon Junction Expression Index as defined in Dong et al. This is
    basically the fraction of reads that fall in a certain junction compared
    to the total number of reads mapping to junctions in that gene.

=cut
