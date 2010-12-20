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
# This script should take a list of mapping files and extract the reads that
# do not map to any place in any of them

use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

my $tmpdir;
my $file_list;
my $genomedir;
my $transdir;
my $juncdir;
my $mapper;
my $splitmapdir;

my %options=%{read_config_file()};
$tmpdir=$options{'LOCALDIR'};
$file_list=$options{'FILELIST'};
$genomedir=$options{'GENOMEDIR'};
$transdir=$options{'TRANSDIR'};
$juncdir=$options{'JUNCTIONSDIR'};
$mapper=$options{'MAPPER'};
$splitmapdir=$options{'SPLITMAPDIR'};

# Get the subroutine for parsing the reads
*maps_parser=get_mapped_fh($mapper);


my %filelist=%{read_file_list($file_list)};
my %lanes=%{get_lanes(\%filelist)};

foreach my $infile (keys %filelist) {
    print STDERR "Processing: $infile\n";
    my $splitroot=$lanes{$infile};
    my $fileroot=$infile;
    $fileroot=~s/.fa(stq)*$//;

    my $unmappedfn=$splitmapdir.'/'.$splitroot.'.unmapped.fa';

    # skip if the file is present
    if (-r $unmappedfn ||
	-r $unmappedfn.'.gz') {
	print STDERR $unmappedfn,"\tPresent. Skipping...\n";
	next;
    } else {
	print STDERR "Building: $unmappedfn\n";
    }

    # Get the necessary filehandles
    my $genomefh=get_fh(`ls $genomedir/$fileroot.gem.map*`);
    my $transfh=get_fh(`ls $transdir/$fileroot.gem.map*`);
    my $juncfh=get_fh(`ls $juncdir/$fileroot.gem.map* | grep -v 'gen.coords'`);
    my $unmappedfh=get_fh($unmappedfn,1);

    print STDERR "Printing unmapped reads from $fileroot to $unmappedfn\n";
    
    while (my $genline=<$genomefh>) {
	my $transline=<$transfh>;
	my $juncline=<$juncfh>;

	# Check we have lines from each file
	unless ($genline && $transline && $juncline) {
	    die "Problem: one read entry seems to be missing\n";
	}
	my $hitsgen=maps_parser($genline);
	my $hitstrans=maps_parser($transline);
	my $hitsjunc=maps_parser($juncline);

	# Check the read we are comparing is the same in the three files
	unless (($hitsgen->[0] eq $hitstrans->[0]) &&
		($hitsgen->[0] eq $hitsjunc->[0])) {
	    print STDERR join("\t",
			      'genome',$hitsgen->[0]),"\n";
	    print STDERR join("\t",
			      'trans',$hitstrans->[0]),"\n";
	    print STDERR join("\t",
			      'junc',$hitsjunc->[0]),"\n";
	    die "Problem: Incorrect read id\n";
	}

	# Determine if it had a hit or not
	if ($hitstrans->[1] || $hitsgen->[1] || $hitsjunc->[1]) {
	    next;
	} else {
	    # This read is unmapped
	    print $unmappedfh ">",$hitsgen->[2],"\n";
	    print $unmappedfh $hitsgen->[0],"\n";
	}

    }
    close($genomefh);
    close($transfh);
    close($juncfh);
    close($unmappedfh);
}


exit;

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$file}=$files->{$file}->[1];
    }

    return(\%lanes);
}

sub get_mapped_fh {
    my $mapper=shift;

    my %parsers;

    $parsers{'GEM'}=sub {
	my $line=shift;

	chomp($line);
	my @line=split(/\t/,$line);
	my $read_id=$line[0];
	my $read_seq=$line[1];
	
	# determine if we are dealing with qualities
	unless ($line[2]=~/^\d[\d:]+\d$/) {
	    splice(@line,2,1);
	}
	
	# if no matches return 0 if matches return1
	if ($line[2]=~/^0(:0)*$/) {
	    return([$read_seq,0,$read_id]);
	} else {
	    return([$read_seq,1]);
	}
    };

    if ($parsers{$mapper}) {
	return($parsers{$mapper});
    } else {
	die "Unknown mapper: $mapper\n";
    }
}
