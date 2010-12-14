#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/.*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script will take the necessary commands for running the recursive mapping
# on the cluster

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh send2cluster);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list get_dbh);
use RNAseq_GEM3 ('check_index','determine_quality_type','get_mapper_routines',
		 'check_input');

my $index;
my $mapper;
my $outdir;
my $filetype;
my $mismatches;
my $paralleldir;
my $readdata;
my $projdir;
my $bindir;
my $splittable;
my $splitdir;
my $recmapdir;
my $queue='mem_6';

# Get some options from the command line
GetOptions('index|I=s' => \$index,
	   'outdir|o=s' => \$outdir);

my %options=%{read_config_file()};
$mapper=$options{'MAPPER'};
$mismatches=$options{'MISMATCHES'};
$paralleldir=$options{'LOCALPARALLEL'};
$readdata=$options{'READDIR'};
$projdir=$options{'PROJECT'};
$bindir=$options{'BIN'};
$splittable=$options{'PREFIX'}.'_split_mapping';
$splitdir=$projdir.'/'.$options{'SPLITMAPDIR'};
$recmapdir=$projdir.'/'.$options{'RECMAPDIR'};
$queue=$options{'CLUSTER'};

# Decide where to put the output
unless ($outdir) {
    die "Where should I put the results ???\n";
}
$outdir=~s/\/$//;
my $usecluster=0;

if ($projdir=~/^\/users/ &&
    $index=~/^\/users/) {
    print STDERR "Mapping in the cluster\n";
    $usecluster=1;
} else {
    print STDERR "I can only use the cluster from /users. Mapping locally\n";
}

# Make sure we have a valid index
unless ($index) {
    die "An index file is required\n";
}

# check for the existence of the index
my $index_ok=check_index($index);
unless ($index_ok) {
    die $index,"\tIs not a valid index\n";
}

# Get the list of files we need to map
my %files=%{split_file_list($splittable,
			    $splitdir)};

# Check each of the files and determine if it is present in the paralleldir
# directory if not unzip it and copy it there.
my %locations=%{check_read_files(\%files,
				 $paralleldir,
				 $splitdir)};

my @filepairs;
foreach my $file (keys %locations) {
    my $qualities;
    # Check the input file
    my $infile=$locations{$file};

    # Make the outfile
    my $outdir=$recmapdir;

    push @filepairs,[$infile,$outdir];

}

unless (@filepairs) {
    print STDOUT "Everything seems to be mapped already\n";
    exit;
}

if ($usecluster) {
    # Build the submission file
    my $subfile=build_run_mapper_submission(\@filepairs,
					    $bindir,
					    $index);
    send2cluster($subfile,
		 $queue);

    # clean up
    my $command="rm $subfile";
    print STDERR "Executing: $command\n";
    system($command);

    } else {
	my %mapper=%{get_mapper_routines()};
	foreach my $pair (@filepairs) {
	    # Map the file
	    # Run the recursive mapper for each of the files
	    my $command="run_recursive_mapper.RNAseq.pl ";
	    my $genomeindex=$options{'GENOMEINDEX'};
	    my $recmapdir=$options{'RECMAPDIR'};
	    $command.="-index $genomeindex ";
	    $command.='-i '.$pair->[0].' ';
	    $command.="-o $recmapdir > $recmapdir/";
	    $command.=$pair->[0].'.rec.mapping.log';
	    print STDERR "Executing: $command\n";
	}
}
exit;

# This subroutine will get a list of the split-mapped files that should be
# present
sub split_file_list {
    my $table=shift;
    my $splitdir=shift;

    my %files;
    my $dbh=get_dbh();

    my ($query,$sth);
    $query ='SELECT filename ';
    $query.="FROM $table";
    $sth=$dbh->prepare($query);
    $sth->execute();

    while (my ($file)=$sth->fetchrow_array()) {
	$files{$file}=1;
    }

    return(\%files);
}

sub build_run_mapper_submission {
    my $pairs=shift;
    my $bidir=shift;
    my $index=shift;

    print STDERR 'Building submission file...';
    my $filenum=@{$pairs};
     
    unless(@{$pairs}) {
	die "No input supplied\n";
    }
    
    # Get the input and output files
    my @infiles;
    my @outfiles;
    foreach my $pair (@{$pairs}) {
	push @infiles,$pair->[0];
	push @outfiles,$pair->[1];
    }

    # Print the submission file
    my $subfile="subfile.$$.job";
    my $outfh=get_fh($subfile,1);
    
    print $outfh <<FORMEND;
# Get the job name
#\$ -N RNAseqRecMap
    
# Set the array jobs
#\$ -t 1-$filenum

# Request 8 cpus this cannot be done, but we can request memmory
#\$ -l h_vmem=16G

# Write in to the current working directory
#\$ -cwd 
export PATH=\$PATH:/soft/bin
infiles=(@infiles)
outfiles=(@outfiles)

export infile=\${infiles[\$SGE_TASK_ID-1]}
export outdir=\${outfiles[\$SGE_TASK_ID-1]}

echo \$HOSTNAME >&2
$bindir/run_recursive_mapper.RNAseq.parallel.pl -index $index -infile \$infile -outdir \$outdir > \$infile.mapping.log
FORMEND
    ;
    close($outfh);
    print STDERR "done\n";

    return($subfile);
}

sub check_read_files {
    my $files=shift;
    my $paralleldir=shift;
    my $readdir=shift;

    my %locations;

    foreach my $file (keys %{$files}) {
	print STDERR "Checking $file...";
	my $filepath1=$paralleldir.'/'.$file;
	my $filepath3=$readdir.'/'.$file;
	if ( -r $filepath1) {
	    print STDERR "Present at $readdir\n";
	    $locations{$file}=$filepath1;
	}  elsif (-r $filepath3) {
	    print STDERR "Present at $readdir\n";
	    $locations{$file}=$filepath3;
	} elsif (-r "$filepath3.gz") {
	    print STDERR "Unzipping in $paralleldir\n";
	    my $command="gunzip -c $filepath3.gz > $filepath1";
	    system($command);
	    $locations{$file}=$filepath1;
	} else {
	    die "I can't find file $file\n";
	}
    }

    return(\%locations);
}
