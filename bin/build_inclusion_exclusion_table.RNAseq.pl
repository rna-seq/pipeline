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
# This script should take a gtf file and extract from it the exon information
# For each exon it will build a list of inclusion and exclusion events.
# Inclusion events are hits that are completely included in the exon, exclusion
# events are those that span the exon, including an exon from either side or
# the exon of interest.
# if a certain number of exons are selected in a separate file it will only
# calculate the exclusion junctions for those

use RNAseq_pipeline3 qw(get_fh parse_gff_line);
use RNAseq_pipeline_settings3 ('read_config_file');

# Declare some variables
my $annotation;
my $exclusionfile;

# Read the options file
my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$exclusionfile=$options{'EXCLUSIONFILE'};

if (-r $exclusionfile) {
    print STDERR "Exclusionfile: $exclusionfile is present. Skipping creation\n";
    exit;
}

# Read the annotation and extract a set og genes with its exons
my %genes=%{get_exons_from_file($annotation)};

# Print out the exons with the inclusion and exclusion events
get_inclusion_exclusion_events(\%genes,
			       $exclusionfile);

exit;

# This sub will take a list of exon ranges and for each one of them decide
# which are the inclusion and exclusion events
sub process_gene_exons {
    my $exons=shift;
    my %exclusion;

    for (my $i=0;$i < @{$exons}; $i++) {
	my @upstream;
	my @downstream;
	my $exon_range=$exons->[$i]->[1];
	my $exon_id=$exons->[$i]->[0];
	$exclusion{$exon_id}={};

	for (my $j=0;$j < @{$exons}; $j++) {
	    if ($i == $j) {
		next;
	    }
	    my $exon2_range=$exons->[$j]->[1];
	    my $exon2_id=$exons->[$j]->[0];

	    # Classify the exon
	    if ($exon_range->overlaps($exon2_range)) {
		next;
	    } elsif ($exon_range->start() > $exon2_range->end()) {
		push @upstream, $exon2_id;
	    } elsif ($exon_range->end() < $exon2_range->start()) {
		push @downstream, $exon2_id;
	    } else {
		warn "Unknown exon position $exon_id, $exon2_id\n";
	    }
	}

	foreach my $upstream (@upstream) {
	    my @up_coords=split('_',$upstream);
	    pop(@up_coords);
	    my $up_coord=pop(@up_coords);
	    pop(@up_coords);
	    my $chr=join('_',@up_coords);
	    foreach my $downstream (@downstream) {
		my @down_coords=split('_',$downstream);
		my $junction=join('_',
				  $chr,
				  $up_coord,
				  'splice',
				  $down_coords[-3]);
		$exclusion{$exon_id}{$junction}++;
	    }
	}
    }
    return(\%exclusion);
}

sub get_inclusion_exclusion_events {
    my $genes=shift;
    my $exclusionfile=shift;
    my %exclusion;

    my $gene_no=keys %genes;
    my $finished=0;
    
    print STDERR "Calculating and printing exclusion event table for $gene_no genes\n";
    print STDERR "This will take a while...\n";

    my $exfh=get_fh($exclusionfile,1);

    foreach my $gene (keys %{$genes}) {
	my @exons;
	foreach my $exon_id (keys %{$genes->{$gene}}) {
	    my @location=split('_',$exon_id);
	    my $strand=pop(@location);
	    my $end=pop(@location);
	    my $start=pop(@location);
	    my $range=Bio::Range->new(-start => $start,
				      -end => $end);
	    push @exons,[$exon_id,$range];
	}
	%exclusion=%{process_gene_exons(\@exons)};

	my @exon_list=keys %exclusion;
	foreach my $exon (@exon_list) {
	    my $junctions='-';
	    if (%{$exclusion{$exon}}) {
		$junctions=join(',',
				keys %{$exclusion{$exon}});
	    }
	    print $exfh join("\t",
			     $gene,
			     $exon,
			     $junctions),"\n";
	}

	# Try to free some memmory to prevent escalation;
	$genes{$gene}='';
	$finished++;
	unless ($finished % 100) {
	    print STDERR $finished,"\tFinished\r";
	}
    }
    close($exfh);

    print STDERR "All\tFinished\n";
}

sub get_exons_from_file {
    my $file=shift;

    print STDERR "Extracting exons from $file...";

    my %genes;
    my $fh=get_fh($file);

    while (my $line=<$fh>) {
	if ($line=~/^#/) {
	    next;
	}
	my %line=%{parse_gff_line($line)};

	if ($line{'type'} ne 'exon') {
	    next;
	}

	my $strand;
	if ($line{'strand'} eq '+') {
	    $strand=1;
	} elsif ($line{'strand'} eq '-') {
	    $strand=-1;
	} else {
	    warn "Unknown strand...\n";
	}

	my $exon_id=join('_',
			 $line{'chr'},
			 $line{'start'},
			 $line{'end'},
			 $strand);

	# Saving things as a hash avoids repeating the same exon many times if
	# it is present in different transcripts. This saves a lot of RAM when
	# we are checking for upstream and downstream exons
	$genes{$line{'feature'}{'gene_id'}}{$exon_id}=1;

    }

    close ($fh);

    print STDERR "done\n";

    return(\%genes);
}
