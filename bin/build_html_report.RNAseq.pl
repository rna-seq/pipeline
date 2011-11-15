#!/soft/bin/perl

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

#    Author : David Gonzalez, david.gonzalez@crg.eu

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
# This script should work as a cgi script

use CGI qw(th td Tr caption table);
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_pipeline_settings3 qw(read_config_file read_file_list);

my $graphsdir;
my $project_id;
my $prefix;
my $reportfn='report.html';

# Read the configuration file
my %options=%{read_config_file()};
$graphsdir=$options{'GRAPHS'};
$project_id=$options{'PROJECTID'};
$prefix=$options{'PREFIX'};

# Build the start page
build_report_main($project_id,
		  $reportfn);

# Build the Read stats page
build_report_reads($project_id,
		   $prefix,
		   $graphsdir);

# Build the mapping page
build_report_mapping($project_id,
		     $prefix,
		     $graphsdir);

# Build the unique mappings page
build_unique_mapping($project_id,
		     $prefix,
		     $graphsdir);

# Build the cluster page
build_report_clusters($project_id,
		      $prefix,
		      $graphsdir);

# Build a page for the expression information
build_report_expression($project_id,
		       $prefix,
		       $graphsdir);

exit;

sub build_report_expression {
    my $project_id=shift;
    my $prefix=shift;
    my $graphsdir=shift;
    my $reportreads='expression_info.html';

    my $reportfh=get_fh($reportreads,1);

    # Get the gene data
    my $genesdetps=$graphsdir.'/'.$prefix.'.Genes.detection.ps';
    my $genesdetjpeg=$graphsdir.'/'.$prefix.'.Genes.detection.jpeg';
    my $genesexpps=$graphsdir.'/'.$prefix.'.Genes.split.ps';
    my $genesexpjpeg=$graphsdir.'/'.$prefix.'.Genes.split.jpeg';
    my $satps=$graphsdir.'/Gene.saturation.ps';
    my $satjpeg=$graphsdir.'/Gene.saturation.jpeg';

    # Get teh exon data
    my $exdetps=$graphsdir.'/'.$prefix.'.Exons.detection.ps';
    my $exdetjpeg=$graphsdir.'/'.$prefix.'.Exons.detection.jpeg';
    my $exexpps=$graphsdir.'/'.$prefix.'.Exons.split.ps';
    my $exexpjpeg=$graphsdir.'/'.$prefix.'.Exons.split.jpeg';
    my $exsatps=$graphsdir.'/Exon.saturation.ps';
    my $exsatjpeg=$graphsdir.'/Exon.saturation.jpeg';

    # Print the html
    my $q=CGI->new();
    print $reportfh $q->start_html('RNAseq Report');
    print $reportfh $q->h1("Mapping stats for $project_id");
    # Print the gene expression stats
    print $reportfh table({-border=>undef},
			  caption('Gene Expression stats'),
			  Tr({-align=>'CENTER',-valign=>'CENTER'},
			     [
			      th(['Description', 'File']),
			      td(['Saturation curve for gene detection' , "<a href=\"$satps\"> <img src=\"$satjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Gene detection' , "<a href=\"$genesdetps\"> <img src=\"$genesdetjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Gene Expression by rank' , "<a href=\"$genesexpps\"> <img src=\"$genesexpjpeg\" heigth=100 width=100 h/> </a>"]),
			     ]
			  )
	);

    # Print the Exon expression stats
    print $reportfh table({-border=>undef},
			  caption('Exon Expression stats'),
			  Tr({-align=>'CENTER',-valign=>'CENTER'},
			     [
			      th(['Description', 'File']),
			      td(['Saturation curve for exon detection' , "<a href=\"$exsatps\"> <img src=\"$exsatjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Exon detection' , "<a href=\"$exdetps\"> <img src=\"$exdetjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Exon Expression by rank' , "<a href=\"$exexpps\"> <img src=\"$exexpjpeg\" heigth=100 width=100 h/> </a>"]),
			     ]
			  )
	);

    print $reportfh $q->end_html();
    close($reportfh);
}

sub build_report_clusters {
    my $project_id=shift;
    my $prefix=shift;
    my $graphsdir=shift;
    my $reportreads='read_clusters.html';

    my $reportfh=get_fh($reportreads,1);
    
    my $clusterstatsps=$graphsdir.'/'.$prefix.'_initial_clusters.ps';
    my $clusterstatsjpeg=$graphsdir.'/'.$prefix.'_initial_clusters.jpeg';
    # Print the html
    my $q=CGI->new();
    print $reportfh $q->start_html('RNAseq Report');
    print $reportfh $q->h1("Read stats for $project_id");
    print $reportfh table({-border=>undef},
			  Tr({-align=>'CENTER',-valign=>'CENTER'},
			     [
			      th(['Description', 'File']),
			      td(['Read Stats' , "<a href=\"$clusterstatsps\"> <img src=\"$clusterstatsjpeg\" heigth=100 width=100 h/> </a>"])
			     ]
			  )
	);
    
    print $reportfh $q->end_html();
    close($reportfh);

}

sub build_unique_mapping {
    my $project_id=shift;
    my $prefix=shift;
    my $graphsdir=shift;
    my $reportreads='unique_mappings.html';

    my $reportfh=get_fh($reportreads,1);

    # Get the transcriptome tables
    my %files=%{read_file_list()};

    my %lanes=%{get_transcript_lanes(\%files)};

    # All transcripts
    my @transcriptfilesall;
    foreach my $lane (keys %lanes) {
	my $type='single';
	my $infilename=$graphsdir.'/'.$lane.'.'.$type.'.unique.gtf.gz.all.stats';
	my $string="<a href=\"$infilename.ps\"> <img src=\"$infilename.jpeg\" heigth=100 width=100 h/> </a>";
	push @transcriptfilesall,$string;
    }

    # Breakdown
    my @transcriptfilesbreak;
    foreach my $lane (keys %lanes) {
	my $type='single';
	my $infilename=$graphsdir.'/'.$lane.'.'.$type.'.unique.gtf.gz.breakdown.stats';
	my $string="<a href=\"$infilename.ps\"> <img src=\"$infilename.jpeg\" heigth=100 width=100 h/> </a>";
	push @transcriptfilesbreak,$string;
    }
    
    my $genomeps=$graphsdir.'/'.$prefix.'_unique_maps_genome.ps';
    my $genomejpeg=$graphsdir.'/'.$prefix.'_unique_maps_genome.jpeg';
    my $juncps=$graphsdir.'/'.$prefix.'_unique_maps_junctions.ps';
    my $juncjpeg=$graphsdir.'/'.$prefix.'_unique_maps_junctions.jpeg';
    my $transps=$graphsdir.'/'.$prefix.'_transcriptome_mapping.ps';
    my $transjpeg=$graphsdir.'/'.$prefix.'_transcriptome_mapping.jpeg';
    my $splitps=$graphsdir.'/'.$prefix.'_split_mapping.ps';
    my $splitjpeg=$graphsdir.'/'.$prefix.'_split_mapping.jpeg';
    # Print the html
    my $q=CGI->new();
    print $reportfh $q->start_html('RNAseq Report');
    print $reportfh $q->h1("Mapping stats for $project_id");
    # Print the general mapping stats
    print $reportfh table({-border=>undef},
			  caption('Unique mapping stats'),
			  Tr({-align=>'CENTER',-valign=>'CENTER'},
			     [
			      th(['Description', 'File']),
			      td(['Unique mappings to the genome' , "<a href=\"$genomeps\"> <img src=\"$genomejpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Unique mappings to the transcriptome (all)',
				 table({},td([@transcriptfilesall]))]),
			      td(['Unique mappings to the transcriptome (breakdown)',
				 table({},td([@transcriptfilesbreak]))]),
			      td(['Unique mappings to the junctions' , "<a href=\"$juncps\"> <img src=\"$juncjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Unique split maps' , "<a href=\"$splitps\"> <img src=\"$splitjpeg\" heigth=100 width=100 h/> </a>"])
			     ]
			  )
	);
    # Print the distribution of the reads across the transcripts

    print $reportfh $q->end_html();
    close($reportfh);
}

sub build_report_mapping {
    my $project_id=shift;
    my $prefix=shift;
    my $graphsdir=shift;
    my $reportreads='mapping_stats.html';

    my $reportfh=get_fh($reportreads,1);
    
    my $genomeps=$graphsdir.'/'.$prefix.'_genome_mapping.ps';
    my $genomejpeg=$graphsdir.'/'.$prefix.'_genome_mapping.jpeg';
    my $transps=$graphsdir.'/'.$prefix.'_transcriptome_mapping.ps';
    my $transjpeg=$graphsdir.'/'.$prefix.'_transcriptome_mapping.jpeg';
    my $juncps=$graphsdir.'/'.$prefix.'_junctions_mapping.ps';
    my $juncjpeg=$graphsdir.'/'.$prefix.'_junctions_mapping.jpeg';
    my $splitps=$graphsdir.'/'.$prefix.'_split_mapping.ps';
    my $splitjpeg=$graphsdir.'/'.$prefix.'_split_mapping.jpeg';
    # Print the html
    my $q=CGI->new();
    print $reportfh $q->start_html('RNAseq Report');
    print $reportfh $q->h1("Mapping stats for $project_id");
    # Print the general mapping stats
    print $reportfh table({-border=>undef},
			  caption('General Mapping stats'),
			  Tr({-align=>'CENTER',-valign=>'CENTER'},
			     [
			      th(['Description', 'File']),
			      td(['Genome Mapping' , "<a href=\"$genomeps\"> <img src=\"$genomejpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Transcriptome Mapping' , "<a href=\"$transps\"> <img src=\"$transjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Junctions Mapping' , "<a href=\"$juncps\"> <img src=\"$juncjpeg\" heigth=100 width=100 h/> </a>"]),
			      td(['Split Mapping' , "<a href=\"$splitps\"> <img src=\"$splitjpeg\" heigth=100 width=100 h/> </a>"])
			     ]
			  )
	);
    # Print the distribution of the reads across the transcripts

    print $reportfh $q->end_html();
    close($reportfh);
}

sub build_report_reads {
    my $project_id=shift;
    my $prefix=shift;
    my $graphsdir=shift;
    my $reportreads='read_stats.html';

    my $reportfh=get_fh($reportreads,1);
    
    my $readstatsps=$graphsdir.'/'.$prefix.'_read_stats.ps';
    my $readstatsjpeg=$graphsdir.'/'.$prefix.'_read_stats.jpeg';
    # Print the html
    my $q=CGI->new();
    print $reportfh $q->start_html('RNAseq Report');
    print $reportfh $q->h1("Read stats for $project_id");
    print $reportfh table({-border=>undef},
			  Tr({-align=>'CENTER',-valign=>'CENTER'},
			     [
			      th(['Description', 'File']),
			      td(['Read Stats' , "<a href=\"$readstatsps\"> <img src=\"$readstatsjpeg\" heigth=100 width=100 h/> </a>"])
			     ]
			  )
	);
    
    print $reportfh $q->end_html();
    close($reportfh);

}

sub build_report_main {
    my $project_id=shift;
    my $reportfn=shift;

    my $reportfh=get_fh($reportfn,1);
    
    # Print the html
    my $q=CGI->new();
    print $reportfh $q->start_html('RNAseq Report');
    print $reportfh $q->h1("RNAseq Report for $project_id");
    print $reportfh $q->li("<a href=\"read_stats.html\"> Read stats </a>");
    print $reportfh $q->li("<a href=\"mapping_stats.html\"> General mapping stats </a>");
    print $reportfh $q->li("<a href=\"unique_mappings.html\"> Unique mappings </a>");
    print $reportfh $q->li("<a href=\"read_clusters.html\"> Read clusters </a>");
    print $reportfh $q->li("<a href=\"expression_info.html\"> Expression information</a>");
    print $reportfh $q->end_html();
    close($reportfh);
}

sub get_transcript_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	$lanes{$files->{$file}->[1]}++;
    }

    return(\%lanes);
}
