package RNAseq_mapping_summary3;
# DGK 2008-2010 CRG

# Export subroutines to caller namespace
# Must be done before strict is used
use Exporter;
@ISA=('Exporter');

push @EXPORT_OK,('get_summary_order','get_summary_figures','get_figure_legends',
		 'build_latex_template','build_postscript','get_table_legends',
		 'get_transcript_figures','get_summary_tables');

use strict;
use warnings;

# Import subs from other modules
use RNAseq_pipeline3 ('get_fh','MySQL_DB_Connect');
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list',
			       'get_lanes_single','get_lanes_paired');

# Set the order in which the figures and tables will be added to the document
sub get_summary_order {
    my @order=('read_stats',
	       'genome_mapping',
	       'unique_maps_genome',
	       'transcriptome_mapping',
	       'junctions_mapping',
	       'unique_maps_junctions',
	       'split_mapping',
	       'unique_maps_split',
	       'clusters',
#	       'gene_detection',
#	       'exon_detection',
#	       'junction_detection'
	);
    return(@order);
}

sub get_summary_figures {
    my $path=shift;
    my $list=shift;
    my $options=shift;
    my $log_fh=shift;
    my %figures;
    my $prefix=$options->{'PREFIX'};

    # This is a description of all the figures and their locations
    my %all=('read_stats' => $path.$prefix.'_read_stats.ps',
	     'genome_mapping' => [$path.$prefix.'_genome_mapping.ps',
				  $path.$prefix.'_unique_maps_genome.ps'],
	     'transcriptome_mapping' => $path.$prefix.'_transcriptome_mapping.ps',
	     'junctions_mapping' => [$path.$prefix.'_junctions_mapping.ps',
				     $path.$prefix.'_unique_maps_junctions.ps'],
	     'split_mapping' => [$path.$prefix.'_split_mapping.ps',
				 $path.$prefix.'_unique_maps_split.ps'],
	     'clusters' => $path.$prefix.'_initial_clusters.ps',
	     'gene_detection' => [$path.$prefix.'.Genes.detection.ps',
				  $path.'Gene.saturation.ps'],
	     'exon_detection' => [$path.$prefix.'.Exons.detection.ps',
				  $path.'Exon.saturation.ps'],
	     'junction_detection' => [$path.$prefix.'.Junctions.detection.ps',
				  $path.'Junction.saturation.ps']
	);
    
    # Select the ones we are interested in
    my @list2;
    foreach my $fig (@{$list}) {
	if ($all{$fig}) {
	    $figures{$fig}=$all{$fig};
	    push @list2,$fig;
	} else {
	    print $log_fh "I have no figures for $fig. So I will exclude it\n";
	}
    }
#    @{$list}=@list2;

    return(\%figures);
}

sub get_summary_tables {
    my $list=shift;
    my $options=shift;
    my $log_fh=shift;
    my %tables;
    my $prefix=$options->{'PREFIX'};

    my $dbh=MySQL_DB_Connect($options->{'DB'},
			     $options->{'HOST'});

    my %table_subs=%{get_table_subs($dbh,
				  $prefix)};

    # This is a description of all the figures and their locations
    my %all=('read_stats' => [$prefix.'_read_stats'],
	     'genome_mapping' => [$prefix.'_genome_mapping'],
	     'unique_maps_genome' => [$prefix.'_unique_maps_genome'],
	     'transcriptome_mapping' => [$prefix.'_transcriptome_mapping'],
	     'junctions_mapping' => [$prefix.'_junctions_mapping'],
	     'unique_maps_junctions' => [$prefix.'_unique_maps_junctions'],
	     'split_mapping' => [$prefix.'_split_mapping'],
	     'unique_maps_split' => [$prefix.'_unique_maps_split'],
	     'clusters' => [$prefix.'_initial_clusters.ps'],
	     'gene_detection' => [$prefix.'.Genes.detection.ps',
				  'Gene.saturation.ps'],
	     'exon_detection' => [$prefix.'.Exons.detection.ps',
				  'Exon.saturation.ps'],
	     'junction_detection' => [$prefix.'.Junctions.detection.ps',
				  'Junction.saturation.ps']
	);
    
    # Select the ones we are interested in
    my @list2;
    foreach my $fig (@{$list}) {
	if ($all{$fig}) {
	    $tables{$fig}=get_table_content($all{$fig},
					    \%table_subs);
	    push @list2,$fig;
	} else {
	    print $log_fh "I have no table for $fig. So I will exclude it\n";
	}
    }
#    @{$list}=@list2;

    return(\%tables);
}

# Get the table information from the database as a tes of arrays, one per line
sub get_table_subs {
    my $dbh=shift;
    my $prefix=shift;

    my %table_subs;
    $table_subs{$prefix.'_read_stats'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT LaneName, ReadLength, TotalReads, NoAmbiguousBases, ';
	$query.='AmbiguousBases, UniqueReads ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Length','Total','Unambiguous','Ambiguous','Unique');
	push @table, [@header];;
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_genome_mapping'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT LaneName, totalReads, mappedReads, ';
	$query.='uniqueReads, 100uniqueReads ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Total','Mapped','Unique','100 Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_unique_maps_genome'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT filename, chromosome, uniqueMaps ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Chrom','Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_transcriptome_mapping'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT LaneName, totalReads, mappedReads, ';
	$query.='uniqueReads, 100uniqueReads ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Total','Mapped','Unique','100 Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_junctions_mapping'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT LaneName, totalReads, mappedReads, ';
	$query.='uniqueReads, 100uniqueReads ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Total','Mapped','Unique','100 Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_unique_maps_junctions'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT filename, type, number ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Type','Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_split_mapping'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT LaneName, totalReads, mappedReads, ';
	$query.='uniqueReads, 100uniqueReads ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Total','Mapped','Unique','100 Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };
    $table_subs{$prefix.'_unique_maps_split'} = sub {
	my $table_id=shift;

	print STDERR $table_id,"\n";
	my @table;
	my ($query,$sth);
	$query ='SELECT filename, type, number ';
	$query.="FROM $table_id\n";
	$sth=$dbh->prepare($query);
	$sth->execute();
	
	my @header=('Lane','Type','Unique');
	push @table, [@header];
	while (my @line=$sth->fetchrow_array()) {
	    $line[1]=~s/&//g;
	    push @table, [@line];
	}
	$sth->finish();
	return(@table);
    };


    return(\%table_subs);
}
sub get_table_content {
    my $tableids=shift;
    my $table_subs=shift;

    my @tables;
    my ($query,$sth);

    foreach my $table_id (@{$tableids} ) {
	print STDERR $table_id,,"\n";
	if ($table_subs->{$table_id}) {
	    my @table=$table_subs->{$table_id}->($table_id);
	    push @tables, [@table];
	}
    }
    return(\@tables);
}


# Get the transcriptome figures for the distribution of the reads along the
# different transcripts
sub get_transcript_figures {
    my $path=shift;
    my $options=shift;
    my $log_fh=shift;

    my $prefix=$options->{'PREFIX'};

    my %files=%{read_file_list()};
    my %lanes=%{get_lanes_single(\%files)};
    my %pairs=%{get_lanes_paired(\%files)};
    my $type;

    foreach my $lane (keys %pairs) {
	if ($pairs{$lane} == 1) {
	    $type='single';
	} elsif ($pairs{$lane} == 2) {
	    $type='paired';
	} else {
	    die "Unknown type\n";
	}
    }
    my %figures;

    # If the reads are single normally the read and pair id should be the same,
    # but this is not necessarilly so right now so we will check to be sure
    if ($type eq 'single') {
	foreach my $fig (keys %lanes) {
	    $figures{$fig}=[$path.$fig.'.'.$type.'.unique.gtf.gz.all.stats.ps',
			    $path.$fig.'.'.$type.'.unique.gtf.gz.breakdown.stats.ps'];
	}
    } elsif ($type eq 'paired') {
	foreach my $fig (keys %pairs) {
	    $figures{$fig}=[$path.$fig.'.'.$type.'.unique.gtf.gz.all.stats.ps',
			    $path.$fig.'.'.$type.'.unique.gtf.gz.breakdown.stats.ps'];
	}
    }

    # Check the figures
    foreach my $fig (keys %figures) {
	unless ($figures{$fig}->[0]) {
	    die "Missing $fig all.stats\n";
	}
	unless ($figures{$fig}->[1]) {
	    die "Missing $fig breakdown.stats\n";
	}
    }

    return(\%figures);
}

sub get_figure_legends {
    my $figures=shift;
    my $log_fh=shift;

    # Here are all the legends for the figures
    my %all=('read_stats' => ['Read analysis summary',
			      'This figure shows some summary information for the reads extracted before mapping. 1: Total number of reads in each lanes as well as the number of reads showing ambiguous nucleotides, no ambiguous nucleotides and the number of unique sequences present in the lane. 2: Fraction of reads with no ambiguous nucleotides and with ambiguous nucleotides. 3: Average quality at each position of the reads, if no qualities are provided the graph shows a flat line and quality 50. 4: Number of Ns at each position of the read.'],
	     'genome_mapping1' => ['Genome mapping summary',
				  'This figure shows the reads mapping to the genome. 1: Total number of reads in each lane and number mapped to the genome, as well as those uniquely mapped and those mapped uniquely qith no similar hit withup to two more mismatches. 2: Fraction of reads in each category for each lane. 3 Distribution of the uniquely mapped reads along the different chromosomes'],
	     'transcriptome_mapping' => ['Transcriptome mapping summary',
					 'This figure shows the reads mapping to the transcriptome. 1: Total number of reads in each lane and number mapped to the genome, as well as those uniquely mapped and those mapped uniquely qith no similar hit withup to two more mismatches. 2: Fraction of reads in each category for each lane'],
	     'junctions_mapping' => ['Junctions mapping summary',
				     'This figure shows the reads mapping to the junctions. 1: Total number of reads in each lane and number mapped to the genome, as well as those uniquely mapped and those mapped uniquely qith no similar hit withup to two more mismatches. 2: Fraction of reads in each category for each lane. 3: FRaction of unique reads falling on known or novel junctions, where novel junctions are defined as those junctions formed by combinations of exons not present in the annotation as spliced'],
	     'split_mapping' => ['Split mapping summary',
				 'This figure shows the results of mapping those reads that could not be mapped to genome transcriptome or junctions using the GEM split mapper. 1 Upper: Total number of reads in each lane and number mapped to the genome, as well as those uniquely mapped and those mapped uniquely qith no similar hit withup to two more mismatches. 2 Upper: Fraction of reads in each category for each lane. 3 Lower: Distribution of the unique mappings according to the location of the two split halves on the genome; Different chromosome, Same chromosome different strand, same chromososme same strand, but with the first part of the read mappping downstream of from the second, smae chromosoem normal orientation'],
	     'clusters' => ['Cluster summary','This figure shows the distribution along the different chromosomes of the clusters of overlapping uniquely mapped reads from the mappings to the genome, the junctions and the split mapping'],
	     'gene_detection' => ['Gene detection summary','This figure shows the number of genes detected in each of the lanes by at least one uniquely mapped read, as well as the saturation curve for the gene detection when combining the different lanes. 1: Gene detection. 2: Gene saturation'],
	     'exon_detection' => ['Exon detection summary','This figure shows the number of exons detected in each of the lanes by at least one uniquely mapped read, as well as the saturation curve for the exon detection when combining the different lanes. 1: Exon detection. 2: Exon saturation'],
	     'junction_detection' => ['Junction detection summary','This figure shows the number of junctions detected in each of the lanes by at least one uniquely mapped read, as well as the saturation curve for the junction detection when combining the different lanes. 1: Junction detection separated by read. 2: Junction saturation'],
);

    # We will return those legends that are requested
    my %legends;
    foreach my $figure (keys %{$figures}) {
	if ($all{$figure}) {
	    $legends{$figure}=$all{$figure};
	} else {
	    print $log_fh "I have no legend for $figure. Setting to -\n";
	    my $text=$figure;
	    $text=~s/_/\\_/g;
	    $legends{$figure}=[$text,
			       'Figure ?'];
	}
    }

    return(\%legends);
}

sub get_table_legends {
    my $tables=shift;
    my $log_fh=shift;

    # Here are all the legends for the figures
    my %all=('read_stats' => ['Read summary',
			      'This table shows some information for the reads extracted before mapping, such as read length, total number of reads and the number of reads showing ambiguous bases. The last column also shows the number of unique reads that are present in the lane'],
	     'genome_mapping' => ['Genome mapping','This table shows the mapping statistics'],
	     'unique_maps_genome' => ['Chromosome distribution',
				  'This table shows the distribution of the uniquely mapping reads along the different chromosomes'],
	     'transcriptome_mapping' => ['Transcriptome mapping summary',
					 'This figure shows the reads mapping to the transcriptome. 1: Total number of reads in each lane and number mapped to the genome, as well as those uniquely mapped and those mapped uniquely qith no similar hit withup to two more mismatches. 2: Fraction of reads in each category for each lane'],
	     'junctions_mapping' => ['Junctions mapping summary',
				     'This table shows a summary of the reads mapping to the junctions'],
	     'unique_maps_junctions' => ['Novel junctions','This table shows thenumber of known and novel junctios that were identified by uniquely mapping reads'],
	     'split_mapping' => ['Split mapping summary',
				 'This table shows the results of mapping those reads that could not be mapped to genome transcriptome or junctions using the GEM split mapper.'],
	     'unique_maps_split' => ['Type of split maps','This is a suumary of the different types of split maps encountered'],
	     'clusters' => ['Cluster summary','This figure shows the distribution along the different chromosomes of the clusters of overlapping uniquely mapped reads from the mappings to the genome, the junctions and the split mapping'],
	     'gene_detection' => ['Gene detection summary','This figure shows the number of genes detected in each of the lanes by at least one uniquely mapped read, as well as the saturation curve for the gene detection when combining the different lanes. 1: Gene detection. 2: Gene saturation'],
	     'exon_detection' => ['Exon detection summary','This figure shows the number of exons detected in each of the lanes by at least one uniquely mapped read, as well as the saturation curve for the exon detection when combining the different lanes. 1: Exon detection. 2: Exon saturation'],
	     'junction_detection' => ['Junction detection summary','This figure shows the number of junctions detected in each of the lanes by at least one uniquely mapped read, as well as the saturation curve for the junction detection when combining the different lanes. 1: Junction detection separated by read. 2: Junction saturation'],
);

    # We will return those legends that are requested
    my %legends;
    foreach my $table (keys %{$tables}) {
	if ($all{$table}) {
	    $legends{$table}=$all{$table};
	} else {
	    print $log_fh "I have no legend for $table. Setting to -\n";
	    my $text=$table;
	    $text=~s/_/\\_/g;
	    $legends{$table}=[$text,
			       'Table ?: No information available'];
	}
    }

    return(\%legends);
}

## These subs will build the latex template and run it
# First getting the template 
sub build_latex_template {
    my $figures=shift;
    my $tables=shift;
    my $legend=shift;
    my $tablegends=shift;
    my $order=shift;
    my $transcript_figs=shift;
    my $lowmem=shift;

    my $templatefn="tmp.$$";
    my $templatefh=get_fh($templatefn.'.tex',1);
    build_latex_header($templatefh);
    # Insert the mapping info
    foreach my $fig (@{$order}) {
	if (exists $tables->{$fig}) {
	    build_latex_table($templatefh,
			      $tables->{$fig},
			      $tablegends->{$fig});
	}
	if ($lowmem) {
	    print $templatefh "\\clearpage\n";
	}
	if (exists $figures->{$fig}) {
	    build_latex_figure($templatefh,
			       $figures->{$fig},
			       $legend->{$fig});
	}
	if ($lowmem) {
	    print $templatefh "\\clearpage\n";
	}
    }
    print $templatefh "\\clearpage\n";

    # Insert the transcript read distribution info
    if ($transcript_figs) {
	foreach my $fig (sort keys %{$transcript_figs}) {
	    build_latex_trans_fig($templatefh,
				  $transcript_figs->{$fig},
				  $fig);
	}
	print $templatefh "\\clearpage\n";
    }

    # Junction mappings: Known and novel
    
    # Gene and exons distribution

    # Splitmapping 

    build_latex_doc_end($templatefh);
    return ($templatefn);
}

# Insert a table into the latex document
sub build_latex_table {
    my $fh=shift;
    my $tables=shift;
    my $legend=shift;
    my @legend=@{$legend};
    # Go through each of the tables
    foreach my $table (@{$tables}) {
	my $header=shift(@{$table});
	my @columns;
	foreach my $col (@{$header}) {
	    push @columns,'c';
	}
	my $columns=join('|',@columns);
	print $fh "\\begin{table}[ht]\n";
	print $fh "\\caption{\\label{read summary} ",$legend[0],"}\n";
	print $fh "\\begin{center}\n";
	print $fh "\\begin{tabular}{|$columns|}\n";

	# Print the table header
	my $head_line=join(' & ',@{$header});
	$head_line=~s/_/\\_/g;
	print $fh '\hline ',$head_line,' \\\\ \hline \hline',"\n";
	while (my $row=shift(@{$table})) {
	    my $line=join(' & ',@{$row});
	    $line=~s/_/\\_/g;
	    print $fh $line,' \\\\ \hline',"\n";
	}
	print $fh "\\end{tabular}\n";

	print $fh "\\parbox[b]{12 cm} {\\footnotesize $legend[1]}\n";
	print $fh "\\end{center}\n";
	print $fh "\\end{table}\n";
    }
 
}


# Insert a figure into the latex document
sub build_latex_figure {
    my $fh=shift;
    my $figure=shift;
    my $legend=shift;
    my @legend=@{$legend};
    if (ref $figure) {
	print $fh <<"FIGEND";
\\begin{figure}[h]
\\begin{minipage}{\\textwidth}
\\begin{center}
\\includegraphics[height=4.8in,width=3.5in,angle=-90]{$figure->[0]}
\\includegraphics[height=4.8in,width=3.5in,angle=-90]{$figure->[1]}
\\caption{\\label{read summary} $legend[0]}
\\parbox[b]{12 cm} {\\footnotesize $legend[1]}
\\end{center}
\\end{minipage}
\\end{figure}
FIGEND

    } else {
    print $fh <<"FIGEND";

\\begin{figure}[h]
\\begin{minipage}{\\textwidth}
\\begin{center}
\\includegraphics[height=7in,width=5.0in,angle=0]{$figure}
\\caption{\\label{read summary} $legend[0]}
\\parbox[b]{12 cm} {\\footnotesize $legend[1]}
\\end{center}
\\end{minipage}
\\end{figure}
FIGEND
    }
}

# Build the transcript distribution figures
sub build_latex_trans_fig {
    my $fh=shift;
    my $figure=shift;
    my $lane=shift;
    $lane=~s/_/-/g;
    print $fh <<"FIGEND";

\\begin{figure}[h]
\\begin{minipage}{\\textwidth}
\\begin{center}
\\includegraphics[height=4.8in,width=3.5in,angle=-90]{$figure->[0]}
\\includegraphics[height=4.8in,width=3.5in,angle=-90]{$figure->[1]}
\\caption{\\label{read summary} $lane Read Distribution}
\\parbox[b]{12 cm} {\\footnotesize Distribution of reads over the transcript length in lane $lane. 1: Overall distribution. 2: Breakdown of the read distribution by transcript length.}
\\end{center}
\\end{minipage}
\\end{figure}
FIGEND
}

# build the latex header
sub build_latex_header {
    my $fh=shift;
    my %options=%{read_config_file()};
    my %files=%{read_file_list()};
    my %lanes=%{get_lanes_paired(\%files)};
    my $paired;

    foreach my $lane (keys %lanes) {
	if ($lanes{$lane} == 1) {
	    $paired='single';
	} elsif ($lanes{$lane} == 2) {
	    $paired='paired';
	} else {
	    die "Unknown type\n";
	}
    }

    print $fh <<'HEAD1END';
\documentclass[11pt]{article}
\usepackage{graphics}
\usepackage[dvips]{graphicx}

HEAD1END
    print $fh '\title{'."Summary of the RNAseq mapping for the X lanes ($paired) from $options{'SPECIES'}}\n";

    print $fh <<'HEAD2END';
\begin{document}
\pagestyle{empty}

\maketitle

\section{Overview}

The first part of this Report gives a summary of the read information as well as the mapping results of the dataset against genome transcriptome and junctions. Those reads that were not mapped to any of these were splitmapped using the GEM split mapper in an attempt to find possible unannotated exon junctions to which these unmapped reads may belong.

The second part of the report provides information regarding the mapping of the reads to the available annotation, showing the number of genes and exons that were found.

In order to identify possible biases in the Data the third part of the report shows the distribution of the read mappings along the transcripts. Both over the complete dataset and broken down into different transcript lengths and different expression levels.

Finally an estimation of the expression levels of the transcripts obtained with the Flux capacitor is shown as well as a measure of the stability

HEAD2END
}

# Build the latex document end
sub build_latex_doc_end {
    my $fh=shift;
    print $fh "\n",'\end{document}',"\n";
}

# Transform the latex template into postscript
sub build_postscript {
    my $file=shift;
    my $prefix=shift;
    my $command="latex $file.tex";
    system($command);
    $command="dvips -f $file.dvi > $prefix.report.ps";
    # Run twice to get all the figures correct
    system($command);
    system($command);

    # Clean up
    $command="rm $file.*";
    system($command);
}

1;
