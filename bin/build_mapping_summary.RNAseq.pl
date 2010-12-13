#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should take the results of the initial evaluation of the reads
# as well as that of the initial mapping and elaborate a summary using Latex

use RNAseq_pipeline_settings3 qw(read_config_file);
use RNAseq_pipeline3 qw(get_fh);
use RNAseq_mapping_summary3 ('get_summary_order','get_summary_figures',
			     'get_figure_legends','build_latex_template',
			     'build_postscript','get_transcript_figures',
			     'get_summary_tables','get_table_legends');

use Getopt::Long;
my $lowmem=0;
GetOptions('lowmem' => \$lowmem);

# Get a log file to print what is being done
my $log_fh;
# Redirect the log to a filehandle
$log_fh=*STDERR; # for debugging
#$log_fh=get_fh('build_mapping_summary.RNAseq.log',1);

# Read the configuration file
print $log_fh 'Reading the configuration file';
my %options=%{read_config_file()};

# Get the location of the figures to include
my $graph_path=$options{'PROJECT'}.'/graphs/';
print $log_fh "The path to the files is set to: $graph_path\n";

# Set the figure order
print $log_fh "Getting the order of the files in the summary...\n";
my @order=get_summary_order();
print $log_fh join("\t",
		   'Order:',
		   @order),"\n";

# Get the main figures
my %figures=%{get_summary_figures($graph_path,
				  \@order,
				  \%options,
				  $log_fh)};

my %tables=%{get_summary_tables(\@order,
				\%options,
				$log_fh)};

# Get the figures for the read distribution along the trasncripts
my %transcript_figs=%{get_transcript_figures($graph_path,
					     \%options,
					     $log_fh)};

# Get the legends to include
my %legends=%{get_figure_legends(\%figures,
				 $log_fh)};
my %tablegends=%{get_table_legends(\%tables,
				    $log_fh)};

# First generate the latex template
my $template=build_latex_template(\%figures,
				  \%tables,
				  \%legends,
				  \%tablegends,
				  \@order,
				  \%transcript_figs);

# Build the ps file
build_postscript($template,
		 $options{'PREFIX'});

exit;



