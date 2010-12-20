= RNASeq Pipeline =

= Dependencies =

* Bio::DB::Fasta

* Bundle::BioPerl

== Installation ==

# Create a new directory
mkdir 009TR.PE40
cd 009TR.PE40/

# Make a directory which will contain the read data and copy the read files in gzipped fasta or fastq format into this directory
mkdir readData
mv ../CLL.PE40/readData/009TR.r1.fastq.gz readData/
mv ../CLL.PE40/readData/009TR.r2.fastq.gz readData/

# Get the scripts form the repository 
svn co svn+ssh://dgonzalez@svn.crg.es/data/svn/big/pipeline/trunk .
A    lib
A    lib/RNAseq_GEM3.pm
...
A    bin
A    bin/build_top_transcripts.RNAseq.pl
...
A    template3.0.txt
A    Readme.txt

= .my.cnf =

Put the database connection information into your

    /Users/maik/.my.cnf

Example:

[client]
host=localhost
user=myuser
password=mypassword

= Setting up a demo database =

RNAseqPipeline
RNAseqPipelineCommon

= Getting help on the start script =

{{{
bin/start_RNAseq_pipeline.3.0.pl --help
$ bin/start_RNAseq_pipeline.3.0.pl --help
Usage:
        start_pair_pipeline.3.0.pl -species ... -genome ... -annotation ... -project ... -experiment ... -template ... -readlength...
    
      Help options:
        -help:           brief help message
        -man:            full documentation
        -debug:
      Mandatory options:
        -species:        Species for which the pipeline is run.
        -genome:         File with the genomic sequence.
        -annotation:     File with the annotation to use.
        -project:        The project to which the experiment will be added.
        -experiment:     The set of reads to be added.
        -template:       File containing the commands that will be executed.
        -readlength:     Nucleotide length of the reads.
      Mapping Options:
        -mapper:         Mapper to be used.
                          Defaults to GEM which is the only one supported currently
        -mismatches:     Number of mismatches with which the mapping will be done.
                          Default 2.
        -paired:         Reads are paired.
        -stranded:       Reads are stranded.
        -threads:        Number of threads to use for mapping.
                          (if mapper allows multiple threads)
    
      Advanced Options:
        -database:       Sets the database to use for the experiment tables.
        -commondb:       Sets the database to use for the common tables.
        -host:           Sets the host where the databases are located.
        -localdir:       Directory in which to store the temporary files generated during the pipeline execution. It is recommandable to set it to a local drive.

Options:
    -help
                Print a brief help message and exits.
    
    -man
                Prints the manual page and exits.
}}}

= Starting on a demo server =

{{{
perl bin/start_RNAseq_pipeline.3.0.pl -species 'Homo sapiens' -genome /Users/maik/Pipeline/ReferenceGenomes/H.sapiens/H.sapiens.genome.hg19.sanger.fa -annotation /Users/maik/Pipeline/ReferenceAnnotations/H.sapiens/H.sapiens.EnsEMBL.55.parsed.gtf -project CLL -experiment 009TR40 -template template3.0.txt -readlength 40 -threads 2 -database RNAseqPipeline -commondb RNAseqPipelineCommon -host localhost 
}}}

In the database, you should now have the following tables:

* annotation_files

* annotation_tables

* exclusion_files

* experiments

* fasta_files

* genome_files

* indices

* mappabilities

* projects

* species_info

= Starting on the production server =

# Run the start script with wahtever options are needed
bin/start_RNAseq_pipeline.3.0.pl -species 'Homo sapiens' -genome /users/rg/dgonzalez/ReferenceGenomes/H.sapiens/H.sapiens.genome.hg19.sanger.fa -annotation /users/rg/dgonzalez/ReferenceAnnotations/H.sapiens/H.sapiens.EnsEMBL.55.parsed.gtf -project CLL -experiment 009TR40 -template template3.0.txt -readlength 40 -threads 2

# Edit the read.list.txt file if necessary to add the pair Id and the read Id separated by tabs and in the fourth field optionally a tag that will be used to group different samples
vi read.list.txt
Paired
001TR.r1.fastq	001TR	001TR.1 Tumor
001TR.r2.fastq	001TR	001TR.2	Tumor
Single (pair id = read_id)
001TR.r1.fastq  001TR.1 001TR.1 Tumor
001TR.r2.fastq  001TR.2	001TR.2 Tumor

# Make sure the soft links to the flux.sh and overlap (if the binary is not there) are in the bin directory or the binaries are in the path (also in the cluster)
http://big.crg.cat/services_and_software
overlap: ln -s /users/rg/sdjebali/bin/overlap bin
flux.sh: ln -s /users/rg/dgonzalez/bin/flux.sh bin

# Run the execute script with the step up to which the pipeline needs to run
bin/execute_RNAseq_pipeline3.0.pl all |tee -a pipeline.log
