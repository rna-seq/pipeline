# Create a new directory
mkdir 009TR.PE40
cd 009TR.PE40/

# Make a directory which will contain the read data and copy the read files in gzipped fasta or fastq format into this directory
mkdir readData
mv ../CLL.PE40/readData/009TR.r1.fastq.gz readData/
mv ../CLL.PE40/readData/009TR.r2.fastq.gz readData/

# Get the scripts form the repository 
svn co svn+ssh://dgonzalez@svn.crg.es/data/svn/big/pipeline/trunk .

# Run the start script with wahtever options are needed
bin/start_RNAseq_pipeline.3.0.pl -species 'Homo sapiens' -genome /users/rg/dgonzalez/ReferenceGenomes/H.sapiens/H.sapiens.genome.hg19.sanger.fa -annotation /users/rg/dgonzalez/ReferenceAnnotations/H.sapiens/H.sapiens.EnsEMBL.55.parsed.gtf -project CLL -experiment 009TR40 -template template3.0.txt -readlength 40 -threads 2

# Edit the read.list.txt file if necessary to add the pair Id and the read Id separated by tabs and in the fourth field optionally a tag that will be used to group different samples
vi read.list.txt

# Make sure the soft links to the flux.sh and overlap (if the binary is not there) are in the bin directory or the binaries are in the path (also in the cluster)
overlap: /users/rg/sdjebali/bin/overlap
flux.sh: /users/rg/dgonzalez/bin/flux.sh

# Run the execute script with the step up to which the pipeline needs to run
bin/execute_RNAseq_pipeline3.0.pl all |tee -a pipeline.log
