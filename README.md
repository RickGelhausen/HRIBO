# SPtools

SPTools can be installed with all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -c bioconda -c conda-forge -n SPTools snakemake 
         
         source activate SPTools

### <u>Basic usage</u>

The retrieval of input files and running the workflow locally and on a server cluster via a queuing system is
working as follows. Create a project directory and change into it:

         mkdir project
         cd project

Retrieve the SPTools from GitHub:

         git clone git@github.com:gelhausr/SPTools.git

The workflow requires a genome sequence (fasta), an annotation file (gtf) and the sequencing results files (fastq).
We recommend retrieving both the genome and the annotation files from [Ensembl Genomes](http://ensemblgenomes.org/).
Copy the genome and the annotation file into the project folder, decompress them and name them genome.fa and annotation.gtf.

Create a folder fastq and copy your compressed fastq.gz files into the fastq folder.

Please copy the template of the sample sheet and the config file into the SPTools folder.

         cp SPTools/templates/config.yaml SPTools/
         cp SPTools/templates/samples.tsv SPTools/
       
Customize the config.yaml with the used adapter sequence and optionally with the path to a precomputed
STAR genome index. For correct removal of reads mapping to ribosomal genes please specify the taxonomic group of
the used organism (Eukarya, Bacteria, Archea).
Now edit the sample sheet corresponding to your project, using one line per sequencing result, stating the used
method (RIBO for ribosome profiling, RNA for RNA-seq), the applied condition (e.g. A, B, CTRL, TREAT), the replicate (e.g. 1, 2,..) and the filename. Following is an example:

|method|	condition |replicate|	fastqFile                 |
|------|-----------|---------|--------------------------------|
|RIBO  |	A         |        1|"fastq/FP-ctrl-1-2.fastq.gz"    |
|RIBO  |	B         |        1|"fastq/FP-treat-1-2.fastq.gz"   |
|RNA   |	A         |        1|"fastq/Total-ctrl-1-2.fastq.gz" |
|RNA   |	B         |        1|"fastq/Total-treat-1-2.fastq.gz"|

Now you can start your workflow.

Run Snakemake locally:

         snakemake --use-conda -s Snakefile --configfile SPTools/config.yaml --directory ${PWD} -j 20 --latency-wait 60 
         

Run Snakemake on the cluster:

Edit cluster.yaml according to your queuing system and cluster hardware. The following example works for Grid Engine:

       snakemake --use-conda -s Snakefile --configfile SPTools/config.yaml --directory ${PWD} -j 20 --cluster-config SPTools/cluster.yaml --cluster "qsub -N {cluster.jobname} -cwd -q {cluster.qname} -pe {cluster.parallelenvironment} -l {cluster.memory} -o {cluster.logoutputdir} -e {cluster.erroroutputdir} -j {cluster.joinlogs} -M <email>" --latency-wait 60 

Once the workflow has finished you can request a automatically generated report.html file with the following command:
         
         snakemake --report report.html

