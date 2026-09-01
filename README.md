<img src="HRIBO.png" width="620">

# High-throughput annotation by Ribo-seq

[![GitHub](https://img.shields.io/github/tag/RickGelhausen/HRIBO.svg)](https://github.com/RickGelhausen/HRIBO)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Documentation Status](https://readthedocs.org/projects/hribo/badge/?version=latest)](http://hribo.readthedocs.io/?badge=latest)
[![PyPI Latest Release](https://img.shields.io/pypi/v/hribo.svg)](https://pypi.org/project/hribo/)

We present HRIBO (High-throughput annotation by Ribo-seq), a workflow to enable reproducible and high-throughput analysis of bacterial Ribo-seq data. The workflow performs all required pre-processing steps and quality control.  Importantly, HRIBO outputs annotation-independent ORF predictions based on two complementary prokaryotic-focused tools, and integrates them with additional computed features. This facilitates both the rapid discovery of ORFs and their prioritization for functional characterization.

For a detailed description of this workflow, the installation, usage and examples, please refer to the [ReadTheDocs documentation](http://hribo.readthedocs.io/?badge=latest).

HRIBO installs all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda create -c bioconda -c conda-forge -n snakemake snakemake

         source activate snakemake

### <u>Basic usage</u>

The retrieval of input files and running the workflow locally and on a server cluster via a queuing system is
working as follows. Create a project directory and change into it:

         mkdir project
         cd project

Retrieve the HRIBO from GitHub:

         git clone git@github.com:RickGelhausen/HRIBO.git

The workflow requires a genome sequence (fasta), an annotation file (gtf) and the sequencing results files (fastq).
We recommend retrieving both the genome and the annotation files from [Ensembl Genomes](http://ensemblgenomes.org/).
Copy the genome and the annotation file into the project folder, decompress them and name them genome.fa and annotation.gtf.

Create a folder fastq and copy your compressed fastq.gz files into the fastq folder.

Please copy the template of the sample sheet and the config file into a `config` folder
in your project directory (not into the HRIBO clone, so that the clone stays clean):

         mkdir -p config
         cp HRIBO/config/config.yaml config/
         cp HRIBO/config/samples.tsv config/

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

         snakemake --sdm conda apptainer -s HRIBO/workflow/Snakefile --directory ${PWD} -j 20 --latency-wait 60


Run Snakemake on the cluster:

Snakemake 8 and later use executor plugins rather than `--cluster`. Edit the bundled
SLURM profile (`HRIBO/workflow/profiles/slurm/config.yaml`) to set your account and
partition, then run:

       snakemake -s HRIBO/workflow/Snakefile --directory ${PWD} --profile HRIBO/workflow/profiles/slurm

This requires `snakemake-executor-plugin-slurm` in your Snakemake environment.

Run only parts of the workflow:

The `workflowSettings.stages` list in `config/config.yaml` decides what a run produces.
Comment out what you do not need; everything a remaining stage depends on is still built,
so asking only for `predictions` still trims, filters and maps the reads on the way there.

A single run can override the list without editing the config file:

       snakemake ... --config stages=mapping             # stop at the BAM files
       snakemake ... --config stages=mapping,tracks,qc   # BAM files, coverage tracks and QC
       snakemake ... --config stages=preprocessing       # trimming, mapping and the QC report

Once the workflow has finished you can request a automatically generated report.html file with the following command:

       snakemake -s HRIBO/workflow/Snakefile --directory ${PWD} --report report.html

