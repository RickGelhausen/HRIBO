Minimal mapping, QC, and tracks tutorial
========================================

This tutorial runs one single-end Ribo-seq library through trimming, mapping,
quality control, and coverage-track generation.  It deliberately does not run
ORF prediction or differential expression.  Use a separate software checkout
and analysis directory so workflow updates do not mix with data or results.

Set up the checkout and project
-------------------------------

The commands below target the current development documentation.  Use a
released 2.0 tag instead of ``development`` when one is available.

.. code-block:: console

   $ tutorial_root="$PWD/hribo-minimal-tutorial"
   $ hribo_checkout="$tutorial_root/software/HRIBO"
   $ project_dir="$tutorial_root/analysis"
   $ mkdir -p "$tutorial_root/software" "$project_dir/config" \
       "$project_dir/data" "$project_dir/fastq"
   $ git clone --branch development --single-branch \
       https://github.com/RickGelhausen/HRIBO.git "$hribo_checkout"
   $ micromamba create --name hribo \
       --file "$hribo_checkout/environment.linux-64.pin.txt"
   $ micromamba activate hribo
   $ cp "$hribo_checkout/config/config.yaml" "$project_dir/config/"
   $ cp "$hribo_checkout/config/samples.tsv" "$project_dir/config/"

Copy a nucleotide reference FASTA, its matching GFF3 or GTF annotation, and a
gzip-compressed FASTQ into the project.  Substitute the three source paths
before executing these commands:

.. code-block:: console

   $ genome_source=/absolute/path/to/genome.fa
   $ annotation_source=/absolute/path/to/annotation.gff3
   $ fastq_source=/absolute/path/to/ribo-control-1.fastq.gz
   $ cp "$genome_source" "$project_dir/data/genome.fa"
   $ cp "$annotation_source" "$project_dir/data/annotation.gff3"
   $ cp "$fastq_source" "$project_dir/fastq/ribo-control-1.fastq.gz"

Configure the sample
--------------------

Replace ``config/samples.tsv`` with one tab-separated data row.  The final
``fastqFile2`` field is empty because this example is single-end:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2
   RIBO	Control	1	fastq/ribo-control-1.fastq.gz	""

Edit the copied ``config/config.yaml``.  Keep the settings not shown below from
the template, set the input paths relative to the project, enter the experiment's
real adapter sequence if one is present, and request only the three tutorial
stages:

.. code-block:: yaml

   biologySettings:
     adapterS3: ""
     adapterS5: ""
     genome: "data/genome.fa"
     annotation: "data/annotation.gff3"
     samples: "config/samples.tsv"

   workflowSettings:
     stages:
       - mapping
       - qc
       - tracks

An empty adapter value means that no adapter sequence is supplied; it is not a
generic value for an unknown adapter.  Resolve adapter uncertainty before a
production analysis.

Validate, then run
------------------

Construct the complete DAG without executing it:

.. code-block:: console

   $ cd "$project_dir"
   $ snakemake \
       --snakefile "$hribo_checkout/workflow/Snakefile" \
       --directory "$project_dir" \
       --configfile "$project_dir/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 8 \
       --resources reparation_instances=1 \
       --dry-run all

Fix every preflight error before removing ``--dry-run``.  Then run with safe
restart and diagnostic options:

.. code-block:: console

   $ snakemake \
       --snakefile "$hribo_checkout/workflow/Snakefile" \
       --directory "$project_dir" \
       --configfile "$project_dir/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 8 \
       --resources reparation_instances=1 \
       --rerun-incomplete \
       --latency-wait 60 \
       --printshellcmds \
       --show-failed-logs all

The first run creates rule environments and therefore requires network access.
After completion, inspect:

* ``maplink/RIBO-Control-1.bam`` and its index;
* ``qc/multi/multiqc_report.html``; and
* the ``globaltracks/``, ``centeredtracks/``, ``fiveprimetracks/``, and
  ``threeprimetracks/`` directories.

Reverse-strand BigWig values are negative by convention so both strands can be
displayed around a shared zero baseline.  See :doc:`../outputs` before using
the normalized tracks quantitatively.

Finally, repeat the dry-run command.  An unchanged successful project should
report that there is nothing to be done.  If it schedules work, inspect which
input, code dependency, configuration value, or incomplete output changed.
