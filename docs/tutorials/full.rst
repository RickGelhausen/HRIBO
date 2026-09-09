Full matched Ribo-seq/RNA-seq tutorial
======================================

This tutorial requests the ``full`` preset for two conditions, with two
Ribo-seq and two RNA-seq replicates per condition.  It includes ORF prediction,
differential expression, metagene profiling, the TIS advisor, and the combined
overview.  Complete :doc:`minimal` first if the checkout, project layout, or
launcher environment is unfamiliar.

.. warning::

   Representative real-data validation is still a release gate for
   |release|.  Biological screenshots and expected numerical or candidate
   results will be added only after that controlled validation.  The commands
   and output contracts below are tested, but they are not a claim that one
   particular biological signal must appear in every dataset.

Prepare a separate analysis directory
-------------------------------------

The commands below keep the HRIBO checkout outside the analysis.  They target
the current development documentation; use a released 2.0 tag when available.

.. code-block:: console

   $ tutorial_root="$PWD/hribo-full-tutorial"
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

Copy the matching genome, annotation, and eight gzip-compressed FASTQ files to
``data/`` and ``fastq/``.  Use biological replicates, not repeated sequencing
of the same library presented as independent replicates.

Define the matched design
-------------------------

Create ``config/samples.tsv`` with exact tab-separated headers and paths.  This
example uses single-end files:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2
   RIBO	Control	1	fastq/ribo-control-1.fastq.gz	""
   RIBO	Control	2	fastq/ribo-control-2.fastq.gz	""
   RIBO	Treated	1	fastq/ribo-treated-1.fastq.gz	""
   RIBO	Treated	2	fastq/ribo-treated-2.fastq.gz	""
   RNA	Control	1	fastq/rna-control-1.fastq.gz	""
   RNA	Control	2	fastq/rna-control-2.fastq.gz	""
   RNA	Treated	1	fastq/rna-treated-1.fastq.gz	""
   RNA	Treated	2	fastq/rna-treated-2.fastq.gz	""

Edit the copied configuration.  Preserve the complete template, set its input
and adapter values, request the ``full`` preset, and make the contrast direction
explicit:

.. code-block:: yaml

   biologySettings:
     genome: "data/genome.fa"
     annotation: "data/annotation.gff3"
     samples: "config/samples.tsv"

   differentialExpressionSettings:
     contrasts: ["Treated-Control"]

   predictionSettings:
     deepribo: "on"
     deepriboASiteOffset: 12

   workflowSettings:
     stages: "full"

Positive log2 fold changes for this contrast indicate higher signal in
``Treated``.  The DeepRibo offset shown is the template default derived from its
published E. coli setup; verify it for the organism and protocol rather than
copying the TIS advisor's P-site offset.  With DeepRibo enabled, the genome
FASTA may contain only uppercase ``A``, ``C``, ``G``, ``T``, and ``N``;
lowercase or other IUPAC ambiguity symbols fail preflight because DeepRibo
cannot encode them safely.  Review all metagene read lengths and filters in
:doc:`../configuration` before the production run.

Plan resources and downloads
----------------------------

The full preset includes every stage, including the expensive predictors and
differential engines.  A default REPARATION job requests 12 threads, 30 GB of
memory, and 30 GB of disk; DeepRibo jobs request up to 20 GB of memory.  The
first run also downloads rule environments, digest-pinned containers, the
DeepRibo model, and a checksum-verified Swiss-Prot release.  Confirm local or
cluster quotas before starting.

Limit concurrent REPARATION instances even when many cores are available.  On
SLURM, configure the bundled profile as described in :doc:`../troubleshooting`.

Dry-run and execute
-------------------

Run preflight and DAG construction first:

.. code-block:: console

   $ cd "$project_dir"
   $ snakemake \
       --snakefile "$hribo_checkout/workflow/Snakefile" \
       --directory "$project_dir" \
       --configfile "$project_dir/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 20 \
       --resources reparation_instances=1 \
       --dry-run all

After resolving every preflight finding classified as an error, execute the
same target:

.. code-block:: console

   $ snakemake \
       --snakefile "$hribo_checkout/workflow/Snakefile" \
       --directory "$project_dir" \
       --configfile "$project_dir/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 20 \
       --resources reparation_instances=1 \
       --rerun-incomplete \
       --latency-wait 60 \
       --printshellcmds \
       --show-failed-logs all

Review the analysis
-------------------

Begin with ``qc/multi/multiqc_report.html`` and the final BAMs, then review the
source tables behind summary plots.  The primary deliverables and normalization
semantics are catalogued in :doc:`../outputs`.  In particular, inspect:

* metagene start and stop profiles for the selected read lengths and both
  transcript strands;
* TIS-advisor evidence, confidence, and per-length offsets rather than only its
  headline recommendation;
* REPARATION and DeepRibo calls separately before the combined updated
  annotation;
* xTail, Riborex, and deltaTE contrast directions and adjusted p-values; and
* ``auxiliary/overview.xlsx`` together with its TSV and GFF3 companions.

Generate a native workflow report only after the requested targets complete:

.. code-block:: console

   $ snakemake \
       --snakefile "$hribo_checkout/workflow/Snakefile" \
       --directory "$project_dir" \
       --configfile "$project_dir/config/config.yaml" \
       --report "$project_dir/report.html"

Finish with an ``all`` dry-run.  It should be a no-op.  For an upgrade or
release decision, retain this result directory and compare it semantically with
the baseline using the protocol in :doc:`../real-data-validation`; do not infer
equivalence from successful execution alone.
