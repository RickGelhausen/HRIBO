Recipe: complete matched analysis
=================================

This recipe requests every HRIBO stage for two conditions with two Ribo-seq
and two matched RNA-seq biological replicates per condition.  It produces QC,
coverage, metagene and TIS reports, ORF predictions, differential results, and
the combined overview.

Complete the smaller :doc:`minimal` recipe first if the checkout, project
layout, or launcher environment is unfamiliar.

Prepare the project
-------------------

Create a separate analysis directory and copy the current templates:

.. code-block:: console

   $ hribo_checkout=/path/to/HRIBO
   $ project_dir=/path/to/hribo-full-analysis
   $ mkdir -p "$project_dir"/{config,data,fastq}
   $ cp "$hribo_checkout/config/config.yaml" "$project_dir/config/"
   $ cp "$hribo_checkout/config/samples.tsv" "$project_dir/config/"

Place a matching genome and annotation in ``data/`` and the eight
gzip-compressed FASTQ files in ``fastq/``.  Use genuine biological replicates;
repeated sequencing of one library is not an independent biological replicate.

Define the matched design
-------------------------

Create ``config/samples.tsv`` with the exact tab-separated headers and paths.
This example uses single-end reads:

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

Edit the copied configuration.  Preserve the rest of the template, set the
correct input and adapter values, request the full preset, and state the
contrast direction explicitly:

.. code-block:: yaml

   biologySettings:
     genome: "data/genome.fa"
     annotation: "data/annotation.gff3"
     samples: "config/samples.tsv"

   differentialExpressionSettings:
     contrasts: ["Treated-Control"]
     padjCutoff: 0.05
     log2fcCutoff: 1.0

   predictionSettings:
     deepribo: "on"
     deepriboASiteOffset: 12

   workflowSettings:
     stages: "full"

For ``Treated-Control``, positive log2 fold changes mean higher signal in the
treated condition.  The shown DeepRibo offset is the template default from its
published *E. coli* setup; verify it for the organism and protocol.  It is not
the P-site offset produced by the TIS advisor.  When DeepRibo is enabled, the
FASTA sequence must use uppercase ``A``, ``C``, ``G``, ``T``, and ``N`` only.

Review all adapter, read-length, metagene-filter, and plotting choices in
:doc:`../configuration` before a production run.

Plan the run
------------

The full preset includes memory- and compute-intensive predictors and three
differential engines.  A REPARATION job can request 12 threads and 30 GB of
memory and disk; DeepRibo can request up to 20 GB of memory.  The first run may
also download Conda environments, containers, a DeepRibo model, and Swiss-Prot
data.  Confirm local or cluster quotas before starting.

HRIBO limits concurrent REPARATION instances through the supplied launchers.
For SLURM, configure the bundled profile as described in
:doc:`../getting-started`.

Validate and run
----------------

Perform the dry-run from the analysis directory:

.. code-block:: console

   $ micromamba activate hribo
   $ cd "$project_dir"
   $ "$hribo_checkout/run_hribo.sh" \
       --configfile "$project_dir/config/config.yaml" \
       --cores 20 \
       --dry-run all

Resolve every preflight error, then run:

.. code-block:: console

   $ "$hribo_checkout/run_hribo.sh" \
       --configfile "$project_dir/config/config.yaml" \
       --cores 20 \
       --rerun-incomplete \
       --show-failed-logs all

Review the analysis
-------------------

After the workflow completes, use the following order:

1. Open ``qc/multi/multiqc_report.html`` and confirm that every library has
   acceptable read quality, mapping, and depletion metrics.
2. Check ``figures/heatmap_SpearmanCorr_readCounts.pdf`` and
   ``pca/PCA_3D.html`` for unexpected sample grouping.
3. Inspect ``metageneprofiling/read_length_fractions.html`` and the per-library
   start/stop profiles, followed by each
   ``tis_advice/<library>/tis_recommendation.html`` report.
4. Load the BAMs, normalized BigWig tracks, and
   ``tracks/updated_annotation.gff`` in a genome browser.
5. Review ``auxiliary/predictions_reparation.xlsx`` and
   ``auxiliary/predictions_deepribo.xlsx`` separately before interpreting their
   combined evidence in ``auxiliary/overview.xlsx``.
6. For each contrast, review the deltaTE and xTail workbooks.  Check both
   effect size and adjusted p-value, and keep the contrast direction in view.
   If RiboRex was enabled, use its workbook as supplementary method evidence,
   not as an additional vote in the deltaTE calls.

The complete result map and interpretation guidance are in
:doc:`../outputs`; exact workbook columns and sheets are in
:doc:`../table-reference`.

Finish by repeating the dry-run.  An unchanged completed analysis should be a
no-op.  A standalone Snakemake workflow report can then be created with the
same Snakefile, directory, and configuration options plus
``--report "$project_dir/report.html"``.
