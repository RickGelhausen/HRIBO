Recipe: mapping, QC, and tracks
===============================

This recipe runs one single-end Ribo-seq library through trimming, mapping,
quality control, and coverage-track generation.  It is a useful first run when
you want to check the installation and inspect one library without running ORF
prediction or differential analysis.

Before starting, install HRIBO as described in :doc:`../getting-started` and
replace the three input paths below with files from your experiment.

Create the project
------------------

.. code-block:: console

   $ hribo_checkout=/path/to/HRIBO
   $ project_dir=/path/to/hribo-mapping-example
   $ mkdir -p "$project_dir"/{config,data,fastq}
   $ cp "$hribo_checkout/config/config.yaml" "$project_dir/config/"
   $ cp "$hribo_checkout/config/samples.tsv" "$project_dir/config/"
   $ cp /path/to/genome.fa "$project_dir/data/genome.fa"
   $ cp /path/to/annotation.gff3 "$project_dir/data/annotation.gff3"
   $ cp /path/to/ribo-control-1.fastq.gz \
       "$project_dir/fastq/ribo-control-1.fastq.gz"

The FASTA and annotation must come from the same assembly, and the FASTQ must
be gzip-compressed.

Describe the library
--------------------

Replace ``config/samples.tsv`` with the following tab-separated content.  The
empty final field marks a single-end library:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2
   RIBO	Control	1	fastq/ribo-control-1.fastq.gz	""

Edit the copied ``config/config.yaml``.  Keep the rest of the template and
change these values:

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

Set the real 3' and/or 5' adapter sequence used by the experiment.  Leave a
field empty only when no adapter should be supplied to trimming.

Validate and run
----------------

Run from the analysis directory.  Start with a dry-run:

.. code-block:: console

   $ micromamba activate hribo
   $ cd "$project_dir"
   $ "$hribo_checkout/run_hribo.sh" \
       --configfile "$project_dir/config/config.yaml" \
       --cores 8 \
       --dry-run all

Fix every preflight error.  Then execute the same target:

.. code-block:: console

   $ "$hribo_checkout/run_hribo.sh" \
       --configfile "$project_dir/config/config.yaml" \
       --cores 8 \
       --rerun-incomplete \
       --show-failed-logs all

The first execution downloads and creates the environments needed by these
stages.  Runtime depends mainly on read count, reference size, available cores,
storage speed, and whether dependencies are already cached.

Review the results
------------------

After ``Done, no error`` appears, inspect:

``qc/multi/multiqc_report.html``
   Compare raw and trimmed read quality, mapping rates, unique alignments, and
   the effect of rRNA/tRNA filtering.

``maplink/RIBO-Control-1.bam``
   Load the final alignment and its ``.bam.bai`` index in a genome browser.

``globaltracks/``, ``centeredtracks/``, ``fiveprimetracks/``, and ``threeprimetracks/``
   Load the BigWig tracks to compare whole-read, central, 5'-end, and 3'-end
   signal.  Reverse-strand values are negative only for browser display.

See :doc:`../outputs` for the filename pattern and normalization choices.
Finally, repeat the dry-run.  An unchanged successful project should report
that there is nothing to be done.
