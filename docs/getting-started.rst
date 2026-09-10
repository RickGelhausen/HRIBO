Install and run HRIBO
=====================

This guide takes a new analysis from installation to its first results.  The
commands use ``/path/to/HRIBO`` for the software checkout and
``/path/to/my-analysis`` for the project that contains the data and results.
Keep these directories separate so that updating HRIBO does not mix source
files with an analysis.

Requirements
------------

The supported execution environment is Linux on x86-64.  You need:

* `Conda or Micromamba
  <https://mamba.readthedocs.io/en/stable/installation/micromamba-installation.html>`_;
* `Apptainer
  <https://apptainer.org/docs/admin/main/installation.html>`_ installed on the
  host or provided as a cluster module;
* a reference genome and matching annotation; and
* gzip-compressed FASTQ files.

The first run needs network access to create rule-specific Conda environments
and download any containers or model data required by the selected stages.
Prediction stages also need substantial disk space; see :doc:`stages` before
requesting every analysis.

1. Install HRIBO
----------------

Clone the HRIBO 2.0.0 release and create the pinned launcher environment.

.. code-block:: console

   $ git clone --branch 2.0.0 --single-branch \
       https://github.com/RickGelhausen/HRIBO.git /path/to/HRIBO
   $ micromamba create --name hribo \
       --file /path/to/HRIBO/environment.linux-64.pin.txt
   $ micromamba activate hribo
   $ snakemake --version
   $ apptainer version

Use ``conda create`` instead of ``micromamba create`` if Conda is preferred.
HRIBO is run directly from its checkout; it is not installed with ``pip``.

2. Create an analysis directory
-------------------------------

Copy the supplied configuration and sample-sheet templates into a new project:

.. code-block:: console

   $ mkdir -p /path/to/my-analysis/{config,data,fastq}
   $ cp /path/to/HRIBO/config/config.yaml \
       /path/to/my-analysis/config/config.yaml
   $ cp /path/to/HRIBO/config/samples.tsv \
       /path/to/my-analysis/config/samples.tsv

The project will look like this before the first run:

.. code-block:: text

   my-analysis/
   ├── config/
   │   ├── config.yaml
   │   └── samples.tsv
   ├── data/
   │   ├── genome.fa
   │   └── annotation.gff3
   └── fastq/
       └── ribo-control-1.fastq.gz

Copy or link your inputs into these directories.  The FASTA and annotation
must describe the same assembly, and annotation sequence identifiers must
match the first identifier on the corresponding FASTA headers exactly.

3. Describe the libraries
-------------------------

Edit ``config/samples.tsv``.  For a single Ribo-seq library, the sheet can be:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2
   RIBO	Control	1	fastq/ribo-control-1.fastq.gz	""

The final field is empty because the library is single-end.  For paired-end
data it contains the read-2 FASTQ.  See :doc:`samples` for accepted methods,
paired-end examples, and the design required for differential analysis.

4. Configure the analysis
-------------------------

Edit the copied ``config/config.yaml`` rather than replacing it with the short
excerpt below.  At minimum, set the three input paths, the correct adapter
sequence or sequences, and the stages you want:

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

An empty adapter setting means that no sequence is supplied to the trimming
step.  Do not leave it empty merely because the adapter is unknown; establish
the library preparation before a production analysis.  The main settings are
explained in :doc:`configuration`; the copied template contains every field.

5. Validate with a dry-run
--------------------------

Activate the launcher environment, change to the analysis directory, and run
the bundled local launcher with ``--dry-run``:

.. code-block:: console

   $ micromamba activate hribo
   $ cd /path/to/my-analysis
   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --cores 8 \
       --dry-run all

The dry-run validates the configuration, sample design, input files, reference
coordinates, and sequence identifiers before scheduling work.  Fix every
reported error before continuing.  It is normal for the dry-run to list all
jobs that would be created.

6. Run the workflow
-------------------

Remove ``--dry-run`` and enable safe restart of incomplete jobs:

.. code-block:: console

   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --cores 8 \
       --rerun-incomplete \
       --show-failed-logs all

``run_hribo.sh`` selects HRIBO's Snakefile, enables Conda and Apptainer,
prints the commands being executed, keeps independent jobs running after a
failure, and limits memory-intensive REPARATION instances.  Arguments supplied
after the script name are forwarded to Snakemake.

A successful run ends with ``Done, no error``.  Repeat the dry-run command; an
unchanged completed analysis should report that there is nothing to be done.

7. Open the first results
-------------------------

For the three-stage example above, start with:

* ``qc/multi/multiqc_report.html`` for read-processing and mapping quality;
* ``maplink/RIBO-Control-1.bam`` and its ``.bai`` index for final uniquely
  mapped alignments; and
* ``globaltracks/``, ``centeredtracks/``, ``fiveprimetracks/``, and
  ``threeprimetracks/`` for strand-aware BigWig coverage.

For a larger analysis, :doc:`outputs` identifies the main result from every
stage and explains what to inspect first.

Run on SLURM
------------

Edit the account and partition in
``workflow/profiles/slurm/config.yaml`` for your cluster, then dry-run and
submit from the analysis directory:

.. code-block:: console

   $ /path/to/HRIBO/slurm_run.sh \
       --configfile "$PWD/config/config.yaml" \
       --dry-run all
   $ /path/to/HRIBO/slurm_run.sh \
       --configfile "$PWD/config/config.yaml" \
       --rerun-incomplete all

The bundled profile uses the Snakemake SLURM executor plugin.  Site modules,
partitions, accounts, and storage paths are cluster-specific; confirm them with
your administrator before submitting a full analysis.

Run Snakemake directly
----------------------

The wrapper is recommended for routine local use.  If you need a direct
Snakemake command, the equivalent core options are:

.. code-block:: console

   $ snakemake \
       --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory "$PWD" \
       --configfile "$PWD/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 8 \
       --resources reparation_instances=1 \
       --rerun-incomplete all

After the selected results finish, a standalone Snakemake report can be
created with the same Snakefile, directory, and configuration options plus
``--report report.html``.
