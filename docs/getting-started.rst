Getting started
===============

Supported execution environment
-------------------------------

The reproducible, tested path is Linux on x86-64.  HRIBO needs:

* `Conda or Micromamba
  <https://mamba.readthedocs.io/en/stable/installation/micromamba-installation.html>`_
  to install the launcher and per-rule environments;
* Apptainer to execute the digest-pinned DeepRibo, Reparation, and deltaTE
  containers; and
* enough local or cluster scratch space for alignments, environments, container
  images, and predictor work directories.

CI currently exercises Snakemake 9.25.2 and Apptainer 1.5.3.  The checked-in
``environment.linux-64.pin.txt`` is the exact launcher environment used for
release validation.  On first use, each selected stage downloads any uncached
rule environments it needs.  Container-backed prediction and differential
stages also download their images; the prediction stage retrieves the DeepRibo
model when enabled and a checksum-verified UniProtKB/Swiss-Prot release.  Plan
network access and disk space for the stages being run.

Apptainer is a host runtime; it is deliberately not installed by HRIBO's
Conda environment.  Use a cluster-provided module or follow the official
`Apptainer installation guide
<https://apptainer.org/docs/admin/main/installation.html>`_, then verify
``apptainer version`` before starting.  On a managed cluster, ask the
administrator which module and container-cache location are supported.

Install HRIBO
-------------

Clone the workflow over HTTPS and create its exact launcher environment.  In
the examples below, ``/path/to/HRIBO`` is the checkout and the analysis lives
in a separate directory.  While 2.0 remains under development, select the
``development`` branch explicitly so Git does not clone the older default
branch; use a released 2.0 tag when one becomes available:

.. code-block:: console

   $ git clone --branch development --single-branch \
       https://github.com/RickGelhausen/HRIBO.git /path/to/HRIBO
   $ micromamba create --name hribo --file /path/to/HRIBO/environment.linux-64.pin.txt
   $ micromamba activate hribo

``conda create --name hribo --file
/path/to/HRIBO/environment.linux-64.pin.txt`` is equivalent when Conda is
preferred.  The readable ``environment.yaml`` documents the direct
requirements; the explicit ``.pin.txt`` file is the reproducible install input.

Create an analysis directory
----------------------------

Keep data and results outside the Git checkout.  From a new project directory,
copy the configuration templates and then edit the copies:

.. code-block:: console

   $ mkdir -p /path/to/my-analysis/config
   $ cp /path/to/HRIBO/config/config.yaml /path/to/my-analysis/config/
   $ cp /path/to/HRIBO/config/samples.tsv /path/to/my-analysis/config/
   $ cd /path/to/my-analysis

Supply a matching reference genome in FASTA, an annotation in GFF3 or GTF, and
gzip-compressed FASTQ files.  Paths in ``config/config.yaml`` and
``config/samples.tsv`` are resolved from the analysis directory.  Absolute
paths are also accepted.

Check the setup before computing
--------------------------------

Use an explicit Snakefile, working directory, and config file.  A dry-run
performs HRIBO's schema and semantic preflight before constructing the DAG:

.. code-block:: console

   $ snakemake \
       --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory "$PWD" \
       --configfile "$PWD/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 20 \
       --resources reparation_instances=1 \
       --dry-run all

Preflight reports configuration, sample-design, FASTQ, reference, coordinate,
and sequence-identifier problems together, before expensive jobs begin.

Run the workflow
----------------

Remove ``--dry-run`` and add failure-friendly execution options:

.. code-block:: console

   $ snakemake \
       --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory "$PWD" \
       --configfile "$PWD/config/config.yaml" \
       --software-deployment-method conda apptainer \
       --cores 20 \
       --resources reparation_instances=1 \
       --rerun-incomplete \
       --latency-wait 60 \
       --printshellcmds \
       --show-failed-logs all

The repository's ``run_hribo.sh`` is a local convenience wrapper for this
pattern.  Invoke it from the analysis directory; it locates the workflow from
its own checkout and forwards additional arguments to Snakemake:

.. code-block:: console

   $ /path/to/HRIBO/run_hribo.sh --configfile "$PWD/config/config.yaml" --cores 20

After a successful run, confirm that an unchanged rerun is a no-op:

.. code-block:: console

   $ snakemake ... --dry-run all

Generate the Snakemake report separately when required:

.. code-block:: console

   $ snakemake \
       --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory "$PWD" \
       --configfile "$PWD/config/config.yaml" \
       --report report.html

Continue with :doc:`samples`, :doc:`configuration`, and :doc:`stages` before a
production run.
