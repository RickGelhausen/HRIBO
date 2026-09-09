Troubleshooting
===============

Start with the first failing Snakemake rule, its file under ``logs/``, and the
preflight summary printed before DAG construction.  ``--show-failed-logs`` and
``--printshellcmds`` make that evidence visible without guessing which tool
failed:

.. code-block:: console

   $ snakemake \
       --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory /path/to/analysis \
       --configfile /path/to/analysis/config/config.yaml \
       --software-deployment-method conda apptainer \
       --cores 20 \
       --resources reparation_instances=1 \
       --rerun-incomplete \
       --printshellcmds \
       --show-failed-logs all

Paths and spaces
----------------

Relative input paths are resolved from the analysis directory, not from the
HRIBO checkout.  Use an absolute ``--configfile`` path when invoking Snakemake
directly, and run ``run_hribo.sh`` or ``slurm_run.sh`` while the analysis
directory is current.

Paths with spaces are supported at the workflow boundary.  Quote shell
variables and quote YAML strings containing spaces:

.. code-block:: yaml

   biologySettings:
     genome: "references/strain A/genome.fa"
     annotation: "references/strain A/annotation.gff3"

Do not add literal quote characters around paths inside ``samples.tsv``.  TSV
fields are not shell fragments; the workflow quotes them when constructing
commands.  If a path is reported missing, compare it with
``realpath -- /path/to/analysis/the/configured/path`` and confirm that
``--directory`` points at the intended project.

Genome and annotation identifiers
---------------------------------

The sequence identifier is the first whitespace-delimited token after ``>`` in
each FASTA header.  Every annotation seqid in column 1 must match one of those
tokens exactly, including case and version suffixes.  For example,
``NC_000913`` and ``NC_000913.3`` are different identifiers.

The preflight also rejects duplicate FASTA identifiers, malformed nine-column
annotation rows, out-of-range coordinates, embedded ``##FASTA`` sections, and
incompatible feature identifiers.  Fix the source files instead of suppressing
these checks; otherwise later mapping, counting, or ORF annotation would lose
records silently.

Requested Ribo stages disappear
-------------------------------

HRIBO adapts the requested stage set to the sample design:

* Without an ordinary ``RIBO`` library, ``predictions``,
  ``differential_expression``, and ``overview`` are removed.
* Without any ``RIBO``, ``TIS``, or ``TTS`` library, ``metagene`` and
  ``tis_advisor`` are also removed.

The resolved stages are printed at startup.  Check method spelling and case in
``samples.tsv`` if a stage was removed unexpectedly.  Differential expression
also needs matched RIBO/RNA conditions and at least two replicates on each
selected side; see :doc:`samples`.

An empty metagene profile is not necessarily a failed stage.  The current
implementation emits valid zero evidence when no read survives the selected
read lengths, mapping method, annotation, overlap, length, or RPKM filters.
Inspect the BAM, selected lengths, and filter thresholds before widening them.
See :doc:`metagene-profiling` for the coordinate and evidence conventions.

Zero mapped reads
-----------------

RPKM, CPM, ``mil``, and ``min`` are undefined for a library with zero mapped
reads, so HRIBO stops instead of emitting infinities or plausible-looking
zeros.  Work backwards through the per-library logs and FastQC/MultiQC report:

* confirm that the FASTQ layout and adapters are correct;
* confirm that the genome is the assembly represented by the reads;
* inspect how many reads remain after rRNA/tRNA depletion and unique mapping;
* check that reference identifiers are consistent; and
* verify that a copied or truncated BAM was not substituted for the workflow
  output.

Do not edit the mapped-read summary or insert a fake denominator.  Correct the
input or configuration and let Snakemake rebuild the affected outputs.

Conda and Apptainer failures
----------------------------

Activate the pinned launcher environment and verify the two deployment tools
before retrying:

.. code-block:: console

   $ snakemake --version
   $ conda --version
   $ apptainer version

Use ``--software-deployment-method conda apptainer`` (or ``--sdm conda
apptainer``).  On first use, a selected stage needs network access for any
uncached rule environments and assets it requires.  Container-backed prediction
and differential stages pull digest-pinned images; prediction additionally
retrieves Swiss-Prot and, when enabled, the DeepRibo model.  A cluster compute
node must either have that access or use caches populated from an allowed node.

When temporary or cache space is small, point Apptainer at a private filesystem
with enough capacity before launching Snakemake.  The directories must already
be writable and executable by the job:

.. code-block:: console

   $ export APPTAINER_CACHEDIR=/scratch/my-user/apptainer-cache
   $ export APPTAINER_TMPDIR=/scratch/my-user/apptainer-tmp
   $ mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

Do not replace digest-pinned container references or install a predictor by
hand to bypass a pull failure.  That changes the scientific environment.  Save
the failing rule, image reference, tool versions, and log when reporting the
problem.

Recovery after interruption
---------------------------

If a machine stopped during a run, first make sure no HRIBO process is still
active.  If Snakemake reports a stale working-directory lock, unlock that exact
analysis directory once, then resume with ``--rerun-incomplete``:

.. code-block:: console

   $ snakemake --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory /path/to/analysis --unlock
   $ /path/to/HRIBO/run_hribo.sh --rerun-incomplete all

The second command must be run from ``/path/to/analysis``.  Do not unlock a
directory while another workflow or cluster job is using it.

Predictor completion receipts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REPARATION owns a validated result tree under
``reparation/<condition>-<replicate>/``.  DeepRibo cutoff estimation owns
``parameters.txt`` and ``s_curve.png`` under
``deepribo/cutoffs/<condition>-<replicate>/``.  Each owner publishes
``.complete`` only after its sibling artifacts have been validated and durably
committed.

If a runner-owned sibling is later missing or corrupt while ``.complete`` still
exists, a downstream rule fails deliberately.  Snakemake tracks the receipt,
so it will not infer that the owner must run again from an undeclared sibling.
After confirming the damaged file and saving any evidence needed for diagnosis,
delete **only** the matching ``.complete`` receipt.  Then request the normal
downstream target with ``--rerun-incomplete``.  Deleting that receipt is an
intentional request for the checked runner to repair or replace its complete
owned publication.

For REPARATION the receipt is, for example,
``reparation/A-1/.complete``; a suitable downstream target is
``reparation/A-1.reparation.gff``.  For cutoff estimation the receipt is
``deepribo/cutoffs/A-1/.complete``; requesting
``deepribo/A-1/predictions.csv`` rebuilds the cutoff if required.

Do not delete arbitrary predictor files, transaction backups, staging
directories, or the complete result directory to force a retry.  Those files
are managed as an atomic group, and manual partial cleanup can discard a valid
rollback snapshot.  If the checked runner reports an ambiguous or unsafe
transaction state, preserve the directory and report the exact error.

SLURM execution
---------------

The bundled profile uses the Snakemake SLURM executor plugin.  It is included
in the pinned launcher environment.  Before the first submission, copy or edit
``workflow/profiles/slurm/config.yaml`` for the site's account and partition;
leave per-rule memory, runtime, and thread requirements with the rules unless a
site-specific override is required.

Activate any site modules and the HRIBO launcher environment before invoking
``slurm_run.sh``.  The launcher intentionally contains no institution-specific
module commands or Conda activation.  Test DAG construction before submitting:

.. code-block:: console

   $ cd /path/to/analysis
   $ /path/to/HRIBO/slurm_run.sh --dry-run all
   $ /path/to/HRIBO/slurm_run.sh --rerun-incomplete --show-failed-logs all

If jobs remain pending, inspect the scheduler reason and compare the requested
account, partition, memory, runtime, and node constraints with local policy.
Increasing Snakemake's ``jobs`` value cannot resolve an invalid account or an
unavailable resource request.
