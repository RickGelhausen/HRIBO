Troubleshooting
===============

Start with the first failing Snakemake rule, not the last downstream message.
The terminal output names the rule and usually points to its file below
``logs/``.  Rerun with failed-log display enabled when needed:

.. code-block:: console

   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --cores 8 \
       --rerun-incomplete \
       --show-failed-logs all

Nothing runs after a preflight error
------------------------------------

This is intentional.  HRIBO validates the requested inputs before constructing
the job graph so that malformed data do not fail hours later inside another
tool.  Read the complete preflight summary; it groups related errors and can
report more than one problem at once.

Common causes include:

* a misspelled configuration key or stage name;
* a missing sample-sheet column or duplicate library identifier;
* a relative path evaluated from the wrong analysis directory;
* a corrupt, truncated, uncompressed, or malformed FASTQ;
* genome and annotation sequence identifiers that do not match; and
* a requested differential contrast that the sample design cannot support.

Fix the source configuration or input file and repeat the dry-run.  Do not
disable the validation or edit generated files to bypass it.

Paths and spaces
----------------

Relative paths in ``config.yaml`` and ``samples.tsv`` are resolved from the
analysis directory selected by ``--directory``.  Run the launchers while that
directory is current and give ``--configfile`` an absolute path when in doubt:

.. code-block:: console

   $ cd "/path/to/my analysis"
   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --dry-run all

Paths containing spaces are supported when shell variables and YAML values are
quoted.  Do not put literal quote characters around paths in ``samples.tsv``;
TSV cells are paths, not shell fragments.

Genome and annotation identifiers do not match
----------------------------------------------

For a FASTA header such as:

.. code-block:: text

   >NC_000913.3 Escherichia coli K-12

the identifier is ``NC_000913.3``.  Column 1 of every GFF3/GTF feature must use
that exact value, including capitalization and version suffix.  The annotation
must also use valid nine-column rows and coordinates within the corresponding
FASTA sequence.

Use the FASTA and annotation from the same assembly release.  Changing only a
header to make validation pass can attach features to the wrong sequence.

A requested stage is missing
----------------------------

HRIBO prints the resolved stage list at startup and removes analyses that
cannot use the available library types:

* without ``RIBO``, it removes ``predictions``,
  ``differential_expression``, and ``overview``;
* without any ``RIBO``, ``TIS``, or ``TTS`` library, it also removes
  ``metagene`` and ``tis_advisor``.

Differential expression additionally requires matched RIBO/RNA conditions and
at least two biological replicates on each selected side.  Check method names,
condition labels, and replicates against :doc:`samples`.

Zero mapped reads
-----------------

RPKM, CPM, and normalized coverage cannot be calculated for a library with no
mapped reads, so HRIBO stops rather than creating plausible-looking zeros.
Work backwards through the report and logs:

1. confirm that the raw FASTQ contains the intended library;
2. check whether trimming retained reads of the expected lengths;
3. verify that the genome is the correct strain and assembly;
4. inspect total and unique mapping rates; and
5. check whether rRNA/tRNA depletion removed nearly all remaining reads.

Correct the input or configuration, then rerun with ``--rerun-incomplete``.
Do not insert a fake mapped-read total into a generated summary.

Empty metagene or TIS results
-----------------------------

An all-zero metagene profile or ``confidence: none`` TIS recommendation can be
a valid result rather than a workflow failure.  Confirm that the rule completed
without an error, then inspect:

* the configured read-length range and mapping methods;
* mapped depth in the final BAM;
* the metagene overlap, length, and RPKM filters;
* annotation start/stop positions; and
* whether the library preparation is expected to produce a sharp initiation
  signal.

See :doc:`metagene-profiling` and :doc:`tis-advisor` before widening filters or
assigning an offset manually.

REPARATION corrects a P-site estimate
-------------------------------------

Plastid can occasionally estimate an offset that is equal to or longer than
the corresponding read length, particularly for a sparsely represented or
noisy footprint length.  Such an offset cannot identify a position inside the
read.  HRIBO therefore replaces only that estimate with the ``default`` value
from Plastid's offset table **before** REPARATION calculates occupancy and ORF
predictions.  The REPARATION log reports the read length, original estimate,
and replacement; the published ``p_site_offsets.txt`` also records the
correction in a comment.

The warning does not by itself mean that the run failed.  Inspect
``reparation/<condition>-<replicate>/p_site_offset.png`` when present and the
corresponding ``p_site_offsets.txt``.  One corrected low-depth length may be a
reasonable fallback.  Corrections across many well-represented lengths suggest
that the initiation signal, annotation, footprint-length range, or library
quality needs closer review.  HRIBO still validates the final table and stops
if no physically valid result is available.

Conda or Apptainer fails
------------------------

Activate the pinned launcher environment and verify both tools:

.. code-block:: console

   $ micromamba activate hribo
   $ snakemake --version
   $ apptainer version

The run must enable both ``conda`` and ``apptainer`` deployment methods; the
bundled launchers already do so.  The first execution of a stage may need to
download Conda packages, a pinned container, a model, or reference data.  A
compute node must either have network access or use caches populated through
your site's supported process.

If temporary storage is too small, ask the cluster administrator for suitable
locations for ``APPTAINER_CACHEDIR`` and ``APPTAINER_TMPDIR``.  Do not replace
the pinned image with an unrelated local installation: that changes the
analysis environment.

Resume an interrupted run
-------------------------

Make sure no HRIBO process or cluster job is still using the analysis
directory, then rerun the normal target with ``--rerun-incomplete``.  Snakemake
will keep completed outputs and schedule missing or incomplete work.

If Snakemake reports a stale directory lock after a confirmed interruption,
unlock that exact analysis once:

.. code-block:: console

   $ cd /path/to/my-analysis
   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --unlock

Then launch the usual command with ``--rerun-incomplete all``.  Never unlock a
directory while another local process or cluster job is active.

SLURM jobs do not start
-----------------------

First confirm that the dry-run succeeds locally.  Then check the account and
partition in ``workflow/profiles/slurm/config.yaml`` and submit a small job.
If it remains pending, inspect the scheduler's reason and compare the requested
memory, runtime, nodes, account, and partition with site policy.  Increasing
Snakemake's job limit cannot fix an invalid account or unavailable resource.

Cluster module loading and environment activation belong in the environment
from which ``slurm_run.sh`` is launched.  The profile does not guess
site-specific module names.

Report a reproducible problem
-----------------------------

When opening an issue, include:

* the HRIBO version or commit;
* the command and effective stage selection;
* the configuration and a redacted sample sheet;
* the first failing rule and complete relevant log;
* Snakemake and Apptainer versions; and
* whether the failure occurs locally or under a scheduler.

Do not attach confidential sequencing data.  A small non-sensitive reproducer
is more useful than an entire result directory.
