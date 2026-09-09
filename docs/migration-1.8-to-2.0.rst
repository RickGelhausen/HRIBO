Migrating from HRIBO 1.8 to 2.0
================================

HRIBO 2.0 is a workflow migration as well as a dependency update.  Keep the
1.8 result directory unchanged and run 2.0 in a separate analysis directory.
This preserves a usable baseline and prevents old intermediate files from
satisfying a rule whose implementation or output contract has changed.

Before starting, record the exact 1.8 tag or commit, archive the old
configuration and sample sheet, and checksum the genome, annotation, and FASTQ
inputs.  Compare like with like: use the same biological inputs, adapter
settings, conditions, replicates, read-length selections, and explicit
contrasts in both runs.

Distribution and repository layout
----------------------------------

HRIBO is now distributed only as a Snakemake workflow.  It is not an
installable Python package, and commands such as ``pip install .`` are not a
supported installation path.  Clone a release or checkout, create the pinned
launcher environment, and let Snakemake create the isolated environments and
containers required by individual rules.

The repository follows the standard Snakemake workflow layout:

.. list-table:: Path changes
   :header-rows: 1
   :widths: 35 35 30

   * - HRIBO 1.8
     - HRIBO 2.0
     - Purpose
   * - ``Snakefile``
     - ``workflow/Snakefile``
     - Workflow entry point
   * - ``rules/``, ``scripts/``, ``envs/``
     - ``workflow/rules/``, ``workflow/scripts/``, ``workflow/envs/``
     - Workflow implementation
   * - ``schemas/``
     - ``workflow/schemas/``
     - Input schemas
   * - ``templates/config.yaml`` and ``templates/samples.tsv``
     - ``config/config.yaml`` and ``config/samples.tsv``
     - User templates
   * - ``slurm_profile/``
     - ``workflow/profiles/slurm/``
     - SLURM executor profile
   * - Separate ``HRIBO_ReadTheDocs`` repository
     - ``docs/`` in this repository
     - Versioned documentation

The maintained pages are rewritten against the current schemas and tested
outputs; the obsolete repository is not a second source of truth.  Its complete
history and namespaced release tags are preserved on
``archive/hribo-readthedocs`` as described in
:ref:`documentation-source-and-legacy-archive`.

Keep the checkout and analysis separate.  All paths in the configuration and
sample sheet are resolved from the analysis directory selected with
``--directory``.  The local and SLURM launchers resolve the workflow from their
own checkout, so they can be invoked from any analysis directory.

Launcher and deployment environment
-----------------------------------

Create the launcher from ``environment.linux-64.pin.txt`` rather than
installing an unpinned Snakemake version.  The readable ``environment.yaml``
describes its direct dependencies, while the explicit pin is the validated
Linux installation input.  Activate that environment before calling
``run_hribo.sh``, ``slurm_run.sh``, or Snakemake directly.

Use both supported deployment methods:

.. code-block:: console

   $ snakemake \
       --snakefile /path/to/HRIBO/workflow/Snakefile \
       --directory /path/to/new-analysis \
       --configfile /path/to/new-analysis/config/config.yaml \
       --software-deployment-method conda apptainer \
       --cores 20 \
       --resources reparation_instances=1 \
       --dry-run all

The old ``--use-conda``/``--use-singularity`` examples and legacy cluster
submission flags no longer describe the supported execution path.  Snakemake
9 uses an executor plugin and the bundled profile for SLURM.

Configuration and sample design
-------------------------------

Start with a fresh copy of the 2.0 templates and transfer values deliberately.
Do not copy the complete 1.8 configuration over the new template: 2.0 validates
the document against ``workflow/schemas/config.schema.yaml`` and rejects stale
or unknown nested keys.

The main changes are:

* ``workflowSettings.workflow`` is replaced by
  ``workflowSettings.stages``.  It accepts ``full``, ``preprocessing``, or a
  YAML list of named stages.  A command-line ``--config stages=...`` override
  takes precedence for one run.
* The old ``differentialExpression`` on/off key is removed.  Differential
  analysis runs when ``differential_expression`` is selected.  The ``full``
  preset includes it and therefore requires a matched RIBO/RNA design.
* ``differentialExpressionSettings.contrasts`` is a YAML list, for example
  ``["Treated-Control"]``.  Positive fold changes mean higher signal in the
  condition on the left.
* ``predictionSettings.deepriboASiteOffset`` makes the DeepRibo 3'-end-to-A-site
  offset explicit.  It is not the 5'-end-to-P-site value reported by the TIS
  advisor.
* TIS-advisor settings and explicit xTail inclusion controls are now present.
  Preserve the 2.0 defaults unless the experiment provides a reason to change
  them.

The sample-sheet header remains exactly:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2

Parsing and design validation are stricter.  Method names are canonical and
case-sensitive; condition names cannot contain the hyphen used to delimit a
contrast; replicates are canonical positive integers without signs, whitespace,
or leading zeros; and each
``method``/``condition``/``replicate`` tuple must be unique.  See
:doc:`samples` for the accepted methods and stage-specific replicate rules.
Both FASTQ headers are required even for single-end data, whose ``fastqFile2``
cells stay empty.  Referenced inputs must contain gzip-compressed bytes.  The
preflight now scans the complete gzip stream and every FASTQ record, so a late
CRC failure, truncated final record, or malformed record stops before trimming.

Stages and output contracts
---------------------------

Stage selection now controls the requested deliverables.  Dependencies are
still built automatically, so selecting ``predictions`` also runs the required
read processing, mapping, and counting.  Conversely, deselecting a stage is the
supported way to stop at an earlier result such as the final BAMs.  Review
:doc:`stages` instead of translating the 1.8 ``full`` switch mechanically.

Several output contracts changed or became explicit:

* Workbooks use consistent headers across generators.  Differential tables use
  canonical names such as ``log2FC``, ``log2FC_SE``, and
  ``pvalue_adjusted`` (with RIBO/RNA/TE prefixes where appropriate).
* Browser-facing files in ``tracks/`` are validated GFF3.  The historical
  ``readcounts/*_annotation.gff`` files are internal count matrices with extra
  sample columns, not browser GFF3 files.
* ``auxiliary/overview.xlsx`` and its ``overview.tsv``, ``overview.gff``, and
  ``overview_misc.gff`` companions are declared workflow outputs rather than
  untracked side effects.
* ``maplink/*.bam`` remains the convenient final-alignment location, but its
  entries are relative links into ``bam/``.  Move the complete result tree or
  copy them with link dereferencing.
* The old report-archive helper is gone.  Request a native Snakemake report
  separately with ``snakemake ... --report report.html`` after the analysis.

Consult :doc:`outputs` before updating a downstream parser.  Prefer column
names and feature identifiers over fixed column positions.

Scientific changes that affect comparisons
-------------------------------------------

Library-wide normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~

HRIBO 1.8 normalized each reference sequence independently.  On a genome with
a chromosome and one or more plasmids, this made the denominator depend on the
feature's contig and inflated values on low-count contigs.  HRIBO 2.0 uses the
sum across every contig in a library:

.. math::

   \mathrm{RPKM} =
   \frac{10^9 \times \mathrm{feature\ count}}
        {\mathrm{feature\ length\ (nt)} \times
         \mathrm{library\ mapped\ reads}}

The same library-wide denominator is used for ``mil`` BigWigs and CPM
metagenes.  ``min`` tracks scale every library to the smallest complete-library
mapped-read total.  Multi-mapping alignments contribute ``1/NH`` to the
effective mapped-read total, matching fractional feature counting.

Single-contig datasets are normally unchanged by the contig-to-library part of
this correction.  Multi-contig workbooks, coverage tracks, metagene CPM values,
translation-efficiency summaries, and RPKM-based filters can change
substantially.  Re-evaluate thresholds rather than applying a factor copied
from a 1.8 result.

Metagene orientation and completeness
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Metagene windows are now consistently oriented in transcript 5'-to-3'
direction.  Start profiles use ``[-outside, inside)`` coordinates.  Stop
profiles use ``[-inside, outside)`` coordinates, so negative positions are in
the coding sequence and positive positions are downstream on both strands.
The 1.8 stop profiles did not apply that orientation consistently.

Global read coverage now fills every strand/anchor combination.  The old
implementation omitted minus-strand start windows and plus-strand stop windows,
which appeared as empty halves even when reads were present.  A current all-zero
profile is instead valid evidence that no selected read survived the configured
window, read-length, mapping, and annotation filters.  Global metagene and
browser-track coverage now also follow the read's CIGAR-aligned blocks instead
of painting deletions, reference skips, or soft clips as observed bases.
Centered browser-track coverage distributes its weight only over aligned
positions remaining after clipping; a centered metagene remains a single
midpoint assignment.

Configured-but-unobserved lengths are retained as explicit zero columns, while
observed lengths outside the configured selection are excluded.  P-site markers
appear only on physical 5' or 3' read-end profiles, on the correct side of the
start codon, and are estimated from raw counts before ``cpm`` or ``window``
presentation normalization.

Prediction boundaries and publication receipts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DeepRibo reverse-strand A-site geometry, predictor parsing, ORF aggregation,
and strict GFF3 publication have been corrected.  REPARATION and DeepRibo
cutoff estimation now publish validated groups of files transactionally.  A
non-empty ``.complete`` file is the Snakemake-tracked receipt for each group;
downstream rules consume siblings only after that receipt exists.  See
:doc:`troubleshooting` for the narrow repair procedure after a damaged sibling.

These changes can alter prediction coordinates or membership.  Review exact
coordinate overlap and a small biological tolerance rather than expecting
byte-identical predictor files.

Compare the old and new runs
----------------------------

The repository includes a semantic comparator for a controlled real-data
migration.  Run it from the pinned development environment, which supplies its
``openpyxl`` and ``pysam`` dependencies:

.. code-block:: console

   $ python .github/scripts/compare_real_data_outputs.py \
       /path/to/hribo-1.8-results \
       /path/to/hribo-2.0-results \
       --report /path/to/comparison.json

It discovers primary tables, GFF3 files, BAMs, BigWigs, and TIS evidence.  It
compares tabular values by stable keys, reports exact and within-three-nucleotide
prediction overlap, and distinguishes semantic equality from byte equality.
Plots are evaluated through their source tables rather than pixels.

Missing or malformed candidate outputs are errors.  Scientific value changes
are reported for review but do not by themselves make the command fail.  Add
``--minimum-key-overlap`` or ``--minimum-correlation`` only after choosing a
defensible project-specific threshold.  ``--allow-missing GLOB`` is repeatable,
but should be used only for a documented, intentional output removal.

Review the JSON report alongside MultiQC and the primary workbooks.  Expected
differences include the normalization and metagene corrections above; they are
not a blanket explanation for unrelated feature loss, changed contrast
direction, or malformed output.
