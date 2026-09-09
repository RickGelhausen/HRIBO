Real-data release validation
============================

Representative biological data have **not yet been run** through the current
development version.  The automated tests exercise fixtures, simulated signal,
workflow DAGs, and production container boundaries, but they do not replace a
side-by-side biological review.  Completing and recording the protocol below
is therefore a release gate.

Prepare comparable runs
-----------------------

Use a dataset that represents the organisms, replicons, read lengths, library
types, replicate structure, and low-coverage cases used in production.  A
multi-contig or chromosome-plus-plasmid assembly is important because several
normalization corrections are invisible on a single contig.

The public datasets used by the former documentation are recorded in
:doc:`historical-example-data`.  They are candidates for this exercise, not
pre-approved fixtures: confirm their current metadata and input preparation,
choose the biologically appropriate subset, and record all decisions with the
run.

#. Put the baseline release and current candidate in separate, clean analysis
   directories.  Never run one version over the other's outputs.
#. Use identical FASTQ files, reference sequences, annotation, sample design,
   adapter settings, read-length ranges, contrasts, and feature filters wherever
   the two configuration formats permit it.  Save both effective
   configurations and record every unavoidable migration difference.
#. Run the biologically relevant complete stage set for both versions.  Keep
   logs, environment/lock information, container digests, commands, and Git
   revisions beside the results.
#. Confirm both runs completed and that an unchanged candidate dry-run is a
   no-op before comparing outputs.

Run the comparator
------------------

From the root of the current HRIBO checkout, run:

.. code-block:: console

   $ python .github/scripts/compare_real_data_outputs.py \
       /path/to/baseline-results \
       /path/to/candidate-results \
       --report /path/to/validation/real-data-comparison.json

The comparator discovers primary workbooks and delimited tables, prediction and
overview GFFs, final BAMs, TIS JSON/TSV evidence, the correlation matrix, and
BigWigs below the four primary track trees.  It compares tables by normalized
headers and stable row keys, reports numeric deltas and correlations, reports
exact and within-3-nt GFF overlap, and uses an order-independent full-alignment
digest for BAM records.

BigWigs are checked without an optional interval-decoding dependency: the tool
validates the BBI header, zoom metadata, chromosome B+ tree, and global numeric
summary, but its final changed/unchanged decision is conservatively based on a
compressed-byte digest.  A changed BigWig therefore always needs the manual
track review below.  HTML and PDF rendering is also a manual check; plots are
not compared pixel by pixel.

Comparator options
------------------

``--report PATH``
   Required destination for the deterministic JSON report.  Parent directories
   are created when needed.

``--allow-missing GLOB``
   Allow a baseline artifact that was intentionally removed from the candidate.
   The option is repeatable and matches paths relative to the result root.
   Quote shell globs.  It does not suppress a malformed artifact and does not
   hide the removal from the report.

``--minimum-key-overlap VALUE``
   Add an error when a comparable table's row-key Jaccard similarity is below
   this value.  The value must be between 0 and 1.

``--minimum-correlation VALUE``
   Add an error when a numeric column with at least three finite paired values
   has Spearman correlation below this value.  A changed column with at least
   three finite pairs also fails when its correlation is undefined; an unchanged
   constant column is exempt.  The value must be between 0 and 1.

For example, after choosing and recording project-specific acceptance limits:

.. code-block:: console

   $ python .github/scripts/compare_real_data_outputs.py \
       baseline candidate \
       --report validation/comparison.json \
       --allow-missing 'auxiliary/retired_legacy_output.xlsx' \
       --minimum-key-overlap 0.90 \
       --minimum-correlation 0.95

There are deliberately no default scientific thresholds.  Choose them before
examining the candidate, justify them for the dataset, and retain them with the
report.

Narrow HRIBO 1.8 compatibility
------------------------------

The baseline parser contains two narrowly scoped compatibility paths for
artifacts that HRIBO 1.8 could emit.  Each use is recorded in artifact metadata,
adds a warning, and requires review:

* In known top-level ``tracks/`` DeepRibo outputs, a baseline CDS record whose
  source is ``deepribo``, whose phase is an ASCII integer, and which has no
  ``deepribo_distance`` attribute is interpreted as the historical use of the
  phase field for distance/rank metadata.  Its comparison phase is normalized
  to ``0`` while the original values and counts remain in the report.
* A zero-byte baseline file matching
  ``tracks/<condition>.deepribo.gff``,
  ``tracks/<condition>.reparation.gff``, or
  ``tracks/<condition>.merged.gff`` is interpreted as an empty feature set.

These allowances never apply to the candidate, arbitrary or nested GFF paths,
or other malformed content.  Candidate empty results must still be valid GFF3:
a header-only file declaring ``##gff-version 3`` is accepted, while a zero-byte,
comment-only, or wrong-version file is invalid.

Exit status
-----------

``0``
   The comparison is structurally valid and no configured threshold failed.
   This can still include changed or candidate-only artifacts and an allowed
   removal.  Read ``summary.review_required`` and review every non-unchanged
   entry; exit 0 is not a biological approval.

``1``
   The JSON report was written but contains structural or threshold errors,
   such as an unallowlisted missing candidate artifact, an invalid output, too
   little key overlap, or too little numeric correlation.

``2``
   Setup, argument, top-level parsing, or report-writing failed.  The requested
   report may not exist and the comparison must be rerun after correcting the
   error.

Expected intentional differences
--------------------------------

Do not waive these differences wholesale.  Confirm that each change has the
expected direction and is confined to affected artifacts.

Metagene geometry
   Current stop profiles run from coding sequence toward downstream sequence on
   both strands.  Legacy stop workbooks can therefore appear mirrored.  The
   old global mapper also lost the minus/start and plus/stop slices; those
   combinations should now contain evidence when matching reads exist.  See
   :doc:`metagene-profiling` for the exact axes.  Configured-but-unobserved
   lengths are now explicit zero columns and unrequested lengths are excluded.

Library-wide normalization
   RPKM, BigWig ``mil``/``min``, and metagene ``cpm`` now use one total across
   every contig in a library.  On multi-contig data, normalized values and RPKM
   filtering can change even when raw alignments do not.  Overview RPKMs and
   TIS evidence can inherit those changes.  A single-contig library should not
   change for this reason alone.

Fractional multi-mapper totals
   Mapped-read summaries now add ``1 / NH`` for each alignment.  Total-mapping
   denominators and weighted mean read lengths can therefore change when reads
   have several reported hits; unique final BAM totals should not.  This makes
   the denominator of ``annotation_total.xlsx`` consistent with the
   fractionally counted multi-mapper numerator.

CIGAR-aware coverage and profile markers
   Current global metagenes and global/centered browser tracks exclude
   deletions, reference skips, and soft-clipped query bases from coverage.
   Results can therefore change for gapped or clipped alignments.  A centered
   metagene remains a single midpoint assignment and is not a per-base coverage
   measure.  Metagene P-site markers are now restricted to physical 5'/3' end
   profiles, placed on the correct side of the start codon, and estimated from
   raw counts independently of the chosen presentation normalization.

Manual scientific review checklist
----------------------------------

* Confirm the candidate produced every expected stage output, and explain every
  ``baseline_only``, ``candidate_only``, or ``invalid`` comparator entry.
* Compare MultiQC summaries, raw/processed read counts, rRNA/tRNA depletion,
  mapping rates, read-length distributions, and per-contig depth.  Large
  upstream changes invalidate downstream numerical comparisons until explained.
* Check BAM reference names, lengths, mapped counts, strand balance, and a
  representative set of alignments in a genome browser.
* Load raw, ``mil``, and ``min`` BigWigs for both strands and all four mapping
  styles.  Inspect chromosome and plasmid loci, confirm reverse values are
  negative by convention, and independently spot-check normalization factors.
* Inspect start and stop metagenes for plus- and minus-strand genes.  Confirm
  transcript orientation, restored global slices, expected read-length peaks,
  and that any zero/``no_evidence`` profile is explained by depth or filtering.
* Recalculate a small set of workbook RPKMs from raw feature counts, feature
  length, and the complete-library effective mapped total.  Cross-check total
  and unique annotation workbooks and translational-efficiency pairs.
* Review the TIS advisor's chosen end, per-length offsets, confidence, warnings,
  and no-recommendation cases against the metagene signal.  Do not transfer its
  P-site offsets to DeepRibo's 3'-to-A-site setting.
* For Reparation and DeepRibo, review exact and within-3-nt overlap, strand,
  start/stop geometry, predicted-only calls, lost calls, score/rank changes, and
  several browser examples rather than relying only on a global overlap rate.
* For every differential-expression contrast, verify the documented
  left-minus-right direction, log2-fold-change signs, adjusted-p-value distributions,
  row-set overlap, top hits, and agreement or justified disagreement among
  xTail, Riborex, and deltaTE.
* Cross-check ``overview.xlsx``, ``overview.tsv``, ``overview.gff``, and
  ``overview_misc.gff`` for row counts, coordinates, feature routing, prediction
  evidence, enabled contrasts, and representative numeric values.
* Open every interactive HTML and representative PDF/static figure, checking
  labels, axes, legends, missing panels, and offline Plotly behavior.
* Record each accepted difference with its cause and reviewer.  Archive the
  inputs or checksums, configurations, revisions, comparator command and JSON,
  manual notes, and final approval together.

The candidate is ready for a release only after this checklist is completed on
representative data and the remaining differences have a documented biological
or implementation explanation.
