Outputs
========

All paths below are relative to the analysis directory.  A library name is
``<method>-<condition>-<replicate>``, for example ``RIBO-treated-1``.  Selecting
a stage makes its listed targets; Snakemake also creates and retains supporting
files required to reach those targets.  See :doc:`stages` for stage selection
and input-aware stage removal.  See :doc:`table-reference` for workbook sheet,
column, contrast-direction, and missing-value semantics.

Stage output map
----------------

``trimming``
   ``trimmed/<library>.fastq`` is the mapping-ready read file.  Raw FastQC
   reports are written below ``qc/1raw/`` and trimmed-read reports below
   ``qc/2trimmed/``.  Single-end reports use ``-raw`` and ``-trimmed`` in their
   names; paired-end reports additionally use the ``_q`` and ``_p`` mate
   suffixes.

``mapping``
   ``maplink/<library>.bam`` and ``maplink/<library>.bam.bai`` contain the final
   uniquely mapped alignments after rRNA/tRNA removal and their indexes.

``qc``
   ``qc/multi/multiqc_report.html`` aggregates raw, trimmed, mapped, uniquely
   mapped, and rRNA/tRNA-removal diagnostics.

``tracks``
   Each library produces 24 BigWigs: four mapping styles, three
   normalizations, and two strands.  Their exact pattern is
   ``<mapping>tracks/<norm>/<library>.<norm>.<strand>.<mapping>.bw``, where
   ``<mapping>`` is ``global``, ``centered``, ``fiveprime``, or ``threeprime``;
   ``<norm>`` is ``raw``, ``mil``, or ``min``; and ``<strand>`` is ``forward``
   or ``reverse``.

``genome_tracks``
   The exact targets are ``tracks/potentialStartCodons.gff``,
   ``tracks/potentialAlternativeStartCodons.gff``,
   ``tracks/potentialStopCodons.gff``, and
   ``tracks/potentialRibosomeBindingSite.gff``.

``readcounts``
   The five user-facing workbooks are ``auxiliary/annotation_total.xlsx``,
   ``auxiliary/annotation_unique.xlsx``,
   ``auxiliary/total_read_counts.xlsx``,
   ``auxiliary/unique_read_counts.xlsx``, and ``auxiliary/samples.xlsx``.
   Their supporting count tables and mapped-read summaries are in
   ``readcounts/``.

``metagene``
   Each ``RIBO``, ``TIS``, or ``TTS`` library gets a directory named
   ``metageneprofiling/<library>/``.  Cross-library read-length outputs are
   ``metageneprofiling/read_length_fractions.html``,
   ``metageneprofiling/read_length_fractions.xlsx``, and
   ``metageneprofiling/read_length_counts.xlsx``.  The per-library directory
   layout is described in :doc:`metagene-profiling`.

``tis_advisor``
   Each Ribo-like library produces
   ``tis_advice/<library>/tis_recommendation.html``,
   ``tis_advice/<library>/tis_recommendation.json``, and
   ``tis_advice/<library>/read_length_evidence.tsv``.  See :doc:`tis-advisor`.

``correlation``
   ``figures/heatmap_SpearmanCorr_readCounts.pdf`` is accompanied by its source
   matrix, ``figures/SpearmanCorr_readCounts.tab``.

``pca``
   The main result is ``pca/PCA_3D.html`` and the second interactive diagnostic
   is ``pca/diffex_QC.html``.  Retained source and diagnostic files include
   ``raw_reads.csv``, ``meta.csv``, ``normalized_counts.tsv``, ``rld.tsv``,
   ``variance_percentages.tsv``, ``rld_cor.tsv``,
   ``raw_count_distributions.pdf``, and ``mean_vs_variance.pdf`` in ``pca/``.

``predictions``
   Reparation produces ``auxiliary/predictions_reparation.xlsx`` and contributes
   to ``tracks/updated_annotation.gff``.  With DeepRibo enabled,
   ``auxiliary/predictions_deepribo.xlsx`` and the accepted DeepRibo calls are
   added.  Intermediate evidence tracks, including
   ``tracks/reparation_annotated.gff`` and, when enabled,
   ``tracks/deepribo_merged.gff`` and ``tracks/deepribo_merged_plus.gff``, are
   retained for audit.

``differential_expression``
   For every configured ``<contrast>``, HRIBO creates the marker
   ``contrasts/<contrast>`` and the filtered workbooks
   ``xtail/<contrast>_sorted.xlsx``, ``riborex/<contrast>_sorted.xlsx``, and
   ``deltate/<contrast>_sorted.xlsx``.  Raw tables and diagnostic PDFs remain in
   their tool directories.  When ``overview`` is also selected, its inputs
   include the cross-contrast tables
   ``xtail/xtail_all.csv``, ``riborex/riborex_all.csv``, and
   ``deltate/deltate_all.csv``.

``overview``
   All four products are declared workflow outputs:
   ``auxiliary/overview.xlsx``, ``auxiliary/overview.tsv``,
   ``auxiliary/overview.gff``, and ``auxiliary/overview_misc.gff``.  The
   workbook and TSV combine annotation, per-library abundance and translation
   efficiency, prediction evidence, and enabled differential analyses.  The
   first GFF represents the combined CDS set; ``overview_misc.gff`` represents
   non-CDS annotation features.

Alignment links and portability
-------------------------------

``maplink/<library>.bam`` is a relative symbolic link to the corresponding
file below ``../bam/``.  Moving or archiving the complete analysis tree keeps
that link valid because the two directories retain their relative positions.
Copying only ``maplink/`` does not: the copied BAM link then has no target.
When exporting only the final BAMs, dereference the links (for example with
``cp --dereference`` or ``rsync --copy-links``) and copy the ``.bam.bai`` files
with them.

Coverage-track conventions
--------------------------

The four track mappings answer different questions:

* ``global`` adds coverage across each CIGAR-aligned reference block, excluding
  soft clips, deletions, and reference skips;
* ``centered`` clips 11 nt from both reference-alignment ends by default and
  distributes one read across the remaining aligned central positions, without
  filling CIGAR gaps;
* ``fiveprime`` assigns a read to its transcript-oriented 5' end; and
* ``threeprime`` assigns it to its transcript-oriented 3' end.

Forward-strand BigWig values are positive and reverse-strand values are
negative.  The sign is the browser convention used to display the strands on
opposite sides of zero; it does not mean that the reverse strand has negative
read abundance.  Use the absolute value when comparing coverage magnitude.

Normalization and mapped-read totals
------------------------------------

Normalization is library-wide, including for assemblies with chromosomes and
plasmids.  HRIBO's three-column mapped-read summaries retain one row per
library and contig, but consumers sum all contig rows before normalizing.

For a library-wide effective mapped total ``N``, feature length ``L``, and
feature count ``C``, workbook abundance is ``RPKM = 1e9 * C / (L * N)``.
Coverage tracks use these factors:

* ``raw``: no library-depth scaling;
* ``mil``: ``raw * 1e6 / N`` (counts per million); and
* ``min``: ``raw * N_min / N``, where ``N_min`` is the smallest complete-library
  total in the experiment.

Metagene ``cpm`` uses the same complete-library ``1e6 / N`` factor.  Its
separate ``window`` normalization is described in
:doc:`metagene-profiling`.

Every mapped alignment with a SAM ``NH`` tag contributes ``1 / NH`` to an
effective total, where ``NH`` is the number of reported hits; an alignment
without that tag contributes one.  A uniquely mapped read therefore contributes
one, while all records for a read with multiple reported hits sum to one when
the records and ``NH`` tag are complete.  This keeps
``annotation_total.xlsx`` denominators consistent with fractional
featureCounts output.  The default final BAMs and their tracks contain unique
alignments, so their ``NH`` contribution is normally one.  Fractional totals
are valid and may appear in ``readcounts/*_mapped_reads.txt``.

Count tables that look like GFF or GTF
--------------------------------------

Several files below ``readcounts/`` have historical ``.gff`` or ``.gtf``
suffixes, including ``*_annotation.gff`` and ``*_annotation.gtf``.  They are
HRIBO internal count tables: nine annotation columns followed by one count
column per library, with an ``#hribo-gff-read-counts-v1`` schema comment.  They
are deliberately **not** valid nine-column GFF3/GTF files and should not be
loaded into a genome browser or passed to a general GFF/GTF parser.

For browser-ready annotation, use the nine-column files under ``tracks/`` or
``auxiliary/overview.gff`` and ``auxiliary/overview_misc.gff``.  Use the
workbooks or ``auxiliary/overview.tsv`` for tabular downstream analysis.
