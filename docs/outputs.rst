Understand the results
======================

HRIBO writes results into the analysis directory from which it was launched.
The selected :doc:`stages` determine which result groups are present.  A
library is named ``<method>-<condition>-<replicate>``; for example,
``RIBO-Treated-1``.

What to inspect first
---------------------

For a typical analysis, review the results in this order:

1. Open ``qc/multi/multiqc_report.html`` and confirm that the reads, trimming,
   mapping, and rRNA/tRNA depletion are plausible for every library.
2. Check the correlation heatmap and PCA when multiple libraries are present.
   Biological replicates should normally resemble each other, while unexpected
   grouping can reveal a sample-label or quality problem.
3. Review read-length distributions and start/stop metagene profiles for each
   Ribo-like library.  Read the TIS-advisor confidence and warnings together
   with its recommended offsets.
4. Load the final BAM, BigWig, and GFF files in a genome browser to inspect
   individual loci.
5. Use ``auxiliary/overview.xlsx`` as the entry point for feature-level
   abundance, prediction evidence, and differential results.  Consult the
   predictor- or method-specific workbooks before accepting a candidate.

Primary result map
------------------

.. list-table:: Results by selected stage
   :header-rows: 1
   :widths: 20 35 45

   * - Stage
     - Main result
     - What the result contains
   * - ``trimming``
     - ``trimmed/<library>.fastq``
     - Adapter-trimmed or assembled reads that continue to mapping.  Raw and
       processed FastQC reports are below ``qc/1raw/`` and ``qc/2trimmed/``.
   * - ``mapping``
     - ``maplink/<library>.bam`` and ``.bam.bai``
     - Final uniquely mapped alignments after rRNA/tRNA removal, ready for a
       genome browser or downstream alignment tools.
   * - ``qc``
     - ``qc/multi/multiqc_report.html``
     - One HTML report aggregating the main read-processing and mapping
       diagnostics for all libraries.
   * - ``tracks``
     - ``globaltracks/``, ``centeredtracks/``, ``fiveprimetracks/``, and
       ``threeprimetracks/``
     - Strand-separated BigWig coverage with raw and normalized views.
   * - ``genome_tracks``
     - ``tracks/potential*.gff``
     - Browser tracks for possible start codons, alternative starts, stop
       codons, and ribosome-binding-site motifs.
   * - ``readcounts``
     - ``auxiliary/annotation_*.xlsx`` and
       ``auxiliary/*_read_counts.xlsx``
     - Per-feature annotation, abundance, direct TE ratios, and summarized
       feature counts using total or unique mappings.
   * - ``metagene``
     - ``metageneprofiling/read_length_fractions.html`` and
       ``metageneprofiling/<library>/``
     - Read-length composition plus start- and stop-centred aggregate profiles,
       source tables, and figures.
   * - ``tis_advisor``
     - ``tis_advice/<library>/tis_recommendation.html``
     - Recommended mapped end, usable read lengths, per-length P-site offsets,
       confidence, warnings, evidence, and machine-readable companions.
   * - ``correlation``
     - ``figures/heatmap_SpearmanCorr_readCounts.pdf``
     - Pairwise Spearman correlation of binned genomic coverage and its source
       matrix.
   * - ``pca``
     - ``pca/PCA_3D.html`` and ``pca/diffex_QC.html``
     - Interactive sample separation and normalized-count diagnostics.
   * - ``predictions``
     - ``auxiliary/predictions_reparation.xlsx`` and optionally
       ``auxiliary/predictions_deepribo.xlsx``
     - ORF coordinates, predictor evidence, abundance, sequence, and annotation
       context.  ``tracks/updated_annotation.gff`` combines accepted calls with
       the supplied annotation.
   * - ``differential_expression``
     - ``xtail/<contrast>_sorted.xlsx``,
       ``riborex/<contrast>_sorted.xlsx``, and
       ``deltate/<contrast>_sorted.xlsx``
     - Per-feature RNA, footprint, and/or translation-efficiency effects,
       p-values, adjusted p-values, and prefiltered up/down sheets.
   * - ``overview``
     - ``auxiliary/overview.xlsx``
     - Consolidated annotation, sequences, abundance, direct TE, predictor
       evidence, and enabled differential statistics.  The same main table is
       available as ``auxiliary/overview.tsv``; two GFF files provide browser
       views.

A full run has a result layout similar to this (supporting directories are
omitted):

.. code-block:: text

   my-analysis/
   ├── auxiliary/
   │   ├── overview.xlsx
   │   ├── overview.tsv
   │   ├── annotation_unique.xlsx
   │   ├── predictions_reparation.xlsx
   │   └── predictions_deepribo.xlsx       # only when enabled
   ├── qc/multi/multiqc_report.html
   ├── maplink/<library>.bam
   ├── globaltracks/ ... threeprimetracks/
   ├── tracks/updated_annotation.gff
   ├── metageneprofiling/<library>/
   ├── tis_advice/<library>/
   ├── figures/heatmap_SpearmanCorr_readCounts.pdf
   ├── pca/PCA_3D.html
   ├── xtail/ ... riborex/ ... deltate/
   └── logs/

Quality control and alignments
------------------------------

The MultiQC report combines results from several points in the workflow.  Use
it to compare raw and trimmed read quality, read counts and lengths, mapping
rates, unique alignments, and the effect of rRNA/tRNA filtering.  A completed
HTML file only means the report was generated; it does not mean that every
library passed biological quality control.

``maplink/<library>.bam`` is the convenient location for the final unique BAM.
It is a relative symbolic link to the corresponding file under ``bam/``.  Move
or archive the entire analysis tree to preserve the link.  If exporting only
the final BAMs, copy with link dereferencing and include each ``.bam.bai``
index.

Coverage and genome-browser tracks
----------------------------------

Each library receives BigWigs for four mapping views, three normalization
choices, and two strands.  Their name pattern is:

.. code-block:: text

   <mapping>tracks/<normalization>/<library>.<normalization>.<strand>.<mapping>.bw

The mapping views answer different questions:

``global``
   Coverage across every aligned reference block of the read.

``centered``
   Coverage from the central aligned portion after clipping the read ends.

``fiveprime`` and ``threeprime``
   One-position views anchored at the transcript-oriented 5' or 3' read end.

The normalizations are ``raw`` (no depth scaling), ``mil`` (counts per
million mapped reads), and ``min`` (all libraries scaled to the smallest
mapped library).  Use ``raw`` to inspect the evidence in one library and a
normalized track when comparing libraries.

Forward-strand BigWig values are positive and reverse-strand values are
negative so genome browsers can display the strands on opposite sides of zero.
The negative sign is a display convention, not negative abundance.

The browser-ready motif files are:

* ``tracks/potentialStartCodons.gff``;
* ``tracks/potentialAlternativeStartCodons.gff``;
* ``tracks/potentialStopCodons.gff``; and
* ``tracks/potentialRibosomeBindingSite.gff``.

Counts, abundance, and direct TE
--------------------------------

``auxiliary/annotation_unique.xlsx`` uses uniquely mapped reads;
``auxiliary/annotation_total.xlsx`` includes fractional contributions from
multi-mapped reads.  Both organize annotated features into sheets and add one
RPKM column per library.  When a Ribo-like library has a matching RNA control
with the same condition and replicate, the workbooks also contain a direct
translation-efficiency ratio.

``auxiliary/unique_read_counts.xlsx`` and
``auxiliary/total_read_counts.xlsx`` provide smaller feature-class summaries.
Use them for an overview of where reads were counted; use the annotation
workbooks when you need individual features and sequences.

RPKM values use the mapped total across the complete library, including all
contigs.  A direct ``*_TE`` value is the Ribo-like RPKM divided by its matched
RNA-like RPKM.  It is a ratio, not a log2 fold change.

Metagene profiles and TIS advice
--------------------------------

``metageneprofiling/read_length_fractions.html`` compares fragment-length
composition across libraries.  Each per-library directory then contains
start- and stop-centred profiles for the requested read ends, normalizations,
and plot formats.  Use :doc:`metagene-profiling` to interpret the axes,
normalizations, peaks, and valid zero profiles.

For TIS advice, begin with
``tis_advice/<library>/tis_recommendation.html``.  It provides the human-facing
recommendation and diagnostic plots.  The adjacent JSON preserves the complete
machine-readable result, and ``read_length_evidence.tsv`` provides one row per
evaluated length.  No recommendation can be a valid result when the data do not
contain a trustworthy initiation peak; see :doc:`tis-advisor`.

ORF predictions
---------------

``auxiliary/predictions_reparation.xlsx`` contains REPARATION calls.
``auxiliary/predictions_deepribo.xlsx`` is present only when DeepRibo is
enabled.  Review coordinates, predictor scores or probabilities, the libraries
contributing evidence, RPKM/TE values, and available annotation metadata.  Use
``auxiliary/overview.xlsx`` for the explicit overlapping-gene context.
Predictor support is evidence to evaluate, not by itself proof of translation.

``tracks/updated_annotation.gff`` combines the supplied annotation with
accepted prediction calls for browser inspection.  Keep the separate
prediction workbooks when reporting which predictor supported an ORF.

Differential results
--------------------

For a contrast named ``<left>-<right>``, positive log2 fold changes mean higher
signal in the left condition and negative values mean higher signal in the
right condition.  For example, a positive value in ``Treated-Control`` means
higher signal in ``Treated``.

Each tool-specific workbook contains an ``all`` sheet and filtered sheets such
as ``TE_up`` and ``TE_down``.  deltaTE additionally separates RNA, RIBO, and TE
changes.  The filtered sheets use ``padjCutoff`` and ``log2fcCutoff`` from the
configuration; always inspect the effect size and adjusted p-value together.
The three tools model translation differently, so review agreement and
disagreement rather than treating one column as interchangeable across tools.

Combined overview
-----------------

``auxiliary/overview.xlsx`` is the most convenient feature-level starting
point.  Its ``all`` sheet joins annotated and predicted CDS-like coordinates
with sequence, per-library RPKM and direct TE, prediction evidence, and any
enabled differential results.  Other sheets provide annotated or feature-type
subsets.

``auxiliary/overview.tsv`` contains the same rows and columns as the ``all``
sheet for scripted analysis.  ``auxiliary/overview.gff`` represents the
combined CDS set and ``auxiliary/overview_misc.gff`` contains non-CDS features
for genome-browser use.

Detailed workbook sheet names, column definitions, missing values, and
coordinate conventions are documented in :doc:`table-reference`.

Primary results versus supporting files
---------------------------------------

Snakemake retains inputs and intermediate files needed to build the selected
results.  Directories such as ``readcounts/``, ``bam/``, and tool-specific work
areas are useful for troubleshooting and custom downstream analysis, but most
users should begin with the primary files listed above.  Do not assume a file
with a historical ``.gff`` or ``.gtf`` suffix below ``readcounts/`` is
browser-ready; use GFF files from ``tracks/`` or ``auxiliary/`` instead.
