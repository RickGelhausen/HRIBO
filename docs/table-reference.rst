Result table reference
======================

Use this page when reading a workbook or importing an HRIBO result into another
program.  For a shorter guide to which result to open first, see
:doc:`outputs`.  Paths below are relative to the analysis directory.  A
library name has the form
``<method>-<condition>-<replicate>`` and a contrast has the form
``<left>-<right>``.

Conventions shared by feature tables
------------------------------------

Coordinates are one-based and inclusive.  ``Identifier`` is normally
``<Genome>:<Start>-<Stop>:<Strand>``; differential-expression tables retain
the identifier supplied to the statistical tool and resolve its coordinates
from the annotation (or from that coordinate form for a predicted ORF).
``Length`` is ``Stop - Start + 1``.  ``Codon_count`` is the integer quotient
``Length / 3`` in the annotation, prediction, and differential workbooks.  The
overview retains the exact quotient and can therefore contain a fractional
value for a non-triplet feature.

Sequence columns are reported in transcript orientation:

``Start_codon`` and ``Stop_codon``
   The first and last three nucleotides of ``Nucleotide_seq``.

``Upstream_15nt``
   The 15 nucleotides immediately upstream of the feature.  Differential
   workbooks do not include this column.

``Nucleotide_seq`` and ``Aminoacid_seq``
   The feature sequence and its bacterial-code (translation table 11)
   translation.  Translation stops at the first stop codon.  ``Aminoacid_seq``
   is blank when the nucleotide length is not divisible by three.

``Locus_tag``, ``Old_locus_tag``, ``Name``, ``Product``, and ``Note`` are
annotation metadata and are blank when unavailable.  ``Genome`` is the FASTA
sequence identifier, ``Feature`` is the annotation feature type, and
``Source`` is the producing annotation source.

Abundance and direct TE columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The annotation, prediction, and overview workbooks add one
``<library>_rpkm`` column per library.  For feature count ``C``, inclusive
feature length ``L``, and complete-library mapped-read total ``N``, its value
is ``1e9 * C / (L * N)``, rounded to two decimal places.  The mapped total is
summed across all contigs, so every contig in a library has the same
denominator.

Matched assay/control pairs add direct translational-efficiency columns:

* ``RIBO`` is paired with ``RNA``;
* ``TIS`` is paired with ``RNATIS``; and
* ``TTS`` is paired with ``RNATTS``.

A pair must have the same condition and replicate.  Its column is named
``<Ribo-like-method>-<condition>-<replicate>_TE`` and contains the Ribo-like
RPKM divided by the matched RNA-like RPKM; it is a ratio, not a logarithm.  If
more than one replicate is matched for a method and condition, an additional
``<Ribo-like-method>-<condition>-avg_TE`` column is the arithmetic mean of the
defined replicate ratios.  Libraries and the resulting dynamic columns are
ordered deterministically by library name.

Sample workbook
---------------

``auxiliary/samples.xlsx`` has one sheet, ``samples``.  It preserves the source
sample-sheet columns and their order.  The required columns are:

``method, condition, replicate, fastqFile, fastqFile2``

The two FASTQ columns contain basenames rather than full source paths.
``fastqFile2`` is blank for a single-end library.  Any additional source
columns accepted by validation are also retained.

Annotation workbooks
--------------------

``auxiliary/annotation_total.xlsx`` and
``auxiliary/annotation_unique.xlsx`` have the same schema.  The former uses
the total/fractional annotation counts and matching mapped totals; the latter
uses uniquely mapped counts and totals.

The sheets, in order, are ``CDS``, ``rRNA``, ``sRNA``, ``transcript``,
``5'-UTR``, ``tRNA``, ``pseudogene``, ``gene``, ``region``,
``miscellaneous``, and ``all``.  ``miscellaneous`` receives feature types not
routed to a named sheet, while ``all`` contains every input feature.  Empty
feature sheets retain their headers.

Every sheet has these columns in this order, with the two dynamic blocks
expanded in place:

.. code-block:: text

   Identifier, Genome, Source, Feature, Start, Stop, Strand,
   Locus_tag, Old_locus_tag, Name, Length, Codon_count,
   <direct TE columns>, <library RPKM columns>,
   Start_codon, Stop_codon, Upstream_15nt, Nucleotide_seq,
   Aminoacid_seq, Product, Note

``Source`` is ``HRIBO`` in these two workbooks.

Read-count summary workbooks
----------------------------

``auxiliary/total_read_counts.xlsx`` and
``auxiliary/unique_read_counts.xlsx`` each contain a ``Main`` sheet with:

.. code-block:: text

   Orientation, Class, Feature_count, <one column per library>

Rows are, in order, ``ncRNA``, ``sRNA``, ``5'-UTR``, ``CDS``, ``rRNA``,
``tRNA``, ``transcript``, ``pseudogene``, and ``total``.  ``Orientation`` is
``sense``.  ``Feature_count`` is the number of counted annotation records in
the class, and each library column is the sum of its feature counts.  The
``total`` row sums the listed classes.  As above, the two workbooks differ by
total/fractional versus unique mappings.

Prediction workbooks
--------------------

Both prediction workbooks contain a single ``CDS`` sheet.  Dynamic TE and
RPKM blocks follow the rules above.

``auxiliary/predictions_reparation.xlsx`` columns are:

.. code-block:: text

   Identifier, Genome, Source, Feature, Start, Stop, Strand,
   Reparation_probability,
   Locus_tag, Old_locus_tag, Name, Length, Codon_count,
   <direct TE columns>, <library RPKM columns>, Evidence,
   Start_codon, Stop_codon, Upstream_15nt, Nucleotide_seq, Aminoacid_seq

``Reparation_probability`` is Reparation's ``prob`` attribute and ``Source``
is ``reparation``.

``auxiliary/predictions_deepribo.xlsx`` columns are:

.. code-block:: text

   Identifier, Genome, Source, Feature, Start, Stop, Strand,
   Deepribo_score, Deepribo_rank, Novel_rank,
   Locus_tag, Old_locus_tag, Name, Length, Codon_count,
   <direct TE columns>, <library RPKM columns>, Evidence,
   Start_codon, Stop_codon, Upstream_15nt, Nucleotide_seq, Aminoacid_seq

``Deepribo_score`` is the prediction value, ``Deepribo_rank`` is the GFF score
field, and ``Novel_rank`` ranks candidates that are not already annotated.

Overview workbook
-----------------

``auxiliary/overview.xlsx`` joins annotation, predictions, abundance, and any
enabled differential results.  ``auxiliary/overview.tsv`` is the same table
as its ``all`` sheet.

Sheets
~~~~~~

``all``
   The union of annotated CDS-like entries and Reparation/DeepRibo prediction
   coordinates.  For this table, CDS-like annotation comprises ``CDS``,
   ``sRNA``, ``rRNA``, ``tRNA``, and ``ncRNA``.

``annotated``
   The subset of ``all`` having at least one of ``Locus_tag``,
   ``Old_locus_tag``, ``Name``, or ``Gene_name``.

Additional populated non-CDS sheets
   Recognized feature types use ``ncRNA``, ``sRNA``, ``rRNA``, ``tRNA``,
   ``start_codons``, ``stop_codon``, ``repeat_region``, ``3'-UTR``, and
   ``5'-UTR``; any other included non-CDS type uses ``misc``.  Unlike the
   annotation workbooks, an empty category is not emitted.  Gene,
   pseudogene, exon, and CDS records are excluded from this non-CDS pass.
   Consequently an ``sRNA``, ``rRNA``, ``tRNA``, or ``ncRNA`` entry can occur
   both in ``all`` and in its feature-specific sheet.

``all`` and ``annotated`` columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The columns are:

.. code-block:: text

   Identifier, Genome, Start, Stop, Strand, Locus_tag, Overlapping_genes,
   Old_locus_tag, Name, Gene_name, Length, Codon_count, Start_codon,
   Stop_codon, Upstream_15nt, Nucleotide_seq, Aminoacid_seq,
   <direct TE columns>, <library RPKM columns>,
   Evidence_reparation, Reparation_probability,
   Evidence_deepribo, Deepribo_rank, Deepribo_score,
   <differential columns>, Product, Note

``Overlapping_genes`` is a comma-separated list of locus tags for overlapping
same-strand annotated CDS-like intervals.  Prediction evidence tokens are
prefixed with ``reparation-`` or ``deepribo-`` when the source did not already
provide that prefix.

Non-CDS sheet columns
~~~~~~~~~~~~~~~~~~~~~

The non-CDS sheets use:

.. code-block:: text

   Identifier, Genome, Start, Stop, Strand, Feature, Locus_tag,
   Old_locus_tag, Name, Gene_name, Length, Codon_count, Start_codon,
   Stop_codon, Upstream_15nt, Nucleotide_seq, Aminoacid_seq,
   <direct TE columns>, <library RPKM columns>,
   <differential columns>, Product, Note

They do not carry the five prediction columns or ``Overlapping_genes``.

Dynamic differential columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For every resolved ``<contrast>``, the overview reserves these columns in
xTail, RiboRex, then deltaTE blocks:

.. code-block:: text

   xtail_<contrast>_TE_log2FC
   xtail_<contrast>_TE_pvalue
   xtail_<contrast>_TE_pvalue_adjusted

   riborex_<contrast>_TE_log2FC
   riborex_<contrast>_TE_pvalue
   riborex_<contrast>_TE_pvalue_adjusted

   deltaTE_<contrast>_RIBO_log2FC
   deltaTE_<contrast>_RIBO_pvalue
   deltaTE_<contrast>_RIBO_pvalue_adjusted
   deltaTE_<contrast>_RNA_log2FC
   deltaTE_<contrast>_RNA_pvalue
   deltaTE_<contrast>_RNA_pvalue_adjusted
   deltaTE_<contrast>_TE_log2FC
   deltaTE_<contrast>_TE_pvalue
   deltaTE_<contrast>_TE_pvalue_adjusted

xTail contributes its final TE estimate and final p-value to this reduced
view.  A configured contrast is kept with its given orientation and order.  If
contrasts are inferred, every pair of eligible conditions is generated in a
deterministic order.  A missing tool/feature/contrast result leaves blank
cells; the columns are not removed.

Cross-condition summary tables
------------------------------

The ``differential_expression`` stage also writes a searchable report at
``diffex_summary/condition_overview.html`` and three tab-separated companion
tables.  See :doc:`differential-summary` for the visual interpretation and
browser-track guide.  The summary is independent of the prediction-dependent
``overview`` stage.

``diffex_summary/condition_matrix.tsv`` has one feature/condition/assay row
for RNA or RIBO detection, with columns in this order:

.. code-block:: text

   feature_id, genome, start, end, strand, feature_type, name,
   condition, assay, state, replicates, passing_replicates,
   mean_count, mean_cpm, sample_counts, sample_cpms

``state`` is ``detected``, ``not_detected``, or ``uncertain``.  A replicate
passes when both its raw count and CPM meet the configured minimums.
``detected`` requires at least ``detectionMinReplicates`` passing replicates;
``not_detected`` requires none passing and at least that many usable
replicates.  Other cases are ``uncertain``.  ``replicates`` counts samples
with a non-zero assay-wide count total; zero-depth samples have undefined
CPM and are not counted.  ``passing_replicates`` counts those meeting both
thresholds.  ``sample_counts`` and ``sample_cpms`` contain semicolon-separated
``<sample>=<value>`` entries so individual replicate values are not hidden by
the means.

``diffex_summary/contrast_matrix.tsv`` has one feature/contrast/assay row,
with columns in this order:

.. code-block:: text

   feature_id, genome, start, end, strand, feature_type, name,
   contrast, assay, state, log2fc, padj, method,
   xtail_te_log2fc, xtail_te_padj,
   riborex_te_log2fc, riborex_te_padj

``assay`` is RNA, RIBO, or TE.  ``state`` is ``up``, ``down``,
``not_significant``, or ``not_tested``.  The ``log2fc`` direction is left
minus right.  Up/down calls require both ``padj <= padjCutoff`` and a fold
change at least ``log2fcCutoff`` in the matching direction.  The three
primary assay states use deltaTE statistics; the xTail and RiboRex columns
are supplemental TE-only results, not combined p-values or extra votes in
the ``state`` field.  ``not_significant`` is the machine-readable name for
"no directional call under both thresholds"; its adjusted p-value may still
be below the configured cutoff when the effect-size boundary is not met.
``not_tested`` means an effect estimate or adjusted p-value is missing.

``diffex_summary/browser_tracks.tsv`` indexes the generated GFF3 tracks:

.. code-block:: text

   path, kind, condition, contrast, assay, direction, feature_count

``kind`` separates condition-detection from contrast-change tracks.
``feature_count`` is the number of selected features in that track.  A
missing feature in a GFF3 track may be ``uncertain``, ``not_significant``, or
``not_tested`` rather than absent; look up its TSV row before interpreting
the omission.  ``diffex_summary/browser/tracks_manifest.json`` also records
the track files and feature IDs that lack browser coordinates.

Differential-expression workbooks
---------------------------------

For a contrast ``<left>-<right>``, every fold-change column has left-minus-
right direction: a positive value means higher signal in ``left`` and a
negative value means higher signal in ``right``.  ``RIBO`` and ``RPF`` denote
footprint abundance, ``RNA``/``mRNA`` denote transcript abundance, and ``TE``
denotes their relative change.  Differential ``TE_log2FC`` columns are log2
effects and must not be confused with the non-logarithmic direct ``*_TE``
ratios described above.

All three workbooks begin with:

.. code-block:: text

   Identifier, Genome, Start, Stop, Strand, Locus_tag, Old_locus_tag, Name

and end with:

.. code-block:: text

   Length, Codon_count, Start_codon, Stop_codon,
   Nucleotide_seq, Aminoacid_seq

RiboRex
~~~~~~~

``riborex/<contrast>_sorted.xlsx`` has ``all``, ``TE_up``, and ``TE_down``
sheets.  Between the common blocks it contains:

.. code-block:: text

   baseMean, log2FC, log2FC_SE, stat, pvalue, pvalue_adjusted

These are the standardized names for the RiboRex/DESeq2 base mean,
translation-efficiency log2 fold change, standard error, Wald statistic,
unadjusted p-value, and adjusted p-value.  ``all`` is sorted first by
``pvalue_adjusted`` and then by location.

xTail
~~~~~

``xtail/<contrast>_sorted.xlsx`` has ``all``, ``TE_up``, and ``TE_down``
sheets.  Between the common blocks it contains:

.. code-block:: text

   mRNA_log2FC, RPF_log2FC,
   log2FC_TE_v1, pvalue_v1,
   log2FC_TE_v2, pvalue_v2,
   log2FC_TE_final, pvalue_final, pvalue_adjusted

``mRNA_log2FC`` and ``RPF_log2FC`` are the two abundance effects.
``log2FC_TE_v1`` compares their fold-change distributions;
``log2FC_TE_v2`` compares the condition-specific RPF/mRNA-ratio
distributions.  xTail selects the estimate paired with the larger, more
conservative, of ``pvalue_v1`` and ``pvalue_v2`` as
``log2FC_TE_final``/``pvalue_final``; a tie selects v2.  HRIBO uses xTail's
default Benjamini-Hochberg adjustment for ``pvalue_adjusted``.  ``all`` is
sorted first by ``pvalue_adjusted`` and then by location.

deltaTE
~~~~~~~

``deltate/<contrast>_sorted.xlsx`` has ``all``, ``RNA_up``, ``RNA_down``,
``RIBO_up``, ``RIBO_down``, ``TE_up``, and ``TE_down`` sheets.  Between the
common blocks it contains:

.. code-block:: text

   RIBO_baseMean, RIBO_log2FC, RIBO_log2FC_SE,
   RIBO_pvalue, RIBO_pvalue_adjusted,
   RNA_baseMean, RNA_log2FC, RNA_log2FC_SE,
   RNA_pvalue, RNA_pvalue_adjusted,
   TE_baseMean, TE_log2FC, TE_log2FC_SE, TE_stat,
   TE_pvalue, TE_pvalue_adjusted

Only identifiers present in all three deltaTE component tables are retained.
The ``all`` sheet is sorted first by ``TE_pvalue_adjusted`` and then by
location.

Up/down sheets
~~~~~~~~~~~~~~

``TE_up`` and ``TE_down`` in the RiboRex and xTail workbooks use ``log2FC``
and ``log2FC_TE_final``, respectively.  Each deltaTE component sheet uses the
matching component's ``*_log2FC`` and ``*_pvalue_adjusted``.  A row is included
when its adjusted p-value is at most ``padjCutoff`` and its fold change is at
least ``log2fcCutoff`` (up) or at most its negative (down).  Both boundaries
are inclusive.  An empty selection remains a header-only sheet.

xTail diagnostic plots
----------------------

For ``<left>-<right>``, HRIBO supplies the right condition to xTail first as
``control`` and the left condition second as ``treated``.  xTail compares the
second condition with the first, which produces the left-minus-right signs
used in its workbook.

``xtail/fc_<contrast>.pdf``
   The x-axis is ``mRNA_log2FC`` and the y-axis is ``RPF_log2FC``; both are
   log2(left/right).  Blue is mRNA-only (labelled ``transcription only``), red
   is RPF-only (labelled ``translation only``), green is above cutoff on both
   axes with the same sign, yellow is above cutoff on both axes with opposite
   signs, and gray is below cutoff on both axes.  Here ``translation only``
   means an RPF abundance fold change, not a direct TE fold change.

``xtail/r_<contrast>.pdf``
   The x-axis is right-condition log2(RPF/mRNA) and the y-axis is
   left-condition log2(RPF/mRNA).  xTail labels these axes
   ``control_log2TE`` and ``treated_log2TE``.  Blue means right-only, red means
   left-only, green means both axes exceed the cutoff with the same sign,
   yellow means both exceed it with opposite signs, and gray means neither
   does.  This is a condition-level TE scatter, not a plot of
   ``log2FC_TE_final``.  The geometric difference ``y - x`` has the same
   left-minus-right interpretation, but xTail's modelled final estimate need
   not equal that simple arithmetic difference.

Both plots use xTail 1.2.0's default absolute cutoff of 1 on each axis.  This
is independent of HRIBO's configured ``log2fcCutoff``, which filters workbook
sheets only.  Genes with incomplete xTail rows are omitted; complete stable
genes remain as gray points.  Point colour represents only the categories
above, not p-value.

Metagene workbooks
------------------

The cross-library ``metageneprofiling/read_length_counts.xlsx`` and
``metageneprofiling/read_length_fractions.xlsx`` have one sheet per observed
contig.  Their first column is ``read_lengths`` and each remaining column is a
library name.  Counts are accepted alignments of each observed configured
query length.  Fractions divide each library/contig length count by that
library/contig's sum; a zero-total column remains zero.

Within ``metageneprofiling/<library>/<normalization>/``, every requested
mapping method has ``<mapping>_readcounts_start.xlsx`` and
``<mapping>_readcounts_stop.xlsx``.  Sheets are retained contigs, or one
``no_evidence`` sheet when none has evidence.  Columns are:

.. code-block:: text

   coordinates, <one column per configured read length>, sum

``sum`` is the row-wise total over read-length columns.  Start coordinates run
from ``-positionsOutsideORF`` through ``positionsInORF - 1``; stop coordinates
run from ``-positionsInORF`` through ``positionsOutsideORF - 1``.  Values use
the containing directory's ``raw``, ``cpm``, or ``window`` normalization.
Configured but unobserved lengths remain explicit zero columns.

Undefined and empty values
--------------------------

Zero is a measured value unless a field is explicitly documented as a
sentinel.  The following cases are undefined or absent:

* A direct TE ratio whose matched RNA-like RPKM is zero is written as the
  ``NaN`` marker.  Replicate averages ignore undefined ratios; if none is
  defined, the average is also ``NaN``.  Spreadsheet readers may display or
  import that marker as an empty/missing value.
* Statistical tools can emit ``NA``/``NaN`` for an unestimable fold change,
  standard error, statistic, or p-value.  Such rows can remain in ``all`` but
  cannot satisfy an up/down filter.  Pooled CSV and overview cells may expose
  them as empty values.
* Missing overview differential results are blank.  In the CDS-like overview,
  absence of a Reparation call is currently represented by probability ``0``;
  absence of a DeepRibo call is represented by rank ``999999`` and score
  ``0``.  Their evidence cells are blank.  These values are compatibility
  sentinels and must not be interpreted as measured predictor output.
* Missing annotation text and unavailable sequence-derived values are blank.
  Header-only differential sheets and all-zero metagene sheets are valid
  outputs, not evidence that a workbook failed to generate.
