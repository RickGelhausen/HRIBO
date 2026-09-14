Compare conditions at a glance
==============================

The cross-condition report provides a visual overview when several conditions
are compared with the same control.  It is part of the
``differential_expression`` stage and does **not** require ``overview`` or
ORF predictions.  Open ``diffex_summary/condition_overview.html`` in a web
browser after the run.  The report is local to your analysis directory; no
data need to be uploaded.

Add ``differential_expression`` to your existing
``workflowSettings.stages`` list.  For several treatments against wildtype,
set ``differentialExpressionSettings.contrasts`` to an explicit list such as
``["Mut1-WT", "Mut2-WT", "Mut3-WT"]``.  The condition names must match the
sample sheet.  An empty contrast list instead requests all eligible pairs.
Only the annotation feature types selected by
``differentialExpressionSettings.features`` are included.  Novel predicted
ORFs are not silently labelled unchanged or absent; they are outside this
report unless they are part of the differential count matrix.

What the report shows
---------------------

The report has two complementary feature-level views:

* The **detection matrix** shows RNA and RIBO evidence separately for each
  condition.  Use it to find features detected in one condition but not
  another, including patterns shared across several conditions.  It retains
  an ``uncertain`` state when the evidence does not justify either conclusion.
* The **contrast matrix** shows RNA, RIBO, and translation-efficiency (TE)
  changes for each configured contrast.  Each cell is ``up``, ``down``,
  ``not_significant`` (displayed as "No directional call"), or ``not_tested``.
  Search and filters help narrow the
  display to a feature or a response pattern without searching across several
  workbooks.

The two views answer different questions.  ``not_detected`` means that a
feature did not meet the report's read-count and normalized-abundance
thresholds at the available sequencing depth; it does not prove biological
absence.  ``not_significant`` means that the result did not meet **both** the
configured adjusted-p-value and directional effect-size cutoffs.  It is not
proof that the adjusted p-value was nonsignificant, that there is no expression,
or that the biological effect is zero.  ``not_tested`` means there was no usable
test result for that cell.  Treat ``uncertain`` detection as unresolved rather
than as absent.

Positive log2 fold changes and ``up`` states refer to the **left** condition
in a contrast named ``<left>-<right>``.  For example, ``Treated-Control``
``up`` means higher signal in ``Treated``.  Inspect the underlying effect
size, adjusted p-value, and replicate counts before interpreting any pattern.
The report keeps RNA, RIBO, and TE separate; it does not average p-values
across methods or contrasts.  The RNA, RIBO, and TE contrast states come from
deltaTE.  xTail and RiboRex TE results remain supplemental method evidence,
not votes in a combined significance call.

Detection and differential thresholds
-------------------------------------

Detection is evaluated separately for RNA and RIBO.  CPM divides a feature's
count by the total counts assigned to the selected feature types in that
sample, then multiplies by one million; it is not a whole-transcriptome
absolute-expression measure.  A replicate passes when its raw count is at
least ``detectionMinCount`` **and** its CPM is at least ``detectionMinCPM``.
A feature is ``detected`` in a
condition if at least ``detectionMinReplicates`` replicates pass.  It is
``not_detected`` if none pass and the condition has at least that many
usable replicates.  A zero-assigned-count sample has undefined CPM and cannot pass.
The remaining cases are ``uncertain``.  The defaults are 10
reads, 1 CPM, and 2 replicates, respectively.  Record the values you used
when sharing a result.  A feature close to a boundary may change state with
sequencing depth or a different threshold.  Differential ``up`` and ``down`` use the
``padjCutoff`` and ``log2fcCutoff`` values in
``differentialExpressionSettings``.  These statistical states are separate
from detection.

Files for analysis and genome browsing
--------------------------------------

The report is accompanied by three tab-separated files:

``diffex_summary/condition_matrix.tsv``
   Per-feature RNA and RIBO detection by condition, with supporting count
   information.  Use this to reproduce or filter the detection matrix in a
   spreadsheet or script.

``diffex_summary/contrast_matrix.tsv``
   Per-feature RNA, RIBO, and TE states by contrast, alongside the available
   differential statistics.  Use it when you need the exact values behind a
   coloured cell.

``diffex_summary/browser_tracks.tsv``
   An index of the generated browser tracks and the condition or contrast,
   assay, and state represented by each file.

``diffex_summary/browser/*.gff3`` contains separate GFF3 tracks for detected
RNA and RIBO features in each condition and for RNA, RIBO, and TE ``up`` and
``down`` features in each contrast.  Load the tracks of interest alongside
``maplink/<library>.bam`` or the normalized BigWig coverage in your genome
browser.  A feature missing from a track is not, by itself, proof of absence:
check the report or TSV for ``uncertain``, ``not_significant``, and
``not_tested`` states and for the thresholds used.  The accompanying
``diffex_summary/browser/tracks_manifest.json`` lists every generated track
and any feature IDs that could not be placed in a browser track because their
coordinates were unavailable.  An empty GFF3 track is a valid result.

GFF3 coordinates are one-based and inclusive.  Each feature carries its
identifier, name, assay, state, and condition or contrast in the attributes;
contrast tracks include available ``log2fc`` and ``padj`` values.  The
filenames include a short hash to keep similarly named conditions distinct,
so use the track index or manifest to choose the right file.

The browser tracks are visual selections, not replacement annotations.  Keep
the original reference annotation and the per-contrast xTail, RiboRex, and
deltaTE workbooks for detailed statistical review.  See :doc:`outputs` for
the other result locations and :doc:`table-reference` for table conventions.
