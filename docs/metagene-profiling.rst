Metagene profiling
==================

HRIBO aggregates uniquely mapped reads around annotated CDS start and stop
codons, separately by library, contig, read length, mapping method, and
normalization.  The profiles are always oriented in transcript direction, so
plus- and minus-strand genes can be interpreted on the same axis.

Coordinate convention
---------------------

Let ``outside`` be ``positionsOutsideORF`` and ``inside`` be
``positionsInORF``.  Both axes increase from the transcript's 5' side toward
its 3' side; genomic coordinates therefore increase for plus-strand genes and
decrease for minus-strand genes.

Start profile
   Coordinates run from ``-outside`` through ``inside - 1``.  Negative values
   are upstream of the CDS, coordinate 0 is the first base of the annotated
   start codon, and positive values proceed into the CDS.  Under this
   convention the start codon occupies 0, 1, and 2.

Stop profile
   Coordinates run from ``-inside`` through ``outside - 1``.  Negative values
   approach the stop through the coding sequence, the annotated stop codon
   occupies -3, -2, and -1, and coordinate 0 is the first base downstream of
   the CDS.  Positive values continue into downstream sequence.

The stop axis is intentionally different from the start axis: zero is the
boundary crossed when leaving the CDS.  It is not a strand-specific reversal.
For both strands, moving right in either plot means moving downstream along the
transcript.

Mapping methods
---------------

``fiveprime``
   Add one count at the read's transcript-oriented 5' end.

``threeprime``
   Add one count at the read's transcript-oriented 3' end.

``centered``
   Add one count at the rounded midpoint of the read interval.

``global``
   Add one count at every CIGAR-aligned reference position where the read
   overlaps the profile window.  Soft-clipped query bases, deletions, and
   reference skips are not coverage.  Each aligned block is clipped to the
   window before it is placed on the transcript-oriented axis.

The ``global`` implementation now fills every valid strand/anchor combination.
Older HRIBO output could contain empty minus-strand start slices and
plus-strand stop slices because a genomic array was sliced in the wrong
direction.  Stop profiles from the old release may also look mirrored relative
to the convention above.  These are expected corrections, not biological
strand asymmetry; see :doc:`real-data-validation` when comparing releases.

Annotation filtering
--------------------

Only CDS features contribute.  The configured filters are applied before
aggregation:

``overlap``
   Remove a CDS when another same-strand CDS lies within
   ``neighboringGenesDistance`` of its interval.

``length``
   Require at least the larger of ``lengthCutoff`` and ``positionsInORF``
   nucleotides, so the requested inside window is supported.

``rpkm``
   Require at least ``rpkmThreshold`` using the selected mapping method.  The
   denominator is the complete accepted-read total across all contigs in the
   library, not the count on the CDS's own contig.

A CDS is also omitted when either its start or stop window would extend beyond
the contig.  This keeps the two anchors comparable and avoids inventing
out-of-range positions.  The profile is an aggregate over all retained CDSs;
it is not divided by the number of retained genes.

Normalization
-------------

``raw``
   Aggregated counts with no library-depth scaling.  A point-mapping read adds
   one at its selected position; ``global`` adds one at each overlapped
   position.

``cpm``
   Multiply every profile value by ``1e6 / N``, where ``N`` is the complete
   accepted-read total for the library summed across all contigs.  This makes
   chromosome and plasmid profiles use one consistent library denominator.

``window``
   Normalize each non-empty contig/read-length column and each anchor
   separately by ``column_sum / window_length``.  The resulting column has a
   mean of one (and a sum equal to the window length), emphasizing profile
   shape rather than library depth.  A zero column remains finite and all-zero
   instead of becoming ``NaN``.

The heatmap applies an additional display-only enrichment transform: each read
length is shown relative to its own median background, falling back to its
mean for a zero median.  This lets a sparse but sharp length remain visible
beside a deep one.  Use the workbooks, rather than heatmap colour, when absolute
``raw`` or ``cpm`` values are required.

P-site markers
--------------

Start-codon figures for ``fiveprime`` and ``threeprime`` profiles can mark a
statistically supported P-site offset for each read length.  A 5' marker is
drawn upstream of the start at ``-offset``; a 3' marker is drawn downstream at
``+offset``, because the reported value is the positive distance from that
mapped end to the P-site in transcript direction.  ``centered`` and ``global``
profiles do not have a physical read-end offset and therefore receive no such
marker.

Marker estimates always come from the raw start-profile counts.  Choosing
``cpm`` or ``window`` changes the presented profile values but cannot move a
marker.  A length whose estimate does not pass the significance test is left
unmarked rather than being assigned a peak from noise.

Biological interpretation
-------------------------

Use metagene profiles as protocol-level quality control, not as proof that
every annotated CDS is translated.  In a successful conventional Ribo-seq
library, protected-fragment lengths are usually enriched over other lengths
and one or more of those lengths show a reproducible read-end peak at a
consistent distance from the annotated start.  The preferred end and expected
distance depend on the library protocol and nuclease, so compare the observed
signal with the experimental design rather than assuming one universal
offset.  Start and stop profiles should be interpreted together and, when
replicates exist, the length-specific pattern should be reproducible between
replicates.

A broad, weak, shifted, or multi-modal aggregate is a reason to inspect the
library, but it is not by itself a diagnosis.  Common explanations include
residual rRNA or tRNA, a small number of extremely abundant features dominating
the aggregate, heterogeneous nuclease protection, mixed library populations,
or inaccurate CDS boundaries.  First compare ``raw`` profiles and the
read-length count workbook; normalization can make a sparse shape visually
prominent.  Then inspect the implicated alignments and genes in the BAM and
coverage tracks, and review trimming, depletion, mapping, and annotation
quality in MultiQC.

The figures from the former documentation repository are deliberately not
reused here: they were produced before the strand/anchor geometry and marker
contracts described above were corrected.  Versioned example figures and
expected values should be added only after the current release candidate has
completed the :doc:`real-data-validation` protocol.

Outputs
-------

For every ``<library>`` and requested ``<normalization>``, the directory
``metageneprofiling/<library>/<normalization>/`` contains:

* ``<mapping>_readcounts_start.xlsx`` and
  ``<mapping>_readcounts_stop.xlsx``;
* ``interactive_metagene_profiling.html`` when ``interactive`` output is
  enabled; and
* one ``<contig>_<mapping>.<format>`` heatmap, plus a per-read-length start
  profile figure when data are available, for each requested static format.

Workbook sheets are named by contig.  Their first column is ``coordinates``,
the following columns are exactly the lengths configured in
``metageneSettings.readLengths``, and ``sum`` is the row-wise total over those
lengths.  On a retained contig, a configured length that was not observed is
kept as an explicit zero column; an observed but unrequested length is
excluded.  This selection
does not redefine the complete-library denominator used by ``cpm`` or the RPKM
annotation filter.  The cross-library files
``metageneprofiling/read_length_fractions.html``,
``read_length_fractions.xlsx``, and ``read_length_counts.xlsx`` provide the
separate read-length composition summary.

Zero and ``no_evidence`` profiles
---------------------------------

An all-zero profile is valid scientific output.  It means that no accepted read
contributed to that anchor/read-length combination after filtering; it does not
on its own indicate a failed job.  When no contig has profile evidence at all,
the workbook contains an explicit ``no_evidence`` sheet with the complete
coordinate axis, every configured read-length column filled with zeros, and a
zero ``sum`` column.  When only one anchor or read length has evidence, the
missing counterpart is represented by a finite zero column so start and stop
workbooks remain structurally comparable.

Before accepting a zero result, inspect the workflow log and then check the
configured read-length range, the filtering summary, ``rpkmThreshold``, CDS
lengths and overlaps, boundary exclusions, and the mapped-read depth on that
contig.
