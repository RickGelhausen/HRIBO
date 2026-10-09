Metagene profiling
==================

HRIBO aggregates uniquely mapped reads around annotated CDS start and stop
codons, separately by library, contig, read length, mapping method, and
normalization.  The profiles are always oriented in transcript direction, so
plus- and minus-strand genes can be interpreted on the same axis.

Begin with ``metageneprofiling/read_length_fractions.html`` to identify the
dominant fragment lengths across libraries.  Then open
``metageneprofiling/<library>/<normalization>/interactive_metagene_profiling.html``
to compare start and stop profiles. For libraries with several profiled
contigs, this page is an index linking to a separate interactive report for
each contig. Each contig page contains its configured mapping methods and plot
views, plus candidate and support counts restricted to that contig. A library
with only one profiled contig opens directly into its plots. Use the adjacent
Excel workbooks when exact values are needed.

The report keeps three complementary views: an enrichment heatmap, the
original-style overlaid read-length lines, and individual start-profile
panels. Overlaid lines show position on the x axis and the selected raw,
CPM, or window-normalized values on the y axis, with a common scale for the
start and stop panels. Each configured read length has the same colour and
dash style in both panels; clicking its legend entry toggles both lines.
Individual panels also show every configured length. Their y axes share a
zero-based linear range by default, allowing amplitudes to be compared.
In the interactive report, the ``Shared y-axis`` and ``Independent y-axes``
buttons switch between comparable amplitudes and a separate scale for each
profile. Weak lengths may appear flatter on the shared scale; use independent
scales to inspect their shapes. Neither choice changes the plotted counts or
normalization. Static exports use the shared scale. The sORF report keeps both
line views, using only the start anchor.

When ``metageneSettings.sorfMaxLength`` is positive, an additional start-only
profile is written under
``metageneprofiling/<library>/sorfs/<normalization>/``. It uses annotated CDSs
shorter than that exclusive nucleotide limit; the shipped template uses
``300``. Setting ``0``, or leaving the field out of an older configuration,
disables this additional group. The general start/stop profiles retain their
existing filters and coordinate conventions.

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

The separate :doc:`tis-advisor` shifts stop-profile coordinates by +3 nt so its
stop/A-site offsets use the first stop-codon base as coordinate 0, matching the
start/P-site convention.  The metagene outputs described here retain the
CDS-boundary axis above.

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

Annotation filtering
--------------------

Only CDS features contribute.

CDS records with identical contig, start, end, and strand are consolidated into
one candidate before filtering and counting. Their feature identifiers are
joined in the candidate table. This prevents duplicate annotation records from
counting the same CDS and its read contributions more than once in either
group.

The configured filters are applied before aggregation:

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

Separate sORF profiles
----------------------

The sORF group uses the same supplied reference annotation as the general
profile. Predictor calls are not added to its input automatically. Eligibility
is based on the annotated CDS span, ``end - start + 1``, including the annotated
stop codon. A span equal to the configured ``sorfMaxLength`` is excluded. The
default ``<300 nt`` definition describes this analysis group; it is not an
exact 100-amino-acid boundary or a universal definition of a small protein.

Within this group, the ``length`` filter is omitted even when selected in
``filteringMethods``. Neither ``lengthCutoff`` nor ``positionsInORF`` therefore
imposes a minimum CDS length. The configured ``overlap`` and ``rpkm`` filters
still apply, using the same annotation context and complete-library read
denominator as the general group.

The profile retains the fixed start axis from ``-positionsOutsideORF`` through
``positionsInORF - 1`` and the same selected read lengths, mapping methods,
normalizations, and output formats. For a CDS shorter than
``positionsInORF``, the positive-coordinate window extends beyond its stop
codon into downstream sequence. Read these profiles as nucleotide windows
around initiation sites, not as coverage constrained to the short coding
sequence. Only start profiles are plotted for this group.
The complete start window must fit within the contig; a stop window is not
required for sORF eligibility. The general profile and advisor continue to
require both start and stop windows to fit.

Each group records eligible and retained CDSs together with the CDSs that
actually contribute reads. The interactive report and figure titles include
these counts, so an aggregate can be interpreted alongside its underlying
annotation and read support.

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

For TTS libraries, inspect the stop enrichment and the advisor's A-site
recommendation.  Start enrichment in a termination library is a diagnostic
feature and does not determine its recommended offsets; see :doc:`tis-advisor`.

A broad, weak, shifted, or multi-modal aggregate is a reason to inspect the
library, but it is not by itself a diagnosis.  Common explanations include
residual rRNA or tRNA, a small number of extremely abundant features dominating
the aggregate, heterogeneous nuclease protection, mixed library populations,
or inaccurate CDS boundaries.  First compare ``raw`` profiles and the
read-length count workbook; normalization can make a sparse shape visually
prominent.  Then inspect the implicated alignments and genes in the BAM and
coverage tracks, and review trimming, depletion, mapping, and annotation
quality in MultiQC.

Outputs
-------

For every ``<library>`` and requested ``<normalization>``, the directory
``metageneprofiling/<library>/<normalization>/`` contains:

* ``<mapping>_readcounts_start.xlsx`` and
  ``<mapping>_readcounts_stop.xlsx``;
* ``interactive_metagene_profiling.html`` when ``interactive`` output is
  enabled, plus sibling ``interactive_metagene_profiling_contig_*.html`` pages
  when several contigs are profiled; and
* one ``<contig>_<mapping>.<format>`` heatmap, plus a per-read-length start
  profile figure and an overlaid read-length figure when data are available,
  for each requested static format. The additional filenames include
  ``(per read length)`` and ``(overlaid read lengths)``, respectively.

Static figures preserve their individual layout heights so that multi-panel
profiles retain space for the data and support labels.

Contig-page filenames contain a sanitized contig name and a stable identifier
to avoid filename collisions. Open the main index to follow the links, and
keep the sibling pages with it when moving or sharing the report directory.
Integrated JavaScript is embedded once in each contig page. ``local`` mode
writes a shared ``plotly.min.js`` beside the reports for offline use.

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

When the sORF group is enabled,
``metageneprofiling/<library>/sorfs/<normalization>/`` contains its start-count
workbooks and requested start figures. The interactive report uses the same
``interactive_metagene_profiling.html`` filename. Its start coordinate and
read-length columns follow the conventions above. It uses the same per-contig
navigation when several contigs are profiled. A mapping with no profile
evidence retains its explicit zero-profile page and the aggregate support
counts that explain it.

At each group root, ``metageneprofiling/<library>/`` and, when enabled,
``metageneprofiling/<library>/sorfs/``, three tab-separated tables provide
support independent of the chosen display normalization:

* ``candidates.tsv`` records annotated CDS eligibility and exclusion reasons;
* ``candidate_counts.tsv`` summarizes input, size-eligible, and retained CDS
  counts by contig and mapping method; and
* ``candidate_support.tsv`` records retained CDSs, CDSs contributing reads, and
  raw profile counts by contig, mapping method, anchor, and configured read
  length.

The count and support tables also include ``[all contigs]`` rows summarizing
the complete library. CDS counts refer to the distinct coordinate candidates
after duplicate annotation records have been consolidated.

These tables distinguish a short-CDS group with few eligible annotations from
a retained group whose selected lengths have little or no read support. See
:doc:`table-reference` for their field conventions.

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
contig. The candidate tables show whether a zero sORF profile follows from
size eligibility, annotation filtering, or a lack of reads at the selected
lengths.
