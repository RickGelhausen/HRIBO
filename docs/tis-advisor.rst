TIS/TTS advisor
===============

The advisor recommends a mapped read end, usable read lengths, and a
read-length-specific ribosomal-site offset for each Ribo-like library:

* ``RIBO`` and ``TIS`` libraries are evaluated at annotated **start codons**
  and receive **P-site** advice for a translation-initiation-site caller such
  as ORFBounder.
* ``TTS`` libraries are evaluated at annotated **stop codons** and receive
  **A-site** advice for translation-termination peaks.  Start profiles remain
  available as diagnostics, but do not determine the TTS recommendation.

Initiation complexes place the start codon in the P-site, while
apidaecin-stalled termination complexes place the stop codon in the A-site;
see `Weaver et al. (2019) <https://pmc.ncbi.nlm.nih.gov/articles/PMC6401488/>`_
and `Florin et al. (2017) <https://www.nature.com/articles/nsmb.3439>`_.
The advisor does not call ORFs itself or apply recommendations automatically.
For RIBO libraries it can also give separate DeepRibo A-site offset guidance
when the 3' read-length evidence supports one shared value.

For routine review, open
``tis_advice/<library>/tis_recommendation.html`` first.  TTS libraries show
``TTS peak advice``; other Ribo-like libraries show ``TIS caller advice``.
The report summarizes the verdict, confidence, warnings, per-length evidence,
and diagnostic plots.  Use the JSON or TSV companions for downstream scripts.
The historical stage, settings, directory, and filenames are retained for
compatibility.

Run it by selecting the ``tis_advisor`` stage and configure the evaluated
lengths and ends under ``tisAdvisorSettings``.  Annotation filtering and profile
windows reuse ``metageneSettings``; see :doc:`configuration` and
:doc:`metagene-profiling`.  The workflow passes the sample-sheet method to the
advisor.  For standalone use, the method is inferred from the BAM's
``<METHOD>-<CONDITION>-<REPLICATE>`` name; an unrecognized prefix defaults to
``RIBO``.  Use ``--library_type TTS`` to evaluate stop peaks when a standalone
BAM has a different naming convention.

The supplied CDS intervals must include the complete stop codon.  Stop
positions are taken from the terminal three CDS bases; a GTF that lists the
stop separately and excludes it from the CDS must be converted to that input
convention before interpreting termination offsets.

How a recommendation is made
----------------------------

The advisor filters the annotation once, then builds transcript-oriented
start- and stop-codon profiles from the final unique BAM for every configured
read end.  All evaluated ends use the same retained CDSs.  The RPKM filter
counts footprints that overlap the CDS for every library type, so initiation
5' ends upstream of the start and termination 3' ends downstream of the stop
do not exclude an otherwise expressed gene.  Profile counts still use the
requested 5' or 3' end.  The selected boundary is evaluated independently for
each configured read length:

* the read length must represent at least 0.5% of reads within the evaluated
  length range;
* the peak must imply a plausible positive distance to the calibrated site;
  and
* the peak must be at least twice its upstream background and at least five
  background standard deviations above it.

Start/P-site distances must be 5--20 nt from a 5' end or 5--28 nt from a
3' end.  Stop/A-site distances must be 8--23 nt from a 5' end or 2--25 nt
from a 3' end for shorter footprints.  For TTS footprints of at least 50 nt,
the 5' search instead targets the leading terminating ribosome near the
footprint's 3' end: its distance must be between ``read_length - 26`` and
``read_length - 3``.  This permits disome-length footprints without establishing
that any particular long read represents a disome.  Include their lengths in
``tisAdvisorSettings.readLengths`` explicitly; the range is not expanded
automatically.  Disome-based TTS profiling is demonstrated by
`Froschauer et al. (2025)
<https://www.nature.com/articles/s41467-025-58329-w>`_.

Background is measured after converting read-end coordinates to calibrated
site coordinates: add the offset for 5' mapping or subtract it for 3' mapping.
The start/P-site background spans -100 through -30 nt; the stop/A-site
background spans -100 through -40 nt, inside the coding region.  Comparing
stop peaks with coding coverage helps avoid treating ordinary elongation
coverage as a termination peak merely because downstream coverage is low.
Ribosome queues can raise this background and make TTS advice more conservative.

At least ten observed background positions must be available.  During pooling,
zero-filled positions introduced by offset correction are excluded from the
background wherever any contributing profile lacks an observed position.
The enrichment divisor, ``background_reference``, is the background median
when positive, otherwise the background mean, and finally the observed-profile
mean if the entire background is zero.  Per-length evidence and pooled metrics
use this same background definition in site coordinates.  Their numerical
ratios can differ from older reports that measured background in unshifted
read-end coordinates.

Selection balances supported read coverage with boundary enrichment. Among
the individually usable lengths, the advisor chooses an anchor that maximizes
``read_fraction * log(1 + enrichment)``. Read fractions come from accepted
mapped reads within the evaluated length range, while enrichment is the
offset-corrected boundary peak relative to its observed background. Ranking
by read count instead of fraction gives the same anchor because the denominator
is shared across lengths. The logarithm reduces the influence of
an extreme enrichment ratio at a rare length. It does not cap the reported
enrichment or interpret all reads of that length as translated footprints.

The enrichment floor is the larger of ``2`` and
``minRelativeEnrichment * anchor_enrichment``. The setting defaults to ``0.5``
and can range from 0 to 1. Each included length must meet this floor, and the
pooled, offset-corrected profile must also meet it, remain at least five
background standard deviations above background, and retain at least ten
observed background positions. Starting with the anchor, the advisor adds the
feasible length carrying the most additional reads, reevaluating the remaining
lengths after each addition. It stops when no addition preserves these quality
requirements. Pooling may lower enrichment while including more supported
reads; an increase in enrichment is no longer required.

This is a deterministic greedy selection heuristic, rather than an exhaustive
search or a validated biological enrichment cutoff. The anchor combines
abundance and enrichment; the relative floor expresses how much anchor
enrichment can be lost in exchange for read coverage. A rare extreme ratio or
an abundant marginal signal can still influence that anchor. Inspect the
reported anchor, floor, selected read fraction, and pooled quality alongside
the metagene profiles and replicates. Increasing the setting makes the
relative criterion stricter; setting it to zero retains the absolute quality
requirements. Reading-frame and periodicity measurements remain separate
diagnostics and confidence evidence.

When both ends are evaluated, the advisor prefers the end whose estimated
offset varies least across usable read lengths.  That consistency identifies
the read end the protocol defines most precisely.  Supported frame-0 fraction,
fraction of evaluated reads covered, and finally a substantial peak-sharpness difference
break ties.  Evidence for the other end remains in the report and JSON.

Three-nucleotide periodicity and reading-frame composition are measured
inside the coding body: downstream of a start or upstream of a stop,
excluding the boundary peak.  Only positions observed within the configured
output axis for every contributing profile are used; alignment padding cannot
create frame or periodicity evidence.  The number of read ends supporting these frame
fractions is reported as ``frame_reads``.  High confidence requires a pooled
peak at least five times background, at least 50% of coding-body read ends
in the annotated frame 0, and at least 30 coding-body read ends contributing
to the frame fractions.  A dominant off-frame signal cannot promote a
recommendation to high confidence.  Weak
periodicity alone does not reject a bacterial Ribo-seq length or lower the
confidence label.  A clear boundary peak with sparse or weak frame evidence
can still receive usable advice, with a warning that its offsets rest on the
boundary peak alone.

The heatmap's scored boundary panel uses the exact ``background_reference``
shown in the evidence table.  The other panel is labeled diagnostic and is
normalized by its whole-window median (or mean when the median is zero).
Hover over either panel for the enrichment value and raw read-end counts.
Compare the scored panel with
the table when judging recommended lengths; the diagnostic panel uses a
different denominator.

Interpreting the 3-nt FFT score
-------------------------------

The ``periodicity`` field is a descriptive spectrum ratio from the raw,
count-weighted metagene, rather than the fraction of correctly positioned
reads.  Each retained CDS contributes its read-end counts; genes are not
given equal weights.  After offset calibration, a 90-nt coding-body window
is selected in P-site coordinates: +15 through +104 relative to a start,
or -105 through -16 relative to a stop.  The termination A-site offset is
converted to a P-site offset first.  Only observed positions within the
configured axis are used, so a clipped window can contain fewer than 90 bins.

For the observed count vector ``x``, the advisor subtracts its mean, computes
``abs(rfft(x - mean(x)))**2``, and divides the power in the frequency bin
nearest 1/3 cycles per nucleotide by the sum of all nonzero-frequency bins.
The implementation sums the returned one-sided squared magnitudes directly,
without doubling interior frequency bins.  It returns zero for fewer than
nine observed bins, zero total counts, a flat vector, or zero non-DC power.

This score has no significance test, RNA-seq control, or replicate calibration.
It discards the phase of the triplet signal and therefore cannot distinguish
annotated frame 0 from frames 1 and 2; the separate frame fractions provide
that information.  Coverage gradients, individual pauses, and spectral
leakage when a clipped window is not a multiple of three affect the ratio.
A score of 13% is not a 13% probability of translation or correct nucleotide
assignment.  The below-15% notice is an internal descriptive heuristic,
not a published bacterial quality or accuracy threshold.  It does not alter
length selection, offsets, or confidence.

Weak metagene periodicity is common in bacterial MNase-based Ribo-seq, and
even a visible triplet pattern can arise from nuclease sequence preference
and codon composition.  Stronger reading-frame signals are possible with
protocols such as RelE-assisted profiling, with their own sequence biases.
See `Mohammad et al. (2019)
<https://pmc.ncbi.nlm.nih.gov/articles/PMC6377232/>`_ and
`Hwang and Buskirk (2017)
<https://academic.oup.com/nar/article/45/1/327/2290904>`_.
These observations do not make weak body periodicity a universal failure
criterion.  In TIS/TTS libraries, boundary peak position and shape,
read-length offset consistency, and agreement between replicates provide
distinct evidence for start/P-site or stop/A-site calibration.  An offset
calibrated on stalled boundary ribosomes is not, by itself, validation of
nucleotide assignments to all elongating ribosomes.

Interpreting offsets
--------------------

Offsets are positive distances from the inclusive mapped read end to the
**first nucleotide of the site codon**, measured in transcript direction.
The calibrated site is the P-site for RIBO/TIS advice and the A-site for TTS
advice:

* with ``fiveprime`` mapping, move downstream from the read's 5' end by the
  reported offset; and
* with ``threeprime`` mapping, move upstream from the read's 3' end by the
  reported offset.

The geometry is strand-aware, so the same interpretation applies to plus- and
minus-strand genes.  Offsets are estimated separately per read length; do not
replace the table with a single value unless the downstream caller explicitly
requires that approximation.  Published offsets depend on the experimental
protocol and coordinate convention, so they are examples rather than universal
defaults.

In the advisor's figures, coordinate 0 is the first base of the selected start
or stop codon.  For stop plots this adds 3 nt to the historical metagene axis,
where the stop occupies -3, -2, and -1.  The separate
``metageneprofiling/`` outputs retain that historical axis; see
:doc:`metagene-profiling`.

The A-site is one codon downstream of the P-site.  TTS reports therefore show
both the directly estimated A-site offset and a derived P-site offset:

* from a 5' end, ``P-site offset = A-site offset - 3``;
* from a 3' end, ``P-site offset = A-site offset + 3``.

For RIBO/TIS estimates, the reverse conversion gives ``A = P + 3`` from a
5' end or ``A = P - 3`` from a 3' end.  These conversions assume simple
aligned footprints and codons separated by three nucleotides. TTS advice
directly estimates A-site distances and derives P-site distances.
ORFBounder JSON exports retain the directly calibrated site
and convert the offset sign as described below.

The confidence label is based on the pooled boundary peak and supported
annotated-frame-0 evidence.
Coverage and periodicity are reported separately and can add warnings, but do
not directly lower the label.  Read those warnings as part of the result.  An
off-frame dominant signal can mean that offsets or annotated boundary
positions are shifted; low covered fraction means that most evaluated reads
were not selected.

TTS interpretation
------------------

Apidaecin profiling can produce upstream ribosome queues and stop-codon
readthrough, as shown by
`Mangano et al. (2020) <https://elifesciences.org/articles/62655>`_.
Start-codon enrichment also occurs in some TTS monosome libraries and should
not be taken as the termination signal; see
`Froschauer et al. (2025)
<https://www.nature.com/articles/s41467-025-58329-w>`_.  Review stop and start
diagnostics together with the protocol, matched RIBO/TIS data when available,
and individual loci.  A pooled stop peak calibrates library geometry; it does
not by itself identify every true termination site or quantify normal
termination efficiency.

No recommendation is a valid result
-----------------------------------

``confidence: none`` with an empty ``read_lengths`` list means that neither
evaluated end supplied a trustworthy peak at the selected boundary.  It is
deliberately not converted into a guessed offset.  Check the read-length
distribution, rRNA/tRNA depletion, annotated start or stop positions,
configured filtering, and whether the protocol enriches the expected
boundary.  A failed biological signal does not by itself mean that the
workflow failed.

Output files
------------

For ``<library>``, HRIBO writes:

``tis_advice/<library>/tis_recommendation.html``
   Human-facing verdict, comparison of evaluated read ends, per-length
   evidence, diagnostic figures, and guidance for the adjacent ORFBounder JSON
   exports. TTS reports include A-site offsets and derived P-site offsets.
   RIBO reports can additionally suggest a separate DeepRibo A-site setting.

``tis_advice/<library>/tis_recommendation.json``
   Complete machine-readable result.  ``library_type`` records the resolved
   method, ``anchor`` is ``start`` or ``stop``, and ``site`` is ``P`` or ``A``.
   ``offset_reference`` states the coordinate convention.
   ``chosen_read_end`` is ``null`` when neither end is usable.
   ``recommendation`` contains selected lengths, generic ``offsets`` to the
   calibrated site, explicit ``p_site_offsets`` and ``a_site_offsets``,
   ``anchor``, ``site``, confidence, rationale, warnings, pooled metrics,
   ``frame_reads``, and covered fraction.  The other site's offsets are derived
   by the three-nucleotide conversion above.  ``read_ends`` retains scores and
   recommendations for every evaluated end.  Per-length scores include
   ``background_reference``, ``background_positions``, and ``frame_reads``
   so enrichment and frame support can be reviewed explicitly.
   ``recommendation.selection`` records the selection strategy
   (``coverage_with_quality_floor``), ``reference_read_length``,
   ``reference_abundance``, ``reference_sharpness``, and ``reference_utility``
   (the reference read fraction times ``log(1 + enrichment)``).
   ``min_relative_enrichment`` records the configured setting;
   ``minimum_sharpness`` is the resulting enrichment floor. ``eligible_fraction``
   counts reads from all individually usable lengths before that relative
   floor is applied. The selected fraction and pooled enrichment remain
   ``recommendation.covered_fraction`` and ``recommendation.sharpness``.
   These fractions use accepted reads within the evaluated length range as
   their denominator. Selection metadata is separate from calibrated offsets.
   ``deepribo_a_site`` contains separate DeepRibo advice or the reason none
   was made; ``applied`` is always
   ``false``.

``tis_advice/<library>/read_length_evidence.tsv``
   One row per evaluated length for the chosen end, or for the first configured
   end when no recommendation exists.  It includes abundance, peak and
   background measurements, ``background_reference`` (the enrichment divisor),
   ``background_positions`` (observed positions supporting the baseline),
   generic ``offset``, frame fractions, ``frame_reads``, periodicity,
   a ``usable`` flag, and rejection reasons.  ``read_end``, ``anchor``, and
   ``site`` identify the calibration; ``p_site_offset`` and ``a_site_offset``
   make the directly estimated and derived distances explicit.

ORFBounder JSON exports
-----------------------

The ``tis_advisor`` stage also prepares JSON inputs for a later ORFBounder
analysis. It does not install or execute ORFBounder, or apply its settings to
any HRIBO analysis.

Each ``tis_advice/<library>/orfbounder/`` directory contains ``manifest.json``.
For each evaluated read end with a supported recommendation, it also contains
``<read_end>/read_lengths.json`` and ``<read_end>/offsets.json``. Both ends are
exported when both have usable advice, each with its own selected lengths and
offsets; the exports are not restricted to the preferred end.

The combined ``tis_advice/orfbounder/`` directory uses the same layout and
collects library entries separately under ``fiveprime/`` and ``threeprime/``.
It preserves every library's calibration rather than averaging conditions or
replicates. Each manifest records the available recommendations, calibrated
sites, confidence, preferred read ends, and reasons for missing advice.

The consumable files contain only sample-keyed values. For example, one
five-prime export may contain:

.. code-block:: json

   {"TIS-WT-1": "28,30", "RIBO-WT-1": "28,30"}

and the corresponding offsets file:

.. code-block:: json

   {
     "TIS-WT-1": {"28": -12, "30": -13},
     "RIBO-WT-1": {"28": -12, "30": -13}
   }

These numbers illustrate the format. Read-length selections are
comma-separated strings, and each canonical decimal length key maps to an
integer nucleotide offset. ORFBounder shifts endpoints using the opposite
sign convention to the advisor's positive distances: ``fiveprime`` exports
``-distance`` and ``threeprime`` exports ``+distance``. This moves a 5' endpoint
downstream or a 3' endpoint upstream on either strand. Exports use the directly
estimated P-site distances for RIBO/TIS and the directly estimated A-site
distances for TTS, measured to the first nucleotide of that codon. The generic
advisor JSON continues to report positive distances and both P-/A-site tables.

Sample keys match the alignment filename stem before its first underscore.
HRIBO's ``METHOD-condition-replicate.bam`` names therefore retain their exact
method, condition, and replicate in both files. ORFBounder accepts one
``mapping_method`` per command or batch row, shared by all supplied assays.
Choose the matching end-specific files and calibrated input libraries manually;
the manifest is guidance, not another ORFBounder input file. Paired TIS/TTS/RIBO
inputs must all have usable advice for the selected end. ORFBounder calling
requires at least one TIS or TTS input; RIBO calibrations can accompany those
assays.

Libraries without a supported recommendation for an end are omitted from that
end's files, with their reasons retained in the manifest. If no library has
advice for an end, its two input files are absent. No fallback ``default``
entries or guessed zero offsets are written. Supplying an omitted library to
ORFBounder would require a separate reviewed calibration, because its input
validation requires every supplied calling sample to be covered by the JSONs.

Existing advisor results can be converted without rerunning metagene profiling
or the advisor. From the analysis directory, run the bundled converter with
the path to your HRIBO checkout:

.. code-block:: bash

   python3 /path/to/HRIBO/workflow/scripts/export_orfbounder_inputs.py \
     --recommendations tis_advice/*/tis_recommendation.json \
     --output-dir tis_advice/orfbounder

This writes the combined end-specific input pairs and manifest. It leaves the
source recommendation files unchanged.

DeepRibo advice is separate from boundary advice
------------------------------------------------

``predictionSettings.deepriboASiteOffset`` has a different scope from the
per-length boundary recommendation.  DeepRibo uses one distance from the
read's 3' end to the ribosomal **A-site** when constructing its occupancy
input.  For a RIBO library, the advisor's separate DeepRibo section uses all
usable *3'-end* P-site estimates.  It subtracts three nucleotides to form an
A-site candidate for each read length, then suggests one value only when those
candidates agree and cover enough reads.  The conversion remains an estimate
to review, not a direct measurement of A-site position.

The report compares any suggestion with the currently configured
``deepriboASiteOffset``.  If they differ, compare the RIBO-library reports
across conditions and replicates before changing that single global setting
and rerunning predictions.  HRIBO never applies the advice automatically.
When 3' mapping was not evaluated, the 3' signal is weak, a strong reading-frame
bias points away from the annotated frame, or lengths disagree, the report gives
no single DeepRibo setting and explains why.  Read-support percentages use reads
accepted by the advisor, which can differ from reads accepted by DeepRibo.
TIS and TTS libraries do not receive a DeepRibo suggestion because DeepRibo
uses RIBO libraries.
