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

For every configured read end, the advisor builds transcript-oriented start-
and stop-codon profiles from the final unique BAM.  It evaluates the selected
boundary independently for each configured read length:

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

For starts, background is measured 100--30 nt upstream of the first start
base.  For stops, it is measured 100--40 nt upstream of the first stop base,
inside the coding region.  Comparing with coding coverage helps avoid treating
ordinary elongation coverage as a termination peak merely because downstream
coverage is low.  Ribosome queues can raise this background and make TTS
advice more conservative.
The peak search window is excluded from the stop background, including when
long 5' footprints place it within the nominal upstream interval.
The TTS annotation RPKM filter counts footprints
that overlap the CDS, even when their selected read end lies downstream of
the stop; profile counts still use the requested 5' or 3' end.

Usable lengths are ranked by peak sharpness, frame bias, and abundance.  The
strongest starts the recommendation; another length is added only when it
improves the pooled, offset-corrected boundary peak by at least 2%.  This
prevents a marginal length from diluting a clear signal.

When both ends are evaluated, the advisor prefers the end whose estimated
offset varies least across usable read lengths.  That consistency identifies
the read end the protocol defines most precisely.  Frame bias, fraction of
evaluated reads covered, and finally a substantial peak-sharpness difference
break ties.  Evidence for the other end remains in the report and JSON.

Three-nucleotide periodicity and reading-frame composition are measured
inside the coding body: downstream of a start or upstream of a stop,
excluding the boundary peak.  Frame bias contributes to the confidence label.
Weak periodicity alone does not reject a bacterial Ribo-seq length or lower
that label.  A clear boundary peak can therefore receive a usable
recommendation with a warning that sub-codon assignment remains uncertain.

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
aligned footprints and codons separated by three nucleotides.  TTS advice
uses an A-site table; it does not present termination offsets as pasteable
ORFBounder ``psiteOffsets``.

The confidence label is based on the pooled boundary peak and frame bias.
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
   evidence, and diagnostic figures.  RIBO/TIS reports include pasteable
   ORFBounder-style configuration; TTS reports include A-site offsets and
   derived P-site offsets.  RIBO reports can additionally suggest a separate
   DeepRibo A-site setting.

``tis_advice/<library>/tis_recommendation.json``
   Complete machine-readable result.  ``library_type`` records the resolved
   method, ``anchor`` is ``start`` or ``stop``, and ``site`` is ``P`` or ``A``.
   ``offset_reference`` states the coordinate convention.
   ``chosen_read_end`` is ``null`` when neither end is usable.
   ``recommendation`` contains selected lengths, generic ``offsets`` to the
   calibrated site, explicit ``p_site_offsets`` and ``a_site_offsets``,
   ``anchor``, ``site``, confidence, rationale, warnings, pooled metrics,
   and covered fraction.  The other site's offsets are derived by the
   three-nucleotide conversion above.  ``read_ends`` retains scores and
   recommendations for every evaluated end.  ``deepribo_a_site`` contains
   separate DeepRibo advice or the reason none was made; ``applied`` is always
   ``false``.

``tis_advice/<library>/read_length_evidence.tsv``
   One row per evaluated length for the chosen end, or for the first configured
   end when no recommendation exists.  It includes abundance, peak and
   background measurements, generic ``offset``, frame fractions, periodicity,
   a ``usable`` flag, and rejection reasons.  ``read_end``, ``anchor``, and
   ``site`` identify the calibration; ``p_site_offset`` and ``a_site_offset``
   make the directly estimated and derived distances explicit.

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
