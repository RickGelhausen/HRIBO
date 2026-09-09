TIS advisor
===========

The TIS advisor evaluates whether one Ribo-like library contains enough
start-codon signal to configure a translation-initiation-site caller such as
ORFBounder.  It recommends a mapped read end, a set of read lengths, and a
P-site offset for each selected length.  It does not call ORFs itself and its
recommendation is not automatically applied to HRIBO's prediction stage.

Run it by selecting the ``tis_advisor`` stage and configure the evaluated
lengths and ends under ``tisAdvisorSettings``.  Annotation filtering and profile
windows reuse ``metageneSettings``; see :doc:`configuration` and
:doc:`metagene-profiling`.

How a recommendation is made
----------------------------

For every configured read end, the advisor builds transcript-oriented
start- and stop-codon profiles from the final unique BAM.  It then evaluates
each configured read length independently:

* the read length must represent at least 0.5% of reads within the evaluated
  length range;
* the start peak must imply a plausible positive P-site distance (5--20 nt
  from a 5' end or 5--28 nt from a 3' end); and
* the peak must be at least twice its upstream background and at least five
  background standard deviations above it.

Usable lengths are ranked by peak sharpness, frame bias, and abundance.  The
strongest starts the recommendation; another length is added only when it
improves the pooled, offset-corrected initiation peak by at least 2%.  This
prevents a marginal length from diluting a clear signal.

When both ends are evaluated, the advisor prefers the end whose estimated
offset varies least across usable read lengths.  That consistency identifies
the read end the protocol defines most precisely.  Frame bias, fraction of
evaluated reads covered, and finally a substantial peak-sharpness difference
break ties.  Evidence for the other end remains in the report and JSON.

Three-nucleotide periodicity and reading-frame composition are reported, and
frame bias contributes to the confidence label.  Weak periodicity alone does
not reject a bacterial Ribo-seq length or lower that label.  A clear initiation
peak can therefore receive a usable recommendation with a warning that
sub-codon assignment remains uncertain.

Interpreting offsets
--------------------

Offsets are positive distances from the selected mapped end to the ribosomal
P-site, measured in transcript direction:

* with ``fiveprime`` mapping, move downstream from the read's 5' end by the
  reported offset; and
* with ``threeprime`` mapping, move upstream from the read's 3' end by the
  reported offset.

The geometry is strand-aware, so the same interpretation applies to plus- and
minus-strand genes.  Offsets are estimated separately per read length; do not
replace the table with a single value unless the downstream caller explicitly
requires that approximation.

The confidence label is based on the pooled initiation peak and frame bias.
Coverage and periodicity are reported separately and can add warnings, but do
not directly lower the label.  Read those warnings as part of the result.  In
particular, an off-frame dominant signal can mean that the offset or annotated
start positions are shifted, and low covered fraction means that most evaluated
reads were not selected.

No recommendation is a valid result
-----------------------------------

``confidence: none`` with an empty ``read_lengths`` list means that neither
evaluated end supplied a trustworthy initiation peak.  It is deliberately not
converted into a guessed offset.  Check the read-length distribution,
rRNA/tRNA depletion, annotation start sites, configured filtering, and whether
the library is actually Ribo-seq.  A failed biological signal does not by
itself mean that the workflow failed.

Output files
------------

For ``<library>``, HRIBO writes:

``tis_advice/<library>/tis_recommendation.html``
   Human-facing verdict, pasteable ORFBounder-style configuration, comparison
   of the evaluated read ends, per-length evidence tables, and diagnostic
   figures.

``tis_advice/<library>/tis_recommendation.json``
   Complete machine-readable result.  ``chosen_read_end`` is ``null`` when
   neither end is usable.  ``recommendation`` contains the selected lengths,
   per-length offsets, confidence, rationale, warnings, pooled metrics, and
   covered fraction.  ``read_ends`` retains the scores and recommendation for
   every evaluated end.

``tis_advice/<library>/read_length_evidence.tsv``
   One row per evaluated length for the chosen end, or for the first configured
   end when no recommendation exists.  It includes abundance, peak and
   background measurements, offset, frame fractions, periodicity, a ``usable``
   flag, and rejection reasons.

TIS offsets are not DeepRibo offsets
------------------------------------

``predictionSettings.deepriboASiteOffset`` has a different target and
coordinate convention.  DeepRibo uses one distance from the read's 3' end to
the ribosomal **A-site** when constructing its occupancy input.  The TIS
advisor reports a read-length-specific distance from its chosen 5' or 3' end
to the **P-site**.  The sites are one codon apart, the measured end may differ,
and read length enters any conversion.  Consequently, copying an advisor
offset into ``deepriboASiteOffset`` is not valid; inspect or calibrate the
DeepRibo A-site setting separately for the organism and protocol.
