Configuration
=============

Start from ``config/config.yaml`` and edit a copy in the analysis directory.
The authoritative machine-readable contract is
``workflow/schemas/config.schema.yaml``; unknown nested settings are rejected.

Biological inputs
-----------------

``biologySettings`` contains adapter sequences and the three input paths:

* ``genome``: nucleotide FASTA, optionally gzip-compressed at input;
* ``annotation``: GFF3 or GTF from the same assembly as the FASTA; and
* ``samples``: the tab-separated sheet described in :doc:`samples`.

``alternativeStartCodons`` is a YAML list such as ``["GTG", "TTG"]``.  These
codons control the alternative-start genome track; canonical ATG starts and
TAG/TGA/TAA stops are reported separately.

Adapter values may be empty or contain one or more comma-separated IUPAC
nucleotide sequences.  Single-end data use ``adapterS3``/``adapterS5``;
paired-end data use the four ``adapterP*`` values.

Differential expression
-----------------------

Differential analysis is enabled by selecting the ``differential_expression``
stage, not by an on/off setting.  Its configuration includes:

``features``
   Annotation feature types included in the count matrix, by default ``CDS``
   and ``sRNA``.

``contrasts``
   A YAML list such as ``["Treated-Control"]``.  The direction is left minus
   right: positive log2 fold changes mean higher signal in ``Treated``.  An
   empty list requests every pairwise combination among conditions that have
   both ``RIBO`` and ``RNA`` libraries.  Explicit contrasts are preferable for
   a release or biological comparison.

``padjCutoff`` and ``log2fcCutoff``
   Positive filtering thresholds used for the sorted result workbooks.  The
   boundaries are inclusive: adjusted p-values equal to ``padjCutoff`` are
   retained, as are fold changes equal to ``log2fcCutoff`` or its negative.

``xtailBins``
   Number of probability-density bins used by xTail.  The default 10,000 is
   computationally expensive but retains upstream behavior.

``xtailMinMeanCount``
   Minimum mean RNA and footprint count retained by xTail.  The default ``1``
   preserves HRIBO's historical inclusion boundary rather than xTail 1.2.0's
   higher default.

ORF prediction
--------------

``predictionSettings.deepribo`` accepts ``on`` or ``off``.  Reparation remains
the base predictor; enabling DeepRibo adds its model, calls, workbook, and
accepted calls to the combined annotation.

``deepriboASiteOffset`` is the nucleotide distance from a read's 3' end to the
ribosomal A-site.  The default ``12`` came from the DeepRibo E. coli setup and
should be checked for a different organism, nuclease, or protocol.  It is not
the P-site offset produced by the TIS advisor; see :doc:`tis-advisor`.

When DeepRibo is enabled, the reference sequence may contain only uppercase
``A``, ``C``, ``G``, ``T``, and ``N``.  HRIBO rejects lowercase sequence and
other IUPAC ambiguity codes during preflight because the pinned DeepRibo parser
is case-sensitive and cannot encode the broader alphabet.  Lowercase and
broader IUPAC nucleotide symbols remain valid when DeepRibo is not selected.

Read lengths and metagenes
--------------------------

Read-length specifications accept comma-separated values and inclusive ranges,
for example ``22,23,27,34-35`` or ``25-34``.

``metageneSettings`` controls the upstream/outside and downstream/inside
windows, annotation filters, mapping methods, read lengths, normalization, and
plot formats.  Valid mapping methods are ``fiveprime``, ``threeprime``,
``centered``, and ``global``.  Valid normalizations are ``raw``, ``cpm``, and
``window``.  These names differ from the BigWig track normalizations
``raw``/``mil``/``min`` described in :doc:`outputs`.

``includePlotlyJS`` selects an embedded, online, or local Plotly JavaScript
resource.  ``integrated`` makes standalone reports and embeds Plotly once per
HTML page; ``online`` makes smaller files that need internet access.

``tisAdvisorSettings`` independently selects the read lengths and the 5' and/or
3' read ends evaluated by the advisor.

Stage selection
---------------

``workflowSettings.stages`` is either ``full``, ``preprocessing``, or a YAML
list of stage names.  The checked-in default lists all ordinary stages but
comments out differential expression.  In contrast, the ``full`` preset does
include differential expression and therefore requires a valid matched design.
See :doc:`stages` before using that preset.
