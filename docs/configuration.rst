Configure an analysis
=====================

Start every project with a fresh copy of ``config/config.yaml``.  Edit the
copy in the analysis directory and keep all top-level sections, even when a
particular stage is not selected.  Paths may be absolute or relative to the
analysis directory.

Settings every user should review
---------------------------------

.. list-table:: Required project choices
   :header-rows: 1
   :widths: 34 66

   * - Setting
     - What to provide
   * - ``biologySettings.genome``
     - Reference nucleotide FASTA, optionally gzip-compressed.
   * - ``biologySettings.annotation``
     - GFF3 or GTF annotation from the same assembly as the FASTA.
   * - ``biologySettings.samples``
     - Tab-separated library sheet described in :doc:`samples`.
   * - ``biologySettings.adapter*``
     - Adapter sequences used by the library preparation, or an empty string
       only when no adapter sequence should be supplied.
   * - ``workflowSettings.stages``
     - The analyses and results to request; see :doc:`stages`.

The input portion of a project might look like this:

.. code-block:: yaml

   biologySettings:
     adapterS3: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
     adapterS5: ""
     adapterP3R1: ""
     adapterP5R1: ""
     adapterP3R2: ""
     adapterP5R2: ""
     genome: "data/genome.fa"
     annotation: "data/annotation.gff3"
     samples: "config/samples.tsv"
     alternativeStartCodons: ["GTG", "TTG"]

Adapter settings
----------------

Use ``adapterS3`` and ``adapterS5`` for single-end libraries.  Paired-end
libraries use ``adapterP3R1``/``adapterP5R1`` for read 1 and
``adapterP3R2``/``adapterP5R2`` for read 2.  A setting may contain one
sequence or a comma-separated list of IUPAC nucleotide sequences.

An empty string means that no adapter is passed for that end.  It is not a
placeholder for an unknown sequence.  Single- and paired-end adapter fields
may coexist in one configuration when the sample sheet contains both layouts.

Alternative start codons are supplied as a YAML list.  They control the
alternative-start genome track; canonical ``ATG`` starts and ``TAG``, ``TGA``,
and ``TAA`` stops are reported separately.

Choose which results to create
------------------------------

``workflowSettings.stages`` accepts a YAML list or one of two presets:

.. code-block:: yaml

   workflowSettings:
     stages:
       - mapping
       - qc
       - tracks

``preprocessing`` requests trimming, mapping, and MultiQC.  ``full`` requests
all stages, including differential expression, and therefore requires a valid
matched Ribo-seq/RNA-seq design.  The checked-in template requests the usual
results but leaves differential expression commented out.

For a one-off run, the command line can override the YAML without changing it:

.. code-block:: console

   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --config stages=mapping,qc,tracks \
       --cores 8 all

See :doc:`stages` for every accepted name and its result files.

ORF prediction
--------------

``predictionSettings.deepribo`` accepts ``on`` or ``off``.  REPARATION is
always the base predictor when ``predictions`` is selected; enabling DeepRibo
adds a second prediction workbook and adds accepted calls to the updated
annotation.

``deepriboASiteOffset`` is the nucleotide distance from a read's 3' end to the
ribosomal A-site.  The template value of ``12`` was derived from the published
DeepRibo *E. coli* setup.  Check it for the organism, nuclease, and protocol
being analysed.  It is not a P-site offset.  The :doc:`tis-advisor` evaluates
both read ends and may provide a separate, advisory-only DeepRibo A-site
suggestion for RIBO libraries.  Compare the advice across libraries before
changing this one global setting and rerunning predictions; HRIBO never
applies it automatically.

When DeepRibo runs, the reference sequence may contain only uppercase ``A``,
``C``, ``G``, ``T``, and ``N``.  Other ambiguity symbols or lowercase sequence
cause preflight to stop because DeepRibo cannot encode them safely.

Differential expression and translation
---------------------------------------

Select the ``differential_expression`` stage to run xTail, RiboRex, and
deltaTE.  Configure it under ``differentialExpressionSettings``:

``features``
   Case-sensitive feature types from annotation column 3 to count, by default
   ``CDS`` and ``sRNA``.  HRIBO warns and skips a requested type that is absent
   while counting the requested types that are present.  It stops with an
   actionable error if none of them occur in the annotation.

``contrasts``
   A list such as ``["Treated-Control"]``.  Results use left-minus-right
   direction: a positive log2 fold change means higher signal in ``Treated``.
   An empty list requests every pairwise comparison between eligible matched
   conditions.

``padjCutoff``
   Adjusted-p-value threshold used to populate the filtered workbook sheets.
   It must be greater than 0 and less than 1.

``log2fcCutoff``
   Non-negative absolute log2-fold-change threshold used for the filtered
   ``up`` and ``down`` sheets.  Values equal to either boundary are included.

``detectionMinCPM``, ``detectionMinCount``, and ``detectionMinReplicates``
   Thresholds for the cross-condition detection report, with defaults of 1
   count per million, 10 raw reads, and 2 biological replicates.  RNA and
   RIBO are evaluated separately.  A replicate passes only when it meets
   both count and CPM thresholds; a feature is ``detected`` if enough
   replicates pass, ``not_detected`` if none pass and enough usable replicates
   were available, and ``uncertain`` otherwise.  Zero-depth samples have
   undefined CPM and are not usable.  These are evidence thresholds,
   not statistical tests.  They do not change xTail, RiboRex, or deltaTE.

``xtailBins`` and ``xtailMinMeanCount``
   xTail density resolution and minimum mean RNA/RPF count.  Higher bin counts
   take longer.  Change these only as part of a documented analysis choice.

The required library design is described in :doc:`samples`.  The visual
cross-condition report and its browser tracks are explained in
:doc:`differential-summary`; workbook fields and sheet names are explained in
:doc:`table-reference`.

Read lengths and metagene profiles
----------------------------------

Read-length fields accept comma-separated values and inclusive ranges, for
example ``22,23,27,34-35`` or ``25-34``.

``readstatSettings.readLengths`` selects the lengths shown in read-length
statistics.  ``metageneSettings`` separately controls the start/stop windows,
read lengths, filters, mapping methods, normalizations, and plot formats used
for metagene profiles.

.. list-table:: Common metagene choices
   :header-rows: 1
   :widths: 32 68

   * - Setting
     - Accepted values or meaning
   * - ``mappingMethods``
     - ``fiveprime``, ``threeprime``, ``centered``, and/or ``global``.
   * - ``readLengths``
     - Lengths included in the profile, for example ``25-34``.
   * - ``positionsOutsideORF`` / ``positionsInORF``
     - Numbers of flanking and within-CDS nucleotides shown around each start
       or stop boundary.
   * - ``normalizationMethods``
     - ``raw``, ``cpm``, and/or ``window``.
   * - ``filteringMethods``
     - ``overlap``, ``length``, and/or ``rpkm``.
   * - ``neighboringGenesDistance``
     - Distance used by the overlap filter.
   * - ``lengthCutoff`` / ``rpkmThreshold``
     - Minimum CDS length and abundance used by the corresponding filters.
   * - ``outputFormats``
     - ``interactive``, ``svg``, ``pdf``, ``png``, and/or ``jpg``.
   * - ``includePlotlyJS``
     - ``integrated`` for standalone HTML, ``online`` for smaller HTML that
       needs internet access, or ``local`` for a separately supplied script.
   * - ``colorList``
     - Optional series colours in read-length order; leave empty for the
       built-in colour-blind-friendly palette.

Review the biological meaning of these choices in
:doc:`metagene-profiling`.  ``tisAdvisorSettings`` independently selects the
read lengths and the ``fiveprime`` and/or ``threeprime`` ends evaluated for
P-site advice:

.. code-block:: yaml

   tisAdvisorSettings:
     readLengths: "22-40"
     mappingMethods: ["fiveprime", "threeprime"]

Validate changes
----------------

Always perform a dry-run after changing the configuration:

.. code-block:: console

   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --dry-run all

HRIBO reports unknown settings inside the named configuration sections,
invalid choices, unavailable contrasts, and incompatible inputs before
executing the workflow.
