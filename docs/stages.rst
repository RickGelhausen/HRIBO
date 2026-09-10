Choose analyses and results
===========================

A stage tells HRIBO which result to produce.  Snakemake automatically runs the
prerequisites, so selecting ``predictions`` also performs the trimming,
mapping, filtering, and counting needed by the predictors.  Stages are not
execution barriers and do not need to be listed in dependency order.

Stage catalogue
---------------

.. list-table:: Available stages
   :header-rows: 1
   :widths: 21 39 40

   * - Stage
     - Question it helps answer
     - Main result
   * - ``trimming``
     - Were adapters removed and are the processed reads usable?
     - Mapping-ready reads in ``trimmed/`` and raw/processed FastQC reports.
   * - ``mapping``
     - Where do the filtered reads align uniquely?
     - ``maplink/<library>.bam`` and its ``.bam.bai`` index.
   * - ``qc``
     - How did read quality, mapping, and rRNA/tRNA depletion perform?
     - ``qc/multi/multiqc_report.html``.
   * - ``tracks``
     - What does strand-aware coverage look like across the genome?
     - Raw and normalized BigWigs in the four ``*tracks/`` directories.
   * - ``genome_tracks``
     - Where are possible start, stop, alternative-start, and RBS motifs?
     - Four browser-ready GFF files in ``tracks/``.
   * - ``readcounts``
     - How many reads and how much normalized signal belong to each feature?
     - Annotation, count-summary, and sample workbooks in ``auxiliary/``.
   * - ``metagene``
     - Which read lengths and start/stop patterns characterize each library?
     - Interactive reports, figures, and workbooks in ``metageneprofiling/``.
   * - ``tis_advisor``
     - Which mapped end, read lengths, and P-site offsets have usable evidence?
     - HTML, JSON, and TSV reports in ``tis_advice/<library>/``.
   * - ``correlation``
     - Do related libraries have similar binned genomic coverage profiles?
     - Spearman correlation heatmap in ``figures/``.
   * - ``pca``
     - Do samples group by condition and assay as expected?
     - Interactive PCA and diagnostics in ``pca/``.
   * - ``predictions``
     - Which annotated or novel ORFs are supported by the predictors?
     - REPARATION and optional DeepRibo workbooks plus
       ``tracks/updated_annotation.gff``.
   * - ``differential_expression``
     - Which features change at RNA, footprint, or translation-efficiency level?
     - Per-contrast xTail, RiboRex, and deltaTE workbooks.
   * - ``overview``
     - How can annotation, abundance, predictions, and differential evidence be
       reviewed together?
     - ``auxiliary/overview.xlsx`` plus TSV and GFF views.

Common selections
-----------------

Use a preset for the two most common broad choices:

.. code-block:: yaml

   workflowSettings:
     stages: "preprocessing"

``preprocessing`` selects ``trimming``, ``mapping``, and ``qc``.

.. code-block:: yaml

   workflowSettings:
     stages: "full"

``full`` selects all thirteen stages.  It includes differential expression and
therefore needs matched Ribo-seq/RNA-seq libraries in at least two conditions.

For a smaller custom analysis, provide a list.  Some useful patterns are:

.. list-table:: Example stage recipes
   :header-rows: 1
   :widths: 33 67

   * - Goal
     - Stage list
   * - Mapping and browser inspection
     - ``mapping``, ``qc``, ``tracks``
   * - Ribo-seq quality and metagene profiling
     - ``qc``, ``tracks``, ``metagene``, ``tis_advisor``
   * - Feature abundance only
     - ``readcounts``
   * - ORF discovery and combined table
     - ``predictions``, ``overview``
   * - Matched differential analysis
     - ``qc``, ``pca``, ``correlation``, ``differential_expression``,
       ``overview``

``overview`` is not a lightweight summary-only stage.  It builds the
prediction results it combines, including DeepRibo when enabled.  When
differential expression is selected, it also collects those contrast results.
Plan the same container downloads and compute resources as for those analyses.

For example:

.. code-block:: yaml

   workflowSettings:
     stages:
       - qc
       - tracks
       - metagene
       - tis_advisor

Override stages for one run
---------------------------

The top-level command-line value ``stages`` takes precedence over the YAML:

.. code-block:: console

   $ /path/to/HRIBO/run_hribo.sh \
       --configfile "$PWD/config/config.yaml" \
       --config stages=mapping,qc,tracks \
       --cores 8 all

Use an override for exploration or a one-off continuation.  Record the command
with the analysis, because the effective stage list will differ from the
configuration file.

Stages and sample types
-----------------------

HRIBO prints the resolved stage list at startup.  If the sample sheet contains
no ``RIBO`` library, it removes ``predictions``,
``differential_expression``, and ``overview``.  If it also contains no
``TIS`` or ``TTS`` library, it removes ``metagene`` and ``tis_advisor``.

Input validation follows the requested work.  For example, ``genome_tracks``
needs the genome but not FASTQ files, while results downstream of mapping need
the genome, annotation, and reads.  See :doc:`samples` for stage-specific
design requirements and :doc:`outputs` for the full result locations.
