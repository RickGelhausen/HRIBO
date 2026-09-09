Workflow stages
===============

A stage describes the deliverables requested from Snakemake.  It is not an
execution barrier: Snakemake also builds every prerequisite.  Requesting only
``predictions`` therefore still stages the references, trims, maps, filters,
and counts the reads needed by the predictors.

.. list-table:: Stage catalogue
   :header-rows: 1
   :widths: 22 78

   * - Stage
     - Requested deliverables
   * - ``trimming``
     - Mapping-ready trimmed or assembled reads plus raw and processed FastQC.
   * - ``mapping``
     - Final uniquely mapped, rRNA/tRNA-depleted BAMs and indexes in ``maplink/``.
   * - ``qc``
     - Aggregated ``qc/multi/multiqc_report.html``.
   * - ``tracks``
     - Global, centered, 5', and 3' BigWig coverage tracks.
   * - ``genome_tracks``
     - GFF tracks for start, stop, alternative-start, and ribosome-binding motifs.
   * - ``readcounts``
     - Annotation, read-count, and sample workbooks in ``auxiliary/``.
   * - ``metagene``
     - Per-library metagene profiles and cross-library read-length summaries.
   * - ``tis_advisor``
     - Per-library P-site/read-end recommendation reports and evidence.
   * - ``correlation``
     - Spearman read-count correlation heatmap.
   * - ``pca``
     - PCA and supporting normalized-count diagnostics.
   * - ``predictions``
     - Reparation, optional DeepRibo, prediction workbooks, and
       ``tracks/updated_annotation.gff``.
   * - ``differential_expression``
     - xTail, Riborex, and deltaTE results for each configured contrast.
   * - ``overview``
     - The combined ``auxiliary/overview.xlsx`` plus TSV and browser-GFF views.

Presets and overrides
---------------------

``preprocessing`` selects ``trimming``, ``mapping``, and ``qc``.  ``full``
selects all 13 stages, including differential expression.  A one-run override
does not require editing the YAML:

.. code-block:: console

   $ snakemake ... all --config stages=mapping,tracks,qc
   $ snakemake ... all --config stages=preprocessing

The top-level command-line ``stages`` value takes precedence over
``workflowSettings.stages``.

Input-aware behavior
--------------------

HRIBO always validates the configuration and sample-sheet structure before
constructing a DAG.  It validates the contents of only the input classes needed
by the requested stages and their dependencies: ``genome_tracks`` needs the
genome but not annotation or FASTQ files, while ``trimming`` needs the FASTQ
files but not the references.  Every other stage descends from mapping and
therefore needs the genome, annotation, and reads.  When the sample sheet has no
``RIBO`` library, HRIBO removes
``predictions``, ``differential_expression``, and ``overview``.  When it also
has no ``TIS`` or ``TTS`` library, it removes ``metagene`` and ``tis_advisor``.
This is reported in the startup stage list.
