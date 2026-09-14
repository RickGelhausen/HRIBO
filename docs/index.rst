HRIBO |release|
================

HRIBO is a Snakemake workflow for analysing bacterial ribosome-profiling
experiments.  Starting from compressed sequencing reads, a reference genome,
and its annotation, HRIBO can produce quality-control reports, mapped reads,
genome-browser tracks, metagene profiles, ORF predictions, differential
translation results, and consolidated result tables.

If this is your first analysis, follow :doc:`getting-started` from top to
bottom.  It covers installation, project setup, a dry-run, execution, and the
first results to inspect.

What goes in
------------

HRIBO needs:

* a reference genome in FASTA format;
* a matching annotation in GFF3 or GTF format;
* one gzip-compressed FASTQ file per single-end library, or two per paired-end
  library;
* a tab-separated sample sheet describing every library; and
* a YAML configuration file selecting the analyses to run.

See :doc:`samples` and :doc:`configuration` for copyable examples.

What comes out
--------------

The exact result set depends on the selected :doc:`stages`.  The main result
types are:

.. list-table:: Main HRIBO results
   :header-rows: 1
   :widths: 28 32 40

   * - Result
     - Main location
     - What it is used for
   * - Quality control
     - ``qc/multi/multiqc_report.html``
     - Review read quality, trimming, mapping, and rRNA/tRNA depletion.
   * - Alignments and coverage
     - ``maplink/`` and ``*tracks/``
     - Inspect reads and strand-aware coverage in a genome browser.
   * - Counts and abundance
     - ``auxiliary/*.xlsx``
     - Compare feature counts, RPKM values, and direct TE ratios.
   * - Metagene and TIS reports
     - ``metageneprofiling/`` and ``tis_advice/``
     - Assess read lengths, start/stop profiles, and P-site offsets.
   * - ORF predictions
     - ``auxiliary/predictions_*.xlsx``
     - Review REPARATION and optional DeepRibo candidates.
   * - Differential translation and condition overview
     - ``diffex_summary/``, ``xtail/``, ``riborex/``, and ``deltate/``
     - Scan condition-specific detection and RNA, RIBO, or TE changes, then
       inspect detailed contrast results.
   * - Combined overview
     - ``auxiliary/overview.xlsx``
     - Explore annotation, abundance, prediction, and differential evidence in
       one workbook.

Start the result review with :doc:`outputs`, which explains what each file
answers and which files are primary results rather than supporting data.

.. toctree::
   :maxdepth: 1
   :caption: Run HRIBO

   getting-started
   samples
   configuration
   stages
   tutorials/minimal
   tutorials/full

.. toctree::
   :maxdepth: 1
   :caption: Understand the results

   outputs
   table-reference
   differential-summary
   metagene-profiling
   tis-advisor

.. toctree::
   :maxdepth: 1
   :caption: Help and reference

   troubleshooting
   migration-1.8-to-2.0
   references

Getting help
------------

See :doc:`troubleshooting` for common input, environment, and cluster issues.
To report a reproducible problem, use the `HRIBO issue tracker
<https://github.com/RickGelhausen/HRIBO/issues>`_ and include the HRIBO version,
configuration, first failing rule, and relevant log file.

HRIBO is distributed under the GNU General Public License version 3.  If HRIBO
contributes to published work, follow the citation guidance in
:doc:`references`.
