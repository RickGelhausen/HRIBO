Prepare the sample sheet
========================

``biologySettings.samples`` points to a tab-separated file with one row per
sequencing library.  Copy ``config/samples.tsv`` into the analysis directory
and preserve these five header names exactly:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2

Single-end example
------------------

Leave ``fastqFile2`` empty for single-end libraries.  Quoted empty fields, as
shown here, are accepted:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2
   RIBO	Control	1	fastq/ribo-control-1.fastq.gz	""
   RIBO	Control	2	fastq/ribo-control-2.fastq.gz	""
   RIBO	Treated	1	fastq/ribo-treated-1.fastq.gz	""
   RIBO	Treated	2	fastq/ribo-treated-2.fastq.gz	""
   RNA	Control	1	fastq/rna-control-1.fastq.gz	""
   RNA	Control	2	fastq/rna-control-2.fastq.gz	""
   RNA	Treated	1	fastq/rna-treated-1.fastq.gz	""
   RNA	Treated	2	fastq/rna-treated-2.fastq.gz	""

Paired-end example
------------------

Put read 1 in ``fastqFile`` and read 2 in ``fastqFile2``:

.. code-block:: text

   method	condition	replicate	fastqFile	fastqFile2
   RIBO	Control	1	fastq/ribo-control-1_R1.fastq.gz	fastq/ribo-control-1_R2.fastq.gz
   RNA	Control	1	fastq/rna-control-1_R1.fastq.gz	fastq/rna-control-1_R2.fastq.gz

Single- and paired-end libraries may occur in the same sheet.  Paired reads
are adapter-trimmed and merged before mapping, so review the four paired-end
adapter settings in :doc:`configuration` and set the ones that apply to the
library preparation.

Column reference
----------------

``method``
   The library type.  Accepted values are ``RIBO``, ``TIS``, ``TTS``, ``RNA``,
   ``RNATIS``, and ``RNATTS``.  ``RIBO`` is standard ribosome profiling and
   ``RNA`` is its matched RNA-seq control.  ``TIS`` and ``TTS`` are
   initiation- and termination-enriched profiling libraries; ``RNATIS`` and
   ``RNATTS`` are their RNA controls.

``condition``
   An alphanumeric biological condition such as ``Control``, ``Treated``, or
   ``Mutant1``.  Spaces, underscores, and hyphens are not accepted.  A hyphen
   is reserved for separating the two conditions in a contrast.

``replicate``
   A positive integer written as ``1``, ``2``, and so on.  Do not use zero,
   signs, decimals, whitespace, or leading zeros.  Each
   ``method``/``condition``/``replicate`` combination must be unique.

``fastqFile``
   The FASTQ for a single-end library or read 1 for paired-end data.  It must
   be gzip-compressed.  Relative paths are resolved from the analysis
   directory.

``fastqFile2``
   Read 2 for paired-end data.  Leave the cell empty for a single-end library.
   Do not write explanatory text such as ``NA`` or ``single-end`` in the cell.

How names appear in the results
-------------------------------

HRIBO names a library ``<method>-<condition>-<replicate>``.  For example, the
first row above becomes ``RIBO-Control-1`` in BAM, BigWig, report, and workbook
names.  Choose short, meaningful condition labels because they are repeated
throughout the result directory.

Design requirements
-------------------

Most stages accept any non-empty sample sheet, but some analyses need a
specific design:

* ``predictions`` and ``overview`` need at least one ``RIBO`` library.
* ``metagene`` and ``tis_advisor`` use ``RIBO``, ``TIS``, or ``TTS``
  libraries.
* ``differential_expression`` needs at least two conditions with matched
  ``RIBO`` and ``RNA`` libraries, at least two biological replicates on each
  side, and equal selected RIBO/RNA replicate counts.
* ``pca`` needs at least two usable libraries.

Run a dry-run after editing the sheet.  HRIBO checks the column names, library
labels, duplicate rows, referenced paths, gzip integrity, and FASTQ structure
before starting the analysis.
