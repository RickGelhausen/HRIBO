Sample sheet
============

``biologySettings.samples`` points to a tab-separated sample sheet.  Use one
row per sequencing library and preserve the exact header names:

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

Columns
-------

``method``
   One of ``RIBO``, ``TIS``, ``TTS``, ``RNA``, ``RNATIS``, or ``RNATTS``.
   ``RIBO`` is ordinary ribosome profiling and ``RNA`` is its matched
   transcriptome control.  ``TIS`` and ``TTS`` identify initiation- and
   termination-enriched profiling libraries; ``RNATIS`` and ``RNATTS`` are
   their RNA controls.

``condition``
   An alphanumeric biological condition name.  Hyphens are not allowed because
   ``-`` separates the two sides of a differential contrast.

``replicate``
   A canonical positive base-ten integer written as text: ``1``, ``2``, and so
   on.  Zero, signs, decimal notation, surrounding whitespace, and leading
   zeros such as ``01`` are rejected.  The tuple ``method``, ``condition``,
   ``replicate`` must be unique.

``fastqFile``
   The gzip-compressed FASTQ for a single-end library, or read 1 for a paired
   library.  Compression is detected from the file content, not its suffix;
   plain-text FASTQ is not supported.  Relative paths are interpreted from the
   analysis directory.

``fastqFile2``
   Read 2 for a paired-end library.  Leave the field empty for single-end data.
   Single- and paired-end libraries can coexist.  Paired reads are adapter-
   trimmed and merged with PEAR; the assembled read continues through mapping.

Design requirements
-------------------

The requested stages determine which design constraints apply:

* ``predictions`` and ``overview`` require at least one ``RIBO`` library.
* ``metagene`` and ``tis_advisor`` accept ``RIBO``, ``TIS``, or ``TTS``.
* ``differential_expression`` requires at least two conditions with matched
  ``RIBO`` and ``RNA`` libraries, at least two replicates per side, and equal
  selected RIBO/RNA replicate counts.
* PCA requires at least two usable libraries when the ``pca`` stage is selected.

The schema is ``workflow/schemas/samples.schema.yaml``.  HRIBO additionally
checks uniqueness, paths, the complete gzip stream and CRC, every four-line
FASTQ record, and the stage-specific design before running.  This preflight is
linear in the compressed inputs: it deliberately reads each distinct FASTQ to
EOF so corruption or malformed records near the end cannot pass unnoticed.
