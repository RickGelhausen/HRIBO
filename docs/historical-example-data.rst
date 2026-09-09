Historical public example data
==============================

The former ``HRIBO_ReadTheDocs`` repository used the public datasets below for
its basic and extended examples.  Their accession mapping is retained here so
that deleting that repository does not discard useful provenance.

.. important::

   These are **candidate validation datasets**, not certified HRIBO 2.0
   fixtures.  The old examples were run with HRIBO 1.x, and their commands,
   screenshots, runtimes, and numerical results have not been carried forward.
   Confirm the current archive metadata and reference files, record checksums,
   and complete :doc:`real-data-validation` before publishing a 2.0 example.

Pseudomonas aeruginosa PAO1
---------------------------

The compact historical example selected two runs from `NCBI BioProject
PRJNA379630 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA379630>`_ and used
the PAO1 reference assembly `GCF_000006765.1
<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006765.1/>`_ (primary
chromosome accession ``NC_002516.2``).  The originating multi-omics study is
described by :cite:t:`grady2017multiomics`.

.. list-table:: Retained PAO1 run mapping
   :header-rows: 1
   :widths: 25 20 20 35

   * - HRIBO method
     - Condition
     - Replicate
     - SRA run
   * - ``RIBO``
     - ``GLY``
     - ``1``
     - ``SRR5356907``
   * - ``RNA``
     - ``GLY``
     - ``1``
     - ``SRR5356908``

This subset has only one condition and one replicate.  It can exercise input
validation, preprocessing, mapping, tracks, read counts, metagene profiling,
and prediction, but it is not a valid differential-expression design.  Do not
use its historical condition label as a substitute for checking the run
metadata.

Salmonella Typhimurium 14028S
-----------------------------

The extended historical example selected LB-grown wild-type and ``csrA``
libraries from `NCBI BioProject PRJNA421559
<https://www.ncbi.nlm.nih.gov/bioproject/PRJNA421559>`_.  It used the
*Salmonella enterica* subsp. *enterica* serovar Typhimurium strain 14028S
assembly `GCF_000022165.1
<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000022165.1/>`_ (primary
chromosome accession ``NC_016856.1``).  The associated experiment is described
by :cite:t:`potts2019csra`.

.. list-table:: Retained Salmonella run mapping
   :header-rows: 1
   :widths: 25 20 20 35

   * - HRIBO method
     - Condition
     - Replicate
     - SRA run
   * - ``RIBO``
     - ``WT``
     - ``1``
     - ``SRR6359966``
   * - ``RIBO``
     - ``WT``
     - ``2``
     - ``SRR6359967``
   * - ``RIBO``
     - ``csrA``
     - ``1``
     - ``SRR6359970``
   * - ``RIBO``
     - ``csrA``
     - ``2``
     - ``SRR6359971``
   * - ``RNA``
     - ``WT``
     - ``1``
     - ``SRR6359974``
   * - ``RNA``
     - ``WT``
     - ``2``
     - ``SRR6359975``
   * - ``RNA``
     - ``csrA``
     - ``1``
     - ``SRR6359978``
   * - ``RNA``
     - ``csrA``
     - ``2``
     - ``SRR6359979``

The last row corrects a legacy sample-sheet typo that pointed replicate 2 at
the replicate-1 FASTQ.  One archive ambiguity remains: sample ``ST_10``
(``SRR6359975``) is currently titled ``RNA WT 1 LB`` in `GEO series GSE107834
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107834>`_, duplicating
the label for ``ST_9`` even though the series design describes biological
replicates.  The table preserves the legacy assignment of ``ST_10`` as
replicate 2; confirm that interpretation against the current archive metadata
before analysis.

If that assignment is confirmed, this subset has two matched RIBO/RNA
replicates per condition and therefore the structure needed to exercise the
current differential engines, subject to the normal read-depth and design
checks.

Preparing a future validation run
---------------------------------

Use the current `NCBI SRA download guidance
<https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/>`_ instead of the former
FTP and Conda commands.  At minimum, record:

* retrieval date, run accessions, archive metadata, and downloaded-file
  checksums;
* exact assembly and annotation accessions, versions, checksums, and sequence
  identifiers;
* the final ``samples.tsv`` and effective ``config.yaml``;
* every conversion, compression, or read-layout decision; and
* HRIBO baseline/candidate commits, commands, environments, logs, and the
  semantic comparison report.

HRIBO requires gzip-compressed, structurally valid FASTQ input, but the
retrieval tool and temporary archive format are outside the workflow contract.
Inspect whether each selected SRA run is single- or paired-end before filling
``fastqFile2``.  The maintained execution pattern is in
:doc:`getting-started`; do not reconstruct commands from the archived 1.x
pages.
