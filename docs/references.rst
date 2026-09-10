References
==========

When publishing an analysis performed with HRIBO, cite the HRIBO workflow
paper :cite:p:`gelhausen2021hribo`.  Record the HRIBO version or commit and the
configuration alongside the citation so the analysis can be reproduced.

For background on ribosome profiling and the interpretation of protected
fragments, see :cite:p:`ingolia2014ribosome`.

HRIBO integrates two prokaryotic ORF predictors.  Cite the DeepRibo paper
:cite:p:`clauwaert2019deepribo` when reporting DeepRibo calls and the
REPARATION paper :cite:p:`ndah2017reparation` when reporting REPARATION calls.

The workflow is implemented with Snakemake :cite:p:`koster2012snakemake` and
uses packages distributed through Bioconda :cite:p:`gruning2018bioconda`.
Relevant processing and output methods include Cutadapt
:cite:p:`martin2011cutadapt`, segemehl :cite:p:`otto2014segemehl`, and the
BigWig format :cite:p:`kent2010bigwig`.  When reporting differential
translation, cite each selected engine as appropriate: xTail
:cite:p:`xiao2016xtail`, RiboRex :cite:p:`li2017riborex`, and/or deltaTE
:cite:p:`chothani2019deltate`.  RiboRex and deltaTE also use the DESeq2
framework :cite:p:`love2014deseq2`.

Bibliography
------------

.. bibliography::
   :style: plain
