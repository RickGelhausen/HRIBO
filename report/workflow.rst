The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`__ and `Cutadapt <http://cutadapt.readthedocs.io>`__, afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>`__. `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ is used after each of these preprocessing steps.
The reads are then mapped with `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`__. New ORFs not contained in the provided annotation are detected with `RiboTish <https://github.com/zhpn1024/ribotish>`__
and a custom modified version of `Reparation <https://github.com/RickGelhausen/REPARATION_blast>`__.
Diffential expression of ORFs is evaluated with `Xtail <https://github.com/xryanglab/xtail>`_ and `Riborex <https://github.com/smithlabcode/riborex>`__.

Results:

Quality control summary by MultiQC (containing FastQC, trimming, featurecount), see Section `Quality control`_ .

Mapped reads, see Section `Mapped tracks`_ .

Open reading frames detected by Ribotish and Reparation for different conditions, see gff files in Section `Novel ORFs`_ .

Differential Expression for each combination of conditions detected by Xtail and Riborex, see Section `Regulation`_ .

Please inspect the included folders of the archive, specifically the genomes, annotation and track folders for genome browser visualisation. 

The combined.gff track summarizes the detected ORF over all used conditions and methods (Ribotish, Reparation).

ORFs are named according to a genome:start-stop:strand scheme, to make them trackable over multiple experiments.

A good strategy is to compare the combined.gff track with the established annotation and find which predicted ORFs are not present. 

Loci of interest can be further investigated by loading the TPM normalized, strand specific big wig files, which show how strong the detected signal was for a specific sample.

Additionally the genome is screened for potential start(ATG,GTG,TTG), stop(TAG,TGA,TAA) codons and ribosome binding sites(AAGG), which can be loaded via the potentialRibosomeBindingSite.gff, potentialStartCodon.gff and potentialStopCodon.gff.
