The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`__ and `Cutadapt <http://cutadapt.readthedocs.io>`__, afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>`__. `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ is used after each of these preprocessing steps.
The reads are then mapped with `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`__. New ORFs not contained in the provided annotation are detected with `RiboTish <https://github.com/zhpn1024/ribotish>`__
and a custom modified version of `Reparation <https://github.com/RickGelhausen/REPARATION_blast>`__.

Results:

Quality control summary by MultiQC (containing FastQC, trimming, featurecount), see Section `Quality control`_ .
ee Section `Global tracks`_ using the full read length for coverage file computation.
See Section `Centered tracks`_ using only the nucleotides in the center of the read for coverage file computation.
See Section `5' single nucleotide mapping tracks`_ using only the nucleotide at the five prime end of the read for coverage file computation.
See Section `3' single nucleotide mapping tracks`_ using only the nucleotide at the three prime end of the read for coverage file computation.
Please inspect the included folders of the archive, specifically the genomes, annotation and track folders for genome browser visualisation.
A excel table containing identified Open reading frames can be found in Section `Novel ORFs`_ , together with genome browser tracks
showing the identified tracks. For a overview load combined.gff showing in which conditions and by which method the ORFs were identified.
ORFs are named according to a genome:start-stop:strand scheme, to make them trackable over multiple experiments.
A good strategy is to compare the combined.gff track with the established annotation and find which predicted ORFs are not present.
Loci of interest can be further investigated by loading the TPM normalized, strand specific big wig files, which show how strong the detected signal was for a specific sample.
Additionally the genome is screened for potential start(ATG,GTG,TTG), stop(TAG,TGA,TAA) codons and ribosome binding sites(AGGAG), which can be loaded via the potentialRibosomeBindingSite.gff, potentialStartCodon.gff and potentialStopCodon.gff. 
If these sequences are different for your organism please contact us for modifing them.

