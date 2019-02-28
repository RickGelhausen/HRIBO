The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`__ and `Cutadapt <http://cutadapt.readthedocs.io>`__, afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>`__. `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ is used after each of these preprocessing steps.
The reads are then mapped with `Segemehl <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`__. New ORFs not contained in the provided annotation are detected with `RiboTish <https://github.com/zhpn1024/ribotish>`__
and a custom modified version of `Reparation <https://github.com/RickGelhausen/REPARATION_blast>`__.
The mapped reads are visualized as wig files. Diffential expression of ORFs is evaluated with `Xtail <https://github.com/xryanglab/xtail>`_ and `Riborex <https://github.com/smithlabcode/riborex>`__.
We used the genome for `Organism <https://www.ensembl.org/index.html>`__  with the following annotation `Gtf <https://www.ensembl.org/index.html>`__ .

