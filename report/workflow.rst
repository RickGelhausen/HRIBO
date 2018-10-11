The data is processed in the following steps.
Adapter removal is performed with `TrimGalore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`_ and `Cutadapt <http://cutadapt.readthedocs.io>`_, afterwards the reads mapping to rRNA genes are
removed with `SortMeRNA <http://bioinfo.lifl.fr/RNA/sortmerna/>`_. `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ is used after each of these preprocessing steps. 
The reads are then mapped with `STAR <https://github.com/alexdobin/STAR>`_. New ORFs not contained in the provided annotation are detected with `RiboTish <https://github.com/zhpn1024/ribotish>`_
and a custom modified version of `Reparation <https://github.com/RickGelhausen/REPARATION_blast>`_.
The mapped reads are visualized as wig files. Diffential expression of ORFs is evaluated with `Xtail <https://github.com/xryanglab/xtail>`_ and `Riborex <https://github.com/smithlabcode/riborex>`_.
