rule preparePCAinput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        samples=config["biologySettings"]["samples"],
        script=str(SCRIPTS / "preparePCAinput.py"),
    output:
        rawreads="pca/raw_reads.csv",
        meta="pca/meta.csv",
    threads: 1
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        sed -e '1s/-/_/g' {input.rawreads:q} > {output.rawreads:q}
        python3 {input.script:q} -s {input.samples:q} -o {output.meta:q}
        """

rule runDeseqPreprocessing:
    input:
        rawreads="pca/raw_reads.csv",
        meta="pca/meta.csv",
        script=str(SCRIPTS / "analyse_variance.R")
    output:
        distr="pca/raw_count_distributions.pdf",
        mv="pca/mean_vs_variance.pdf",
        norm="pca/normalized_counts.tsv",
        rld="pca/rld.tsv",
        pvar="pca/variance_percentages.tsv",
        cor="pca/rld_cor.tsv"
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        Rscript {input.script:q} -r {input.rawreads:q} -m {input.meta:q} -o pca/
        """
        
rule plotPCA:
    input:
        rld="pca/rld.tsv",
        pvar="pca/variance_percentages.tsv",
        cor="pca/rld_cor.tsv",
        script=str(SCRIPTS / "plot_PCA.py")
    output:
        plot="pca/PCA_3D.html",
        plot2="pca/diffex_QC.html"
    threads: 1
    conda:
        "../envs/plotly.yaml"
    shell:
        """
        python3 {input.script:q} -r {input.rld:q} -p {input.pvar:q} -c {input.cor:q} -o pca/
        """
