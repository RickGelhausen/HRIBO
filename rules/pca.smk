rule preparePCAinput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        samples=config["biologySettings"]["samples"],
    output:
        rawreads="pca/raw_reads.csv",
        meta="pca/meta.csv",
    threads: 1
    conda:
        "../envs/pytools.yaml"
    shell:
        """
        mkdir -p pca;
        sed -e '1s/-/_/g' {input.rawreads} > {output.rawreads};
        HRIBO/scripts/preparePCAinput.py -s {input.samples} -o {output.meta};
        """

rule runDeseqPreprocessing:
    input:
        rawreads="pca/raw_reads.csv",
        meta="pca/meta.csv"
    output:
        distr="pca/raw_count_distributions.pdf",
        mv="pca/mean_vs_variance.pdf",
        norm="pca/normalized_counts.tsv",
        rld="pca/rld.tsv",
        pvar="pca/variance_percentages.tsv"
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        mkdir -p pca;
        HRIBO/scripts/analyse_variance.R -r {input.rawreads} -m {input.meta} -o pca/;
        """

rule plotPCA:
    input:
        rld="pca/rld.tsv",
        pvar="pca/variance_percentages.tsv"
    output:
        plot="pca/PCA_3D.html",
        plot2="diffex_QC.html"
    threads: 1
    conda:
        "../envs/plotly.yaml"
    shell:
        """
        mkdir -p pca;
        HRIBO/scripts/plot_PCA.py -r {input.rld} -p {input.pvar} -o pca/;
        """