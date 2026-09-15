rule contrastInput:
    output:
        "contrasts/{contrast}"
    shell:
        "touch {output:q}"

rule prepareXtailInput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}",
        script=str(SCRIPTS / "prepare_diffex_input.py")
    output:
        ribo="diffex_input/xtail/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/xtail/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/xtail/{contrast}_condition_vector.csv"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        """
        python3 {input.script:q} -r {input.rawreads:q} -c {wildcards.contrast:q} -t xtail -o diffex_input/xtail/
        """
