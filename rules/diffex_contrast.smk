rule contrastInput:
    output:
        "contrasts/{contrast}"
    run:
        if not os.path.exists("contrasts"):
            os.makedirs("contrasts")
        for f in getContrast(wildcards):
            print(f)
            open((f), 'a').close()

rule prepareRiborexInput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}"
    output:
        ribo="diffex_input/riborex/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/riborex/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/riborex/{contrast}_condition_vector.csv"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        """
        mkdir -p diffex_input/riborex/;
        python3 HRIBO/scripts/prepare_diffex_input.py -r {input.rawreads} -c {wildcards.contrast} -o diffex_input/riborex/
        """


rule prepareXtailInput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}"
    output:
        ribo="diffex_input/xtail/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/xtail/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/xtail/{contrast}_condition_vector.csv"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        """
        mkdir -p diffex_input/xtail/;
        python3 HRIBO/scripts/prepare_diffex_input.py -r {input.rawreads} -c {wildcards.contrast} -o diffex_input/xtail/
        """
