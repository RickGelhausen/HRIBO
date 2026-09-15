
rule prepareRiborexInput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}",
        script=str(SCRIPTS / "prepare_diffex_input.py")
    output:
        ribo="diffex_input/riborex/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/riborex/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/riborex/{contrast}_condition_vector.csv"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        """
        python3 {input.script:q} -r {input.rawreads:q} -c {wildcards.contrast:q} -t riborex -o diffex_input/riborex/
        """


rule riborex:
    input:
        ribo="diffex_input/riborex/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/riborex/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/riborex/{contrast}_condition_vector.csv",
        script=str(SCRIPTS / "riborex.R")
    output:
        table="riborex/{contrast}_deseq2.csv"
    conda:
        "../envs/riborex.yaml"
    threads: 1
    shell:
        """
        Rscript {input.script:q} -r {input.ribo:q} -m {input.rna:q} -c {input.cv:q} -x {output.table:q}
        """

rule riborexxlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        riborex_out="riborex/{contrast}_deseq2.csv",
        script=str(SCRIPTS / "generate_excel_riborex.py"),
        script_deps=[
            str(SCRIPTS / "excel_utils.py"),
            str(SCRIPTS / "gff_utils.py"),
            str(SCRIPTS / "lib" / "__init__.py"),
            str(SCRIPTS / "lib" / "misc.py"),
        ]
    output:
        xlsx_sorted="riborex/{contrast}_sorted.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        padj_cutoff=config["differentialExpressionSettings"]["padjCutoff"],
        log2fc_cutoff=config["differentialExpressionSettings"]["log2fcCutoff"]
    shell:
        """
        python3 {input.script:q} -a {input.annotation:q} -g {input.genome:q} -i {input.riborex_out:q} -o {output.xlsx_sorted:q} --padj_cutoff {params.padj_cutoff:q} --log2fc_cutoff {params.log2fc_cutoff:q}
        """

rule poolriborex:
    input:
        riborex=expand("riborex/{contr}_sorted.xlsx", contr=CONTRASTS),
        script=str(SCRIPTS / "merge_differential_expression.py")
    output:
        "riborex/riborex_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 {input.script:q} {input.riborex:q} -o {output:q} -t riborex
        """
