
rule xtail:
    input:
        ribo="diffex_input/xtail/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/xtail/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/xtail/{contrast}_condition_vector.csv",
        script=str(SCRIPTS / "xtail.R")
    output:
        table="xtail/{contrast}.csv",
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf",
    conda:
        "../envs/xtail.yaml"
    threads: 10
    params:
        bins=config["differentialExpressionSettings"].get("xtailBins", 10000),
        min_mean_count=config["differentialExpressionSettings"].get(
            "xtailMinMeanCount", 1
        ),
    shell:
        """
        Rscript {input.script:q} \
            -r {input.ribo:q} \
            -m {input.rna:q} \
            -c {input.cv:q} \
            -x {output.table:q} \
            -f {output.fcplot:q} \
            -p {output.rplot:q} \
            --threads {threads} \
            --bins {params.bins:q} \
            --min_mean_count {params.min_mean_count:q}
        """

rule xtailxlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        xtail_out="xtail/{contrast}.csv",
        script=str(SCRIPTS / "generate_excel_xtail.py"),
        script_deps=[
            str(SCRIPTS / "excel_utils.py"),
            str(SCRIPTS / "gff_utils.py"),
        ]
    output:
        xlsx_sorted="xtail/{contrast}_sorted.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        padj_cutoff=config["differentialExpressionSettings"]["padjCutoff"],
        log2fc_cutoff=config["differentialExpressionSettings"]["log2fcCutoff"]
    shell:
        """
        python3 {input.script:q} -a {input.annotation:q} -g {input.genome:q} -i {input.xtail_out:q} -o {output.xlsx_sorted:q} --padj_cutoff {params.padj_cutoff:q} --log2fc_cutoff {params.log2fc_cutoff:q}
        """

rule poolxtail:
    input:
        xtail=expand("xtail/{contr}_sorted.xlsx", contr=CONTRASTS),
        script=str(SCRIPTS / "merge_differential_expression.py")
    output:
        "xtail/xtail_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 {input.script:q} {input.xtail:q} -o {output:q} -t xtail
        """
