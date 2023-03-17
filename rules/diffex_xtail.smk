
rule xtail:
    input:
        ribo="diffex_input/xtail/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/xtail/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/xtail/{contrast}_condition_vector.csv"
    output:
        table="xtail/{contrast}.csv",
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf",
    conda:
        "../envs/xtail.yaml"
    threads: 10
    shell:
        """
        mkdir -p xtail;
        HRIBO/scripts/xtail.R -r {input.ribo} -m {input.rna} -c {input.cv} -x {output.table} -f {output.fcplot} -p {output.rplot};
        """

rule xtailxlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        xtail_out="xtail/{contrast}.csv",
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
        python3 HRIBO/scripts/generate_excel_xtail.py -a {input.annotation} -g {input.genome} -i {input.xtail_out} -o {output.xlsx_sorted} --padj_cutoff {params.padj_cutoff} --log2fc_cutoff {params.log2fc_cutoff}
        """

rule poolxtail:
    input:
        xtail=expand("xtail/{contr}_sorted.xlsx", contr=CONTRASTS)
    output:
        "xtail/xtail_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.xtail} -o {output} -t xtail
        """
