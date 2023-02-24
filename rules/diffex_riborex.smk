
rule riborex:
    input:
        ribo="diffex_input/riborex/{contrast}_ribo_readcount_table.tsv",
        rna="diffex_input/riborex/{contrast}_rna_readcount_table.tsv",
        cv="diffex_input/riborex/{contrast}_condition_vector.csv"
    output:
        table="riborex/{contrast}_deseq2.csv"
    conda:
        "../envs/riborex.yaml"
    threads: 1
    shell:
        """
        mkdir -p riborex;
        HRIBO/scripts/riborex.R -r {input.ribo} -m {input.rna} -c {input.cv} -x {output.table};
        """

rule riborexxlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        riborex_out="riborex/{contrast}_deseq2.csv"
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
        python3 HRIBO/scripts/generate_excel_riborex.py -a {input.annotation} -g {input.genome} -i {input.riborex_out} -o {output.xlsx_sorted} --padj_cutoff {params.padj_cutoff} --log2fc_cutoff {params.log2fc_cutoff}
        """

rule poolriborex:
    input:
        riborex=expand("riborex/{contr}_sorted.xlsx", contr=CONTRASTS)
    output:
        "riborex/riborex_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.riborex} -o {output} -t riborex
        """