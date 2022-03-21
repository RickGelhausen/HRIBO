
rule riborex:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}"
    output:
        table="riborex/{contrast}_deseq2.csv"
    conda:
        "../envs/riborex.yaml"
    threads: 1
    shell:
        """
        mkdir -p riborex;
        HRIBO/scripts/riborex.R -c {input.contrastfile} -t HRIBO/samples.tsv -r {input.rawreads} -x {output.table};
        """

rule riborexxlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        riborex_out="riborex/{contrast}_deseq2.csv"
    output:
        xlsx_sorted="riborex/{contrast}_sorted.xlsx",
        xlsx_signif="riborex/{contrast}_significant.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/generate_excel_xtail.py -a {input.annotation} -g {input.genome} -i {input.riborex_out} -o {output.xlsx_sorted}
        """

cur_contrast=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(samples["condition"].unique(),2)))] for item in sublist]

rule poolriborex:
    input:
        riborex=expand("riborex/{contr}_sorted.csv", contr=cur_contrast)
    output:
        "riborex/riborex_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.riborex} -o {output} -t riborex
        """