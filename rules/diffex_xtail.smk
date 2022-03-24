
rule xtail:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}"
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
        HRIBO/scripts/xtail_classic.R -c {input.contrastfile} -t HRIBO/samples.tsv -r {input.rawreads} -x {output.table} -f {output.fcplot} -p {output.rplot};
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
    shell:
        """
        python3 HRIBO/scripts/generate_excel_xtail.py -a {input.annotation} -g {input.genome} -i {input.xtail_out} -o {output.xlsx_sorted}
        """

cur_contrast=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(sorted(samples["condition"].unique(), key=lambda s: s.lower()),2)))]  for item in sublist]

rule poolxtail:
    input:
        xtail=expand("xtail/{contr}_sorted.xlsx", contr=cur_contrast)
    output:
        "xtail/xtail_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.xtail} -o {output} -t xtail
        """
