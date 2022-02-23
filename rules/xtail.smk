rule contrastInput:
    output:
        "contrasts/{contrast}"
    run:
        if not os.path.exists("contrasts"):
            os.makedirs("contrasts")
        for f in getContrast(wildcards):
            print(f)
            open((f), 'a').close()

rule xtail:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}"
    output:
        table=report("xtail/{contrast}.csv", caption="../report/xtail_table.rst", category="Regulation"),
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf",
        tablesorted="xtail/{contrast}_sorted.csv",
        tablesignificant="xtail/{contrast}_significant.csv"
    conda:
        "../envs/xtail.yaml"
    threads: 1
    shell:
        """
        mkdir -p xtail;
        HRIBO/scripts/xtail_classic.R -c {input.contrastfile} -t HRIBO/samples.tsv -r {input.rawreads} -x {output.table} -f {output.fcplot} -p {output.rplot};
        (head -n 2 {output.table} && tail -n +3 {output.table} | sort -r -n -t',' -k 10) > {output.tablesorted};  awk -F ',' 'NR==1; (NR>1) && ($10 < 0.05 )' {output.tablesorted} > {output.tablesignificant};
        """

rule xtailxlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        xtail_sorted="xtail/{contrast}_sorted.csv",
        xtail_signif="xtail/{contrast}_significant.csv"
    output:
        xlsx_sorted="xtail/{contrast}_sorted.xlsx",
        xlsx_signif="xtail/{contrast}_significant.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/differential_expression_xlsx.py -a {input.annotation} -g {input.genome} --tool xtail -i {input.xtail_sorted} -o {output.xlsx_sorted}
        python3 HRIBO/scripts/differential_expression_xlsx.py -a {input.annotation} -g {input.genome} --tool xtail -i {input.xtail_signif} -o {output.xlsx_signif}
        """

cur_contrast=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(sorted(samples["condition"].unique()),2)))] for item in sublist]

rule poolxtail:
    input:
        xtail=expand("xtail/{contr}_sorted.csv", contr=cur_contrast)
    output:
        "xtail/xtail_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.xtail} -o {output} -t xtail
        """
