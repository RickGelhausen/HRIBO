
rule longestTranscript:
    input:
        rules.retrieveAnnotation.output
    output:
        "xtail/longest_protein_coding_transcripts.gtf"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell:
        "mkdir -p xtail; HRIBO/scripts/longest_orf_transcript.py -a {input} -o {output}"

rule sizeFactors:
    input:
        rules.longestTranscript.output,
        expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        "normalization/sfactors.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p normalization; HRIBO/scripts/generate_size_factors.R -t HRIBO/samples.tsv -b maplink/ -a {input[0]} -s {output};")

rule cdsNormalizedCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="tracks/updated_annotation.gff",
    output:
        raw="normalization/raw_reads.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p normalization; HRIBO/scripts/generate_raw_counts.R -b maplink/ -a {input.annotation} -t HRIBO/samples.tsv -r {output.raw};")

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
        rawreads="normalization/raw_reads.csv",
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
        annotation=rules.retrieveAnnotation.output,
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

rule riborex:
    input:
        rawreads="normalization/raw_reads.csv",
        contrastfile="contrasts/{contrast}"
    output:
        tabledeseq2="riborex/{contrast}_deseq2.csv",
    conda:
        "../envs/riborex.yaml"
    threads: 1
    shell: ("mkdir -p riborex; HRIBO/scripts/riborex.R -c {input.contrastfile} -t HRIBO/samples.tsv -r {input.rawreads} -x {output.tabledeseq2};")

rule riborexresults:
    input:
        tabledeseq2="riborex/{contrast}_deseq2.csv"
    output:
        tablesorted="riborex/{contrast}_sorted.csv",
        tablesignificant=report("riborex/{contrast}_significant.csv", caption="../report/riborex.rst", category="Regulation")
    threads: 1
    shell: ("mkdir -p riborex; (head -n 2 {input.tabledeseq2} && tail -n +3 {input.tabledeseq2} | sort -r -n -t',' -k 7) > {output.tablesorted};  awk -F ',' 'NR==1; (NR>1) && ($7 < 0.05 )' {output.tablesorted} > {output.tablesignificant};")

rule riborexxlsx:
    input:
        annotation=rules.retrieveAnnotation.output,
        genome=rules.retrieveGenome.output,
        riborex_sorted="riborex/{contrast}_sorted.csv",
        riborex_signif="riborex/{contrast}_significant.csv"
    output:
        xlsx_sorted="riborex/{contrast}_sorted.xlsx",
        xlsx_signif="riborex/{contrast}_significant.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/differential_expression_xlsx.py -a {input.annotation} -g {input.genome} --tool riborex -i {input.riborex_sorted} -o {output.xlsx_sorted}
        python3 HRIBO/scripts/differential_expression_xlsx.py -a {input.annotation} -g {input.genome} --tool riborex -i {input.riborex_signif} -o {output.xlsx_signif}
        """

cur_contrast=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(samples["condition"].unique(),2)))] for item in sublist]
rule poolriborex:
    input:
        riborex=expand("riborex/{contr}_sorted.xlsx", contr=cur_contrast)
    output:
        "riborex/riborex_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.riborex} -o {output} -t riborex
        """

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
