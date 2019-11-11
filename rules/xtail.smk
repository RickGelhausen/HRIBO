rule longestTranscript:
    input:
        rules.retrieveAnnotation.output
    output:
        "xtail/longest_protein_coding_transcripts.gtf"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell:
        "mkdir -p xtail; SPtools/scripts/longest_orf_transcript.py -a {input} -o {output}"

rule sizeFactors:
    input:
        rules.longestTranscript.output,
        expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        "normalization/sfactors.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p normalization; SPtools/scripts/generate_size_factors.R -t SPtools/samples.tsv -b maplink/ -a {input[0]} -s {output};")

rule cdsNormalizedCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="xtail/newAnnotation.gff",
    output:
        raw="normalization/raw_reads.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p normalization; SPtools/scripts/generate_raw_counts.R -b maplink/ -a {input.annotation} -t SPtools/samples.tsv -r {output.raw};")

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
        SPtools/scripts/xtail_classic.R -c {input.contrastfile} -t SPtools/samples.tsv -r {input.rawreads} -x {output.table} -f {output.fcplot} -p {output.rplot};
        head -n 2 {output.table} && tail -n +3 {output.table} | sort -r -n -t',' -k 10 > {output.tablesorted};  awk -F ',' 'NR==1; (NR>1) && ($10 < 0.05 )' {output.tablesorted} > {output.tablesignificant};
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
    shell: ("mkdir -p riborex; SPtools/scripts/riborex.R -c {input.contrastfile} -t SPtools/samples.tsv -r {input.rawreads} -x {output.tabledeseq2};")

rule riborexresults:
    input:
        tabledeseq2="riborex/{contrast}_deseq2.csv"
    output:
        tablesorted="riborex/{contrast}_sorted.csv",
        tablesignificant=report("riborex/{contrast}_significant.csv", caption="../report/riborex.rst", category="Regulation")
    threads: 1
    shell: ("mkdir -p riborex; (head -n 2 {input.tabledeseq2} && tail -n +3 {input.tabledeseq2} | sort -r -n -t',' -k 7) > {output.tablesorted};  awk -F ',' 'NR==1; (NR>1) && ($7 < 0.05 )' {output.tablesorted} > {output.tablesignificant};")
