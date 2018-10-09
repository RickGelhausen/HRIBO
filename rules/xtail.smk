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
        expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples())
    output:
        "normalization/sfactors.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p normalization; SPtools/scripts/generate_size_factors.R -t SPtools/samples.tsv -b maplink/ -a {input[0]} -s {output};")

rule cdsNormalizedCounts:
    input:
        bam=expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        annotation=rules.longestTranscript.output,
        sizefactor="normalization/sfactors.csv"
    output:
        norm="normalization/norm_CDS_reads.csv",
        raw="normalization/raw_CDS_reads.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p normalization; SPtools/scripts/generate_normalized_counts_CDS.R -b maplink/ -a {input.annotation} -s {input.sizefactor} -t SPtools/samples.tsv -n {output.norm} -r {output.raw};")

rule contrastInput:
    #params:
    #    contrasts=expand("{contrastpair}", contrastpair=lambda wildcards: getContrast(wildcards))
    output:
        "contrasts/{contrast}"
    run:
        if not os.path.exists("contrasts"):
            os.makedirs("contrasts")
        for f in getContrast(wildcards):
            print(f)
            open((f), 'a').close()

rule cdsxtail:
    input:
        normreads="normalization/norm_CDS_reads.csv",
        sizefactor="normalization/sfactors.csv",
        contrastfile="contrasts/{contrast}"
    output:
        table="xtail/{contrast}.csv",
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf"
    conda:
        "../envs/xtail.yaml"
    #params:
        #contrast="-c {}".format(lambda wildcards: getContrast(wildcards))
        #contrast=expand("{contrastpair}", contrastpair=lambda wildcards: getContrast(wildcards))
    threads: 1
    shell: ("mkdir -p xtail; SPtools/scripts/xtail_normalized_counts.R -c {input.contrastfile} -t SPtools/samples.tsv -r {input.normreads} -x {output.table} -f {output.fcplot} -p {output.rplot};")

rule riborex:
    input:
        rawreads="normalization/raw_CDS_reads.csv",
        sizefactor="normalization/sfactors.csv",
        contrastfile="contrasts/{contrast}"
    output:
        tabledeseq2="riborex/{contrast}.deseq2.csv",
        tableedger="riborex/{contrast}.edger.csv",
        tablevoom="riborex/{contrast}.voom.csv"
    conda:
        "../envs/riborex.yaml"
    #params:
        #contrast="-c {}".format(lambda wildcards: getContrast(wildcards))
        #contrast=expand("{contrastpair}", contrastpair=lambda wildcards: getContrast(wildcards))
    threads: 1
    shell: ("mkdir -p riborex; SPtools/scripts/xtail_normalized_counts.R -c {input.contrastfile} -t SPtools/samples.tsv -r {input.rawreads} -x {output.tabledeseq2} -y {output.tableedger} -z {output.tablevoom};")
