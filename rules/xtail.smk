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
        "xtail/sfactors.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p xtail; SPtools/scripts/generate_size_factors.R -t SPtools/samples.tsv -b maplink/ -a {input[0]} -s {output};")

rule cdsNormalizedCounts:
    input:
        bam=expand("maplink/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples()),
        annotation=rules.longestTranscript.output,
        sizefactor="xtail/sfactors.csv"
    output:
        "xtail/norm_CDS_reads.csv"
    conda:
        "../envs/normalization.yaml"
    threads: 1
    shell: ("mkdir -p xtail; SPtools/scripts/generate_normalized_counts_CDS.R -b maplink/ -a {input.annotation} -s {input.sizefactor} -t SPtools/samples.tsv -n {output};")

rule cdsxtail:
    input:
        "xtail/norm_CDS_reads.csv"
    output:
        table=report("xtail/xtail.csv", caption="../report/xtail_cds_fc.rst", category="CDS"),
        fcplot="xtail/xtail_cds_fc.pdf",
        rplot="xtail/xtail_cds_r.pdf"
    conda:
        "../envs/xtail.yaml"
    threads: 1
    shell: ("mkdir -p xtail; SPtools/scripts/xtail_normalized_counts.R -t SPtools/samples.tsv -r {input} -x {output.table} -f {output.fcplot} -p {output.rplot};")
