import itertools as iter

def getcontrast(wildcards):
  conditions=wildcards.samples['condition'].tolist()
  contrasts=iter.combinations(conditions, 2)
return contrasts

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
        normreads="xtail/norm_CDS_reads.csv"
        contrast=getcontrast
    output:
        table="xtail/{contrast}.csv",
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf"
    conda:
        "../envs/xtail.yaml"
    threads: 1
    shell: ("mkdir -p xtail; SPtools/scripts/xtail_normalized_counts.R -c {input.contrast} -t SPtools/samples.tsv -r {input.normreads} -x {output.table} -f {output.fcplot} -p {output.rplot};")
