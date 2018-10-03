def getfastq(wildcards):
    return samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["fastqFile"]].dropna()

rule trim:
    input:
        reads=getfastq
    output:
        fastq="trimmed/{method}-{condition}-{replicate}.fastq",
        report="trimmed/{method}-{condition}-{replicate}.fastq.gz_trimming_report.txt"
    params:
        ada=lambda wildcards, output: ("" if not ADAPTERS else (" -a " + ADAPTERS)),
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.reads[0]))[0])[0])
    conda:
        "../envs/trimgalore.yaml"
    threads: 1
    shell:
        "mkdir -p trimmed; trim_galore {params.ada} --phred33 -q 20 --length 15 --output_dir trimmed/ --trim-n --suppress_warn --clip_R1 1 --dont_gzip fastq/{params.prefix}.fastq.gz; mv trimmed/{params.prefix}_trimmed.fq {output.fastq}; mv trimmed/{params.prefix}.fastq.gz_trimming_report.txt {output.report}; sed -i -- 's%trimmed/{params.prefix}%{output.fastq}%g' {output.report}"
