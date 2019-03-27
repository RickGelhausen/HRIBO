def getfastq(wildcards):
    print(wildcards)
    return samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["fastqFile"]].dropna()

rule linktrim:
    input:
        reads=getfastq
    output:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq"
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.reads[0]))[0])[0])
    threads: 1
    shell:
        "mkdir -p trimlink; ln -s fastq/{params.prefix}.fastq {output.fastq};"

rule trim:
    input:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq"
    output:
        fastq="trimmed/{method}-{condition}-{replicate}.fastq",
    params:
        adapter=lambda wildcards, output: ("" if not ADAPTERS else ' '.join([" -a " + s for s in ADAPTERS.split(",")]),
        quality=" -q 20 --trim-n "
    conda:
        "../envs/cutadapt.yaml"
    threads: 1
    shell:
        "mkdir -p trimmed; cutadapt {params.ada} {params.quality}  -o {output.fastq} {input.fastq}"
