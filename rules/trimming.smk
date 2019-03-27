def getfastq(wildcards):
    print(wildcards)
    return samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["fastqFile"]].dropna()

#def adaptersstring(wildcards):  

rule linktrim:
    input:
        fastq=getfastq
    output:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
    params:
        prefix=lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.basename(input.fastq[0]))[0]),
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output.fastq))
    threads: 20
    shell:
        "mkdir -p trimlink; ln -s {params.inlink} {params.outlink};"

rule trim:
    input:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
    output:
        fastq="trimmed/{method}-{condition}-{replicate}.fastq"
    params:
        adapter=lambda wildcards, output: ("" if not ADAPTERS else (" ".join([" -a %s" % adapter for adapter in ADAPTERS.split(",")]))),
        quality=" --trim-n "
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    shell:
        "mkdir -p trimmed; cutadapt -j {threads}  {params.adapter} {params.quality}  -o {output.fastq} {input.fastq}"
