def getpairedfastq(wildcards):
    print(wildcards)
    return (samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["fastqFile","fastqFile2"]].dropna())

#def adaptersstring(wildcards):

rule pairlinktrim:
    input:
        fastq=getpairedfastq
    output:
        fastq1="trimlink/{method}-{condition}-{replicate}_Q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_P.fastq.gz"
    params:
        prefix1=lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.basename(input.fastq[0]))[0]),
        prefix2=lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.basename(input.fastq[1]))[0]),
        inlink1=lambda wildcards, input:(os.getcwd() + "/" + str(input[0])),
        inlink2=lambda wildcards, input:(os.getcwd() + "/" + str(input[1])),
        outlink1=lambda wildcards, output:(os.getcwd() + "/" + str(output.fastq1)),
        outlink2=lambda wildcards, output:(os.getcwd() + "/" + str(output.fastq2))
    threads: 1
    shell:
        "mkdir -p trimlink; ln -s {params.inlink1} {params.outlink1}; ln -s {params.inlink2} {params.outlink2};"

rule pairtrim:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_Q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_P.fastq.gz"
    output:
        fastq1="trimmed/{method}-{condition}-{replicate}_Q.fastq",
        fastq2="trimmed/{method}-{condition}-{replicate}_P.fastq"
    params:
        adapter=lambda wildcards, output: ("" if not ADAPTERS else (" ".join([" -a %s" % adapter for adapter in ADAPTERS.split(",")]))),
        quality=" -q 20 --trim-n ",
        filtering=" -m 10 "
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    shell:
        "mkdir -p trimmed; cutadapt -j {threads} {params.adapter} {params.quality} {params.filtering} -o {output.fastq1} {output.fastq2} {input.fastq1} {input.fastq2}"
