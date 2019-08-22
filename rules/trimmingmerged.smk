def getpairedfastq(wildcards):
    print(wildcards)
    return (samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["fastqFile","fastqFile2"]].dropna())

# rule mergePairs:
#     input:
#         fastq=getpairedfastq
#     output:
#         fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
#     params:
#         reads1=lambda wildcards, input: input.fastq[0],
#         reads2=lambda wildcards, input: input.fastq[1]
#     threads: 5
#     conda:
#         "../envs/dc.yaml"
#     shell:
#         """
#         mkdir -p trimlink
#         flash2 {params.reads1} {params.reads2} -t {threads} -z -o trimlink/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} -M 100
#         mv trimlink/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}.extendedFrags.fastq.gz {output.fastq}
#         """
#
# rule trim:
#     input:
#         fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
#     output:
#         fastq="trimmed/{method}-{condition}-{replicate}.fastq"
#     params:
#         adapter=lambda wildcards, output: ("" if not ADAPTERS else (" ".join([" -a %s" % adapter for adapter in ADAPTERS.split(",")]))),
#         quality=" -q 20 --trim-n ",
#         filtering=" -m 10 "
#     conda:
#         "../envs/cutadapt.yaml"
#     threads: 20
#     shell:
#         "mkdir -p trimmed; cutadapt -j {threads} {params.adapter} {params.quality} {params.filtering} -o {output.fastq} {input.fastq}"

rule linktrim:
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

rule trim:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_Q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_P.fastq.gz"
    output:
        fastq1="tmp_trimmed/{method}-{condition}-{replicate}_Q.fastq",
        fastq2="tmp_trimmed/{method}-{condition}-{replicate}_P.fastq"
    params:
        adapter=lambda wildcards, output: ("" if not ADAPTERS else (" ".join([" -a %s" % adapter for adapter in ADAPTERS.split(",")]))),
        quality=" -q 20 --trim-n ",
        filtering=" -m 10 "
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    shell:
        "mkdir -p tmp_trimmed; cutadapt -j {threads} {params.adapter} {params.quality} {params.filtering} -o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq2}"

rule mergeTrim:
    input:
        read1="tmp_trimmed/{method}-{condition}-{replicate}_Q.fastq",
        read2="tmp_trimmed/{method}-{condition}-{replicate}_P.fastq"
    output:
        fastq="trimmed/{method}-{condition}-{replicate}.fastq"
    threads: 5
    conda:
        "../envs/dc.yaml"
    shell:
        """
        mkdir -p trimmed
        flash2 {input.read1} {input.read2} -t {threads} -o trimmed/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} -M 100
        mv trimmed/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}.extendedFrags.fastq {output.fastq}
        """
