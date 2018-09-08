rule genomeSegemehlIndex:
    input:
        genome=rules.retrieveGenome.output
    output:
        index="genomeSegemehlIndex/genome.idx"
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    params:
        indexpath=lambda wildcards: ("NOTSET" if not INDEXPATH else (INDEXPATH))
    log:
        "logs/genomeIndex.log"
    shell:
        "if [ -d {params.indexpath} ]; then ln -T -s {params.indexpath} {output.index}; echo \"Index linked\"; else mkdir -p genomeSegemehlIndex; echo \"Computing Segemehl index\"; segemehl.x --threads {threads} -x {output.index} -d {input.genome} 2> {log}; fi"

rule map:
    input:
        genome=rules.retrieveGenome.output,
        genomeSegemehlIndex="genomeSegemehlIndex/genome.idx",
        fastq="norRNA/{method}-{condition}-{replicate}.fastq",
    output:
        sam="sam/{method}-{condition}-{replicate}.sam",
        unmapped="sam/unmapped/{method}-{condition}-{replicate}.sam"
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log:
        "logs/{method}-{condition}-{replicate}_segemehl.log"
    shell:
        "mkdir -p sam; mkdir -p sam/unmapped; segemehl.x -s -d {input.genome} -i {input.genomeSegemehlIndex} -q {input.fastq} --threads {threads} -o {output.sam} -u {output.unmapped} 2> {log}"

rule samtobam:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        bam="bam/{method}-{condition}-{replicate}.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        "mkdir -p bam; samtools view -bh {input.sam} | samtools sort -o {output.bam} -O bam"

rule maplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method}-{condition}-{replicate}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink; ln -s {params.inlink} {params.outlink}"

rule ribomaplink:
    input:
        "bam/{method}-{condition}-{replicate}/Aligned.sortedByCoord.out.bam"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink}"

rule rnamaplink:
    input:
        "bam/{method}-{condition}-{replicate}/Aligned.sortedByCoord.out.bam"
    output:
        "maplink/{method, RNA}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNA/; ln -s {params.inlink} {params.outlink}"
