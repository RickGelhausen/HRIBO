rule genomeSegemehlIndex:
    input:
        genome=rules.retrieveGenome.output
    output:
        index=temp("genomeSegemehlIndex/genome.idx")
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=120
    log:
        "logs/genomeIndex.log"
    shell:
        "echo \"Computing Segemehl index\"; segemehl.x --threads {threads} -x {output.index} -d {input.genome} 2> {log}"



rule map:
    input:
        fastq="trimmed/{method}-{condition}-{replicate}.fastq",
        genome=rules.retrieveGenome.output,
        genomeSegemehlIndex="genomeSegemehlIndex/genome.idx"
    output:
        sammulti=temp("sammulti/{method}-{condition}-{replicate}.sam"),
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=240
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log:
        "logs/{method}-{condition}-{replicate}_segemehl.log"
    shell:
        """
        segemehl.x -e -d {input.genome} -i {input.genomeSegemehlIndex} -q {input.fastq} --threads {threads} -o {output.sammulti} 2> {log}
        """

rule samuniq:
    input:
        sammulti="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        sam=temp("sam/{method}-{condition}-{replicate}.rawsam"),
        #unmapped=temp("sammulti/{method}-{condition}-{replicate}.sam.unmapped")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        """
        set +e
        awk '$2 != "4"' {input.sammulti} > {input.sammulti}.mapped
        samtools view -H <(cat {input.sammulti}) | grep '@HD' > {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V >> {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@RG' >> {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@PG' >> {output.sam}
        cat {input.sammulti}.mapped |grep -v '^@' | grep -w 'NH:i:1' >> {output.sam}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

# The strand inverter is for protocols that sequence the antisense strand. No
# method tag currently selects it, so this is a straight copy; the rule is kept
# as the single place to reintroduce that behaviour.
rule samstrandswap:
    input:
        sam="sam/{method}-{condition}-{replicate}.rawsam"
    output:
        sam=temp("sam/{method}-{condition}-{replicate}.sam")
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    shell:
        "cp {input.sam} {output.sam}"

rule sammultitobam:
    input:
        sam="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        temp("bammulti/{method}-{condition}-{replicate}.bam")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "samtools view -@ {threads} -bh {input.sam} | samtools sort -@ {threads} -o {output} -O bam"

rule samtobam:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        temp("rRNAbam/{method}-{condition}-{replicate}.bam")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "samtools view -@ {threads} -bh {input.sam} | samtools sort -@ {threads} -o {output} -O bam"

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
        "ln -s {params.inlink} {params.outlink}"
