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
        "echo \"Computing Segemehl index\"; segemehl.x --threads {threads} -x {output.index:q} -d {input.genome:q} 2> {log:q}"



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
        segemehl.x -d {input.genome:q} -i {input.genomeSegemehlIndex:q} -q {input.fastq:q} --threads {threads} -o {output.sammulti:q} 2> {log:q}
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
        samtools view -H {input.sammulti:q} | awk '$1 == "@HD"' > {output.sam:q}
        samtools view -H {input.sammulti:q} | awk '$1 == "@SQ"' | sort -t$'\t' -k1,1 -k2,2V >> {output.sam:q}
        samtools view -H {input.sammulti:q} | awk '$1 == "@RG"' >> {output.sam:q}
        samtools view -H {input.sammulti:q} | awk '$1 == "@PG"' >> {output.sam:q}
        awk '$1 !~ /^@/ && int($2 / 4) % 2 == 0 {{for (i=12; i<=NF; i++) if ($i == "NH:i:1") {{found=1; print; break}}}} END {{if (!found) exit 1}}' {input.sammulti:q} >> {output.sam:q}
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
        "cp {input.sam:q} {output.sam:q}"

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
        "samtools view -@ {threads} -bh {input.sam:q} | samtools sort -@ {threads} -o {output:q} -O bam"

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
        "samtools view -@ {threads} -bh {input.sam:q} | samtools sort -@ {threads} -o {output:q} -O bam"

rule maplink:
    input:
        bam="bam/{method}-{condition}-{replicate}.bam",
        stager=str(SCRIPTS / "stage_input.py")
    output:
        bam="maplink/{method}-{condition}-{replicate}.bam"
    threads: 1
    shell:
        "python3 {input.stager:q} portable-link {input.bam:q} {output.bam:q}"
