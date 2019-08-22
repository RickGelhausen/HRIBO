rule pairedmap:
    input:
        genome=rules.retrieveGenome.output,
        genomeSegemehlIndex="genomeSegemehlIndex/genome.idx",
        fastq1="trimmed/{method}-{condition}-{replicate}_Q.fastq",
        fastq2="trimmed/{method}-{condition}-{replicate}_P.fastq"
    output:
        sammulti="sammulti/{method}-{condition}-{replicate}.sam"
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log:
        "logs/{method}-{condition}-{replicate}_segemehl.log"
    shell:
        """
        mkdir -p sammulti; segemehl.x -e -d {input.genome} -i {input.genomeSegemehlIndex} -q {input.fastq1} -p {input.fastq2} --threads {threads} -o {output.sammulti} 2> {log}
        """

rule samuniq:
    input:
        sammulti="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        sam="sam/{method}-{condition}-{replicate}.rawsam"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    shell:
        """
        set +e
        mkdir -p sam
        samtools view -f 0x4 {input.sammulti} > {input.sammulti}.unmapped
        samtools view -H <(cat {input.sammulti}) | grep '@HD' > {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V >> {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@RG' >> {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@PG' >> {output.sam}
        samtools view -f 0x2 {input.sammulti} |grep -v '^@' | grep -w 'NH:i:1' >> {output.sam}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """
