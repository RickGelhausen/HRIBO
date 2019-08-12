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
        "mkdir -p sammulti; segemehl.x -e -d {input.genome} -i {input.genomeSegemehlIndex} -q {input.fastq1} -p {input.fastq2} --threads {threads} -o {output.sammulti} 2> {log}"
