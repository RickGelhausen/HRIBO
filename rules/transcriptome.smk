rule transcripts:
    input:
        rules.maplink.output,
        rules.genomeSize.output
    output:
        "transcriptome/RNA-{condition}-{replicate}.gtf"
    conda:
        "../envs/scallop.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; scallop -i maplink/{params.prefix}.bam -o {output}"
