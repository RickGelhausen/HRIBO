rule xtailclassicnormalize:
    input:
        bam="maplink/{method}-{condition}-{replicate}.bam",
        annotation=rules.retrieveAnnotation.output
    output:
        xtailclassic="xtailclassic/{method}-{condition}-{replicate}.csv"
    conda:
        "../envs/xtailcounts.yaml"
    threads: 1
    shell: ("mkdir -p xtailclassic; python RPF_count_CDS.py {input.bam} {input.annotation} > {output}")
