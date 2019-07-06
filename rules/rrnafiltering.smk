rule rrnaannotation:
    input:
        annotation=rules.retrieveAnnotation.output
    output:
        annotation="annotation/rrna.gtf"
    conda:
        "../envs/gawk.yaml"
    threads: 1
    shell:
        "mkdir -p index/annotation; cat {input.annotation} | awk '{{if ($3 == \"rRNA\") print $0;}}' > {output.annotation}"

rule rrnafilter2:
    input:
        mapuniq="rRNAbam/{method}-{condition}-{replicate}.bam",
        annotation="annotation/rrna.gtf"	
    output:
        bam="bam/{method}-{condition}-{replicate}.bam"
    conda:
        "../envs/bedtools.yaml"
    threads: 20
    shell:
        "mkdir -p norRNA; mkdir -p mapuniqnorrna; bedtools intersect -v -a {input.mapuniq} -b {input.annotation} > {output.bam}"

