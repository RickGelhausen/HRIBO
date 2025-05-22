rule rrnaannotation:
    input:
        annotation=rules.checkAnnotation.output
    output:
        annotation="annotation/rrna.bed"
    conda:
        "../envs/gawk.yaml"
    threads: 1
    shell:
        """
        mkdir -p annotation; awk -F'\\t' '$3 == "rRNA" || $3 == "tRNA"' {input.annotation} | awk -F'\\t' '{{print $1 FS $4-1 FS $5 FS "." FS "." FS $7}}' > {output.annotation}
        """

rule rrnafilter2:
    input:
        mapuniq="rRNAbam/{method}-{condition}-{replicate}.bam",
        annotation="annotation/rrna.bed"
    output:
        bam="bam/{method}-{condition}-{replicate}.bam"
    conda:
        "../envs/bedtools.yaml"
    threads: 20
    shell:
        "mkdir -p norRNA; mkdir -p mapuniqnorrna; bedtools intersect -v -a {input.mapuniq} -b {input.annotation} > {output.bam}"
