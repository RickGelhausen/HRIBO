rule rrnaannotation:
    input:
        annotation=rules.retrieveAnnotation.output
    output:
        annotation="annotation/rrna.bed"
    conda:
        "../envs/gawk.yaml"
    threads: 1
    shell:
        """
        mkdir -p annotation; awk -F'\\t' '$3 == "rRNA"' annotation/annotation.gtf | awk -F'\\t' '{{print $1 FS $4 FS $5 FS $6 FS $8 FS $7}}' > {output.annotation}
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
