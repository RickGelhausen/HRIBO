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
        awk -F'\\t' '$3 == "rRNA" || $3 == "tRNA"' {input.annotation:q} | awk -F'\\t' '{{print $1 FS $4-1 FS $5 FS "." FS "." FS $7}}' > {output.annotation:q}
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
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "bedtools intersect -v -a {input.mapuniq:q} -b {input.annotation:q} > {output.bam:q}"
