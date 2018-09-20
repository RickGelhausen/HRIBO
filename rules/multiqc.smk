rule multiqc:
    input:
        expand("tracks/{sample.condition}.ribotish.gff", sample=samples.itertuples()),
        expand("tracks/{sample.method}-{sample.condition}-{sample.replicate}.bw", sample=samples.itertuples())
    output: 
        "QC/Multi/multiqc_report.html"
    params: 
        dir="QC/Multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc -f --exclude picard --exclude gatk -k json -z -o {params.dir} . 2> {log}"
