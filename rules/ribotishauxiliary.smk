rule ribotishGFF:
    input:
        "ribotish/{condition}-newORFs.tsv_all.txt"
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotish.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/create_ribotish_gff.py {input} --condition {wildcards.condition} --output_gff3_filepath {output}"

rule ribotishAnnotation:
    input:
        annotation="qc/featurecount/annotation.gtf",
        sizes="genomes/sizes.genome"
    output:
        "ribotish/annotation_processed.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p ribotish; HRIBO/scripts/create_ribotish_annotation.py -a {input.annotation} --genome_sizes {input.sizes} --annotation_output {output}"
