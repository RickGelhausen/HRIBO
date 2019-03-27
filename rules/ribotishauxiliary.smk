rule ribotishGFF:
    input:
        "ribotish/{condition}-newORFs.tsv_all.txt"
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotish.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/ribotish.py {input} --condition {wildcards.condition} --output_gff3_filepath {output}"

rule ribotishAnnotation:
    input:
        annotation="qc/featurecount/annotation.gtf"
    output:
        "ribotish/annotation_processed.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p ribotish; SPtools/scripts/createRiboTISHannotation.py -a {input.annotation} -o {output}"
