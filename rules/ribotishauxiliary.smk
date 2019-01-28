rule ribotishGFF:
    input:
        expand("ribotish/{sample.condition}-newORFs.tsv_all.txt", sample=samples.itertuples())
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotish.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/ribotish.py {input} --output_gff3_filepath {output}"
