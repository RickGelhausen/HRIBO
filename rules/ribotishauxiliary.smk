rule ribotishGFF:
    input:
        expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples["condition"])
    output:
        "tracks/{condition, [a-zA-Z]+}.ribotish.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/ribotish.py {input} --output_gff3_filepath {output}"
