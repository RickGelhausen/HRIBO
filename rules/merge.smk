rule mergeConditions:
    input:
        ribotish="tracks/{condition}.ribotish.gff",
        reparation="tracks/{condition}.reparation.gff"
    output:
        report("tracks/{condition, [a-zA-Z]+}.merged.gff", caption="../report/novelmerged.rst", category="Novel ORFs")
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.ribotish} > {output}.unsorted; cat {input.reparation} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", zip, condition=set(samples["condition"])))
    output:
        report("tracks/all.gff", caption="../report/novelall.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        report("tracks/combined.gff", caption="../report/combined.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/noverlapper.py -i {input} -o {output}"

rule newAnnotation:
    input:
        newOrfs="tracks/combined.gff",
        currentAnnotation="annotation/annotation.gtf"
    output:
        "xtail/newAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input.newOrfs} {input.currentAnnotation} -o {output}"
