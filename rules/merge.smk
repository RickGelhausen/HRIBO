rule mergeConditions:
    input:
        ribotish="tracks/{condition}.ribotish.gff",
        reparation="tracks/{condition}.reparation.gff"
    output:
        "tracks/{condition, [a-zA-Z]+}.merged.gff"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.ribotish} > {output}.unsorted; cat {input.reparation} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};"

rule noverlap:
    input:
        mergedGff="tracks/{condition}.merged.gff"
    output:
        "tracks/{condition}.filtered.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/noverlapper.py -i {input.mergedGff} -o {output}"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{sample.condition}.filtered.gff", sample=samples.itertuples())
    output:
        "tracks/all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        "tracks/filtered.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/noverlapper.py -i {input} -o {output}"

rule newAnnotation:
    input:
        newOrfs="tracks/filtered.gff",
        currentAnnotation="xtail/longest_protein_coding_transcripts.gtf"
    output:
        "xtail/newAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input.newOrfs} {input.currentAnnotation} -o xtail/tmp.gff; SPtools/scripts/merge_orfs.py -i xtail/tmp.gff -o {output};"
