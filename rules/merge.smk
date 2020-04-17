rule mergeConditions:
    input:
        reparation="tracks/{condition}.reparation.gff"
    output:
        "tracks/{condition}.merged.gff"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.reparation} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", zip, condition=set(samples["condition"]))
    output:
        "tracks/all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/concatenate_gff.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        "tracks/reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/merge_duplicates_reparation.py -i {input} -o {output}"

rule reannotatedORFs:
    input:
        annotation=rules.retrieveAnnotation.output,
        reparation="tracks/reparation.gff"
    output:
        "tracks/reparation_annotated.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/reannotate_orfs.py -a {input.annotation} -c {input.reparation} -o {output}"

rule uniteAnnotation:
    input:
        "tracks/totalAnnotation.gff"
    output:
        "tracks/updated_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/annotation_unite.py -a {input} -o {output}"
