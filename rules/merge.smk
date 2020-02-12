rule mergeConditions:
    input:
        reparation="tracks/{condition}.reparation.gff"
    output:
        report("tracks/{condition}.merged.gff", caption="../report/novelmerged.rst", category="Novel ORFs")
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.reparation} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", zip, condition=set(samples["condition"]))
    output:
        report("tracks/all.gff", caption="../report/novelall.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/concatenate_gff.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        report("tracks/combined.gff", caption="../report/combined.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/merge_duplicates_reparation.py -i {input} -o {output}"

rule reannotatedORFs:
    input:
        annotation=rules.retrieveAnnotation.output,
        combined="tracks/combined.gff"
    output:
        report("tracks/combined_annotated.gff", caption="../report/combined_annotation.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/reannotate_orfs.py -a {input.annotation} -c {input.combined} -o {output}"

rule newAnnotation:
    input:
        newOrfs="tracks/combined_annotated.gff",
        currentAnnotation=rules.retrieveAnnotation.output
    output:
        "xtail/totalAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/concatenate_gff.py {input.newOrfs} {input.currentAnnotation} -o {output}"

rule uniteAnnotation:
    input:
        "xtail/totalAnnotation.gff"
    output:
        "xtail/newAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/annotation_unite.py -a {input} -o {output}"
