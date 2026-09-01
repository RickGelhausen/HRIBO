rule mergeConditions:
    input:
        reparation="tracks/{condition}.reparation.gff"
    output:
        "tracks/{condition}.merged.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/concatenate_gff.py {input.reparation:q} -o {output:q}"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", condition=conditions)
    output:
        "tracks/all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/concatenate_gff.py {input.mergedGff:q} -o {output:q}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        "tracks/reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/merge_duplicates_reparation.py -i {input:q} -o {output:q}"

rule reannotatedORFs:
    input:
        annotation=rules.checkAnnotation.output,
        reparation="tracks/reparation.gff"
    output:
        "tracks/reparation_annotated.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/reannotate_orfs.py -a {input.annotation:q} -c {input.reparation:q} -o {output:q}"

rule uniteAnnotation:
    input:
        "tracks/totalAnnotation.gff"
    output:
        "tracks/updated_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/annotation_unite.py -a {input} -o {output}"
