def updated_annotation_predictions(_wildcards):
    predictions = ["tracks/reparation_annotated.gff"]
    if DEEPRIBO:
        predictions.append("tracks/deepribo_merged_plus.gff")
    return predictions


rule mergeConditions:
    input:
        reparation="tracks/{condition}.reparation.gff",
        script=str(SCRIPTS / "concatenate_gff.py")
    output:
        "tracks/{condition}.merged.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} {input.reparation:q} -o {output:q}"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.merged.gff", condition=conditions),
        script=str(SCRIPTS / "concatenate_gff.py")
    output:
        "tracks/all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} {input.mergedGff:q} -o {output:q}"

rule filterAll:
    input:
        gff="tracks/all.gff",
        script=str(SCRIPTS / "merge_duplicates_reparation.py"),
        script_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        "tracks/reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -i {input.gff:q} -o {output:q}"

rule reannotatedORFs:
    input:
        annotation=rules.checkAnnotation.output,
        reparation="tracks/reparation.gff",
        script=str(SCRIPTS / "reannotate_orfs.py"),
        script_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        "tracks/reparation_annotated.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -a {input.annotation:q} -c {input.reparation:q} -o {output:q}"

rule updatedAnnotation:
    input:
        annotation=rules.checkAnnotation.output,
        predictions=updated_annotation_predictions,
        script=str(SCRIPTS / "build_updated_annotation.py"),
        script_deps=[str(SCRIPTS / "concatenate_gff.py")]
    output:
        "tracks/updated_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} --annotation {input.annotation:q} --predictions {input.predictions:q} --output {output:q}"
