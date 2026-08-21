rule newAnnotationReparationOnly:
    input:
        reparation_orfs="tracks/reparation_annotated.gff",
        currentAnnotation=rules.checkAnnotation.output
    output:
        "tracks/totalAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p tracks;
        {SCRIPTS}/concatenate_gff.py {input.reparation_orfs} {input.currentAnnotation} -o {output}
        """
