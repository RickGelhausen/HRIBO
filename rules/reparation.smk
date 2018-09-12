rule reparation:
    input:
        sam="{method}-{condition}-{replicate}.sam",
        genome=rules.retrieveGenome.output,
        gtf=rules.retrieveAnnotation.output,
        db=
    output:

    conda:
        "../envs/reparation.yaml"
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_reparation.log"
    shell:
        reparation.pl -sam {input.sam} -g {input.genome} -gtf {input.gtf} -db {input.db} -out "reparation"
