rule uniprotDBRetrieve:
    input:
        HTTP.remote("TODO.fasta",keep_local=True,allow_redirects=True)
    output:
        "uniprotDB/TODO.fasta"
    threads: 1
    run:
        outputname = os.path.basename(input[0])
        shell("mkdir -p uniprotDBM; mv {input} uniprotDB/{outputName}")

rule reparation:
    input:
        sam="sam/RIBO-{condition}-{replicate}.sam",
        genome=rules.retrieveGenome.output,
        gtf=rules.retrieveAnnotation.output,
        db="uniprotDB/TODO.fasta"
    output:
        bed="reparation/{condition}-{replicate}/Predicted_ORFs.bed"
    conda:
        "../envs/reparation.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]))
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_reparation.log"
    shell:
        mkdir -p reparation; reparation.pl -sam {input.sam} -g {input.genome} -gtf {input.gtf} -db {input.db} -out {params.prefix} -threads {threads}
