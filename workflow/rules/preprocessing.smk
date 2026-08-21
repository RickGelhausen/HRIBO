rule retrieveGenome:
    input:
        genome=config["biologySettings"]["genome"]
    output:
        "genomes/genome.fa"
    threads: 1
    shell:
        "cp {input.genome} genomes/genome.fa"

rule retrieveAnnotation:
    input:
        annotation=config["biologySettings"]["annotation"]
    output:
        "annotation/annotation.gff"
    threads: 1
    shell:
        "cp {input.annotation} annotation/annotation.gff"

rule checkAnnotation:
    input:
        rules.retrieveAnnotation.output
    output:
        "annotation/annotation_processed.gff"
    threads: 1
    shell:
        "{SCRIPTS}/gtf2gff3.py -a {input} -o {output}"
