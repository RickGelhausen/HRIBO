rule retrieveGenome:
    input:
        genome=config["biologySettings"]["genome"],
        stager=str(SCRIPTS / "stage_input.py")
    output:
        genome="genomes/genome.fa"
    threads: 1
    shell:
        "python3 {input.stager:q} text {input.genome:q} {output.genome:q}"

rule retrieveAnnotation:
    input:
        annotation=config["biologySettings"]["annotation"],
        stager=str(SCRIPTS / "stage_input.py")
    output:
        annotation="annotation/annotation.gff"
    threads: 1
    shell:
        "python3 {input.stager:q} text {input.annotation:q} {output.annotation:q}"

rule checkAnnotation:
    input:
        annotation=rules.retrieveAnnotation.output,
        converter=str(SCRIPTS / "gtf2gff3.py"),
        converter_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        annotation="annotation/annotation_processed.gff"
    threads: 1
    shell:
        "python3 {input.converter:q} -a {input.annotation:q} -o {output.annotation:q}"
