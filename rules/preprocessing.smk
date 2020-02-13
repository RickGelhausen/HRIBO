rule retrieveGenome:
    input:
        "genome.fa"
    output:
        report("genomes/genome.fa", caption="../report/genome.rst", category="Genome")
    threads: 1
    shell:
        "mkdir -p genomes; mv genome.fa genomes/"

rule retrieveAnnotation:
    input:
        "annotation.gff"
    output:
        report("annotation/annotation.gff", caption="../report/annotation.rst", category="Annotation")
    threads: 1
    shell:
        "mkdir -p annotation; mv annotation.gff annotation/"

