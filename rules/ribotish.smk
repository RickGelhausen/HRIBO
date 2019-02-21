rule ribomaplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink}"

rule ribobamindexlink:
    input:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, RIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RIBO/; ln -s {params.inlink} {params.outlink}"

rule rnaribomaplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RNARIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNARIBO/; ln -s {params.inlink} {params.outlink}"

rule rnaribobamindexlink:
    input:
        "bam/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, RNARIBO}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNARIBO/; ln -s {params.inlink} {params.outlink}"

rule ribotishQuality:
    input:
        fp="maplink/RIBO/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.retrieveAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex="maplink/RIBO/{condition}-{replicate}.bam.bai"
    output:
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.pdf",
        reporttxt=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.txt", caption="../report/ribotishquality.rst", category="Ribotish"),
        offsetdone="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
    params:
        offsetparameters="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py"
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishquality.log"
    shell:
        "mkdir -p ribotish; ribotish quality -v --th 0.2 -p {threads} -b {input.fp} -g {input.annotation} -o {output.reporttxt} -f {output.reportpdf} 2> {log}; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"


rule ribotish:
    input:
        fp= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.bam", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        genome=rules.retrieveGenome.output,
        annotation=rules.retrieveAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.bam.bai", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        offsetparameters= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.qualdone", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        report="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt",
        #report=report("ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt", caption="../report/ribotish.rst", category="Ribotish"),
        filtered="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv"
    params:
        fplist= lambda wildcards, input: ','.join(list(set(input.fp))),
        codons= lambda wildcards: ("" if not CODONS else (" --alt --altcodons " + CODONS)),
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}_ribotish.log"
    shell:
        "mkdir -p ribotish; ribotish predict --longest -v {params.codons} -p {threads} -b {params.fplist} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"

rule ribotishreport:
    input:
        "ribotish/{condition}-{replicate}-qual.pdf"
    output:
        "figures/{condition}-{replicate}-qual.jpg"
        #report("figures/{condition}-{replicate}-qual.jpg", caption="../report/ribotishquality.rst", category="Novel ORFs - Ribotish")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input} -quality 100 -flatten -sharpen 0x1.0 {output};")
