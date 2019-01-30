rule tismaplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, TIS}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/TIS/; ln -s {params.inlink} {params.outlink}"

rule tisbamindexlink:
    input:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, TIS}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/TIS/; ln -s {params.inlink} {params.outlink}"

rule rnatismaplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method, RNATIS}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNATIS/; ln -s {params.inlink} {params.outlink}"

rule rnatisbamindexlink:
    input:
        "bam/{method}-{condition}-{replicate}.bam.bai"
    output:
        "maplink/{method, RNATIS}/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.bai"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink/RNATIS/; ln -s {params.inlink} {params.outlink}"

rule ribotishQuality:
    input:
        fp="maplink/TIS/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.retrieveAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex="maplink/TIS/{condition}-{replicate}.bam.bai"
    output:
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.pdf",
        reporttxt=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.txt", caption="../report/ribotishquality.rst", category="Ribotish"),
        offsetdone="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
    params:
        offsetparameters="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py"
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishquality.log"
    shell:
        "mkdir -p ribotish; ribotish quality -v --th 0.2 -p {threads} -b {input.fp} -t -g {input.annotation} -o {output.reporttxt} -f {output.reportpdf} 2> {log}; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"

rule ribotish:
    input:
	tis=expand("maplink/TIS/{{condition}}-{sample.replicate}.bam", sample=samples.itertuples()),
        genome=rules.retrieveGenome.output,
        annotation=rules.retrieveAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex=expand("maplink/RIBO/{{condition}}-{sample.replicate}.bam.bai", sample=samples.itertuples()),
        offsetparameters=expand("maplink/RIBO/{{condition}}-{sample.replicate}.qualdone", sample=samples.itertuples())
        #offsetparameters="maplink/RIBO/{condition}-{replicate}.bam.para.py"
    output:
        report="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt",
        #report=report("ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt", caption="../report/ribotish.rst", category="Ribotish"),
        filtered="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv"
    params:
        fplist= lambda wildcards, input: ','.join(list(set(input.fp))),
	tislist= lambda wildcards, input: ','.join(list(set(input.tis))),
        codons= lambda wildcards: ("" if not CODONS else (" --alt --altcodons " + CODONS)),
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}_ribotish.log"
    shell:
        "mkdir -p ribotish; ribotish predict --longest -v {params.codons} -p {threads} -t {params.tislist} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"
