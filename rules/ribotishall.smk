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

rule ribotishQualityRIBO:
    input:
        fp="maplink/RIBO/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex="maplink/RIBO/{condition}-{replicate}.bam.bai",
    output:
        offsetdone="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
    params:
        offsetparameters="maplink/RIBO/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py",
        reporttxt="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.txt",
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual.pdf"
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishqualityribo.log"
    shell:
        "mkdir -p ribotish; ribotish quality -v -p {threads} -b {input.fp} -g {input.annotation} -o {params.reporttxt} -f {params.reportpdf} 2> {log}; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"

rule ribotishQualityTIS:
    input:
        tis="maplink/TIS/{condition}-{replicate}.bam",
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        tisindex="maplink/TIS/{condition}-{replicate}.bam.bai",
    output:
        offsetdone="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
    params:
        offsetparameters="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py",
        reporttxt="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual_tis.txt",
        reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual_tis.pdf"
    conda:
        "../envs/ribotish.yaml"
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishqualitytis.log"
    shell:
        "mkdir -p ribotish; ribotish quality -v -p {threads} -b {input.tis} -b {input.tis} -g {input.annotation} -o {params.reporttxt} -f {params.reportpdf} 2> {log}; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"

rule ribotish:
    input:
        fp= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.bam", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        tis= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.bam", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"]),
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
        bamindex= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.bam.bai", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        tisindex= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.bam.bai", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"]),
        ribooffsetparameters= lambda wildcards: expand("maplink/RIBO/{{condition}}-{replicate}.qualdone", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        tisoffsetparameters= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.qualdone", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        report="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt",
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
        "mkdir -p ribotish; ribotish predict --harr -v {params.codons} -p {threads} -t {params.tislist} -b {params.fplist} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"
