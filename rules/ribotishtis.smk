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

rule ribotishQualityTIS:
   input:
       fp="maplink/TIS/{condition}-{replicate}.bam",
       genome=rules.retrieveGenome.output,
       annotation=rules.ribotishAnnotation.output,
       samindex=rules.genomeSamToolsIndex.output,
       bamindex="maplink/TIS/{condition}-{replicate}.bam.bai"
   output:
       reportpdf="ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual_tis.pdf",
       reporttxt=report("ribotish/{condition, [a-zA-Z]+}-{replicate,\d+}-qual_tis.txt", caption="../report/ribotishquality.rst", category="Ribotish"),
       offsetdone="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.qualdone"
   params:
       offsetparameters="maplink/TIS/{condition, [a-zA-Z]+}-{replicate,\d+}.bam.para.py"
   threads: 10
   log:
       "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_ribotishquality.log"
   shell:
       "conda activate /scratch/bi03/egg/miniconda3/envs/ribotish; mkdir -p ribotish; ribotish quality -v -p {threads} -b {input.fp} -t -g {input.annotation} -o {output.reporttxt} -f {output.reportpdf} 2> {log} || true; if grep -q \"offdict = {{'m0': {{}}}}\" {params.offsetparameters}; then mv {params.offsetparameters} {params.offsetparameters}.unused; fi; touch {output.offsetdone}"

rule ribotish:
    input:
        tis=lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.bam", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"]),
        genome=rules.retrieveGenome.output,
        annotation=rules.ribotishAnnotation.output,
        samindex=rules.genomeSamToolsIndex.output,
	tisindex= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.bam.bai", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"]),
        offsetparameters= lambda wildcards: expand("maplink/TIS/{{condition}}-{replicate}.qualdone", zip, replicate=samples.loc[(samples["method"] == "TIS") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        report="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv_all.txt",
        filtered="ribotish/{condition, [a-zA-Z]+}-newORFs.tsv"
    params:
        tislist= lambda wildcards, input: ','.join(list(set(input.tis))),
        codons= lambda wildcards: ("" if not CODONS else (" --alt --altcodons " + CODONS)),
    threads: 10
    log:
        "logs/{condition, [a-zA-Z]+}_ribotish.log"
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/ribotish; mkdir -p ribotish; ribotish predict --harr -v {params.codons} -p {threads} -t {params.tislist} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"
