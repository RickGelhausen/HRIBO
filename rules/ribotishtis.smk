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

rule ribotish:
    input:
        fp=expand("maplink/RIBO/{{condition}}-{sample.replicate}.bam", sample=samples.itertuples()),
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
        "mkdir -p ribotish; ribotish predict --longest -v {params.codons} -p {threads} -t {params.tislist} -b {params.fplist} -g {input.annotation} -f {input.genome} -o {output.filtered} 2> {log}"