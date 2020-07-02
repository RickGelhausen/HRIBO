rule transcripts:
    input:
        rules.maplink.output,
    output:
        "transcriptome/RNA-{condition}-{replicate}.gtf"
    conda:
        "../envs/scallop.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; scallop -i maplink/{params.prefix}.bam -o {output}"

rule mergedTranscripts:
    input:
        lambda wildcards: expand("transcriptome/RNA-{condition}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RNA") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "transcriptome/all-combined.gtf"
    conda:
        "../envs/gffcompare.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; gffcompare -i {input} -o transcriptome/all"

rule mergedAnnotationTranscripts:
   input:
        rules.mergedTranscripts.output,
        rules.retrieveAnnotation.output
    output:
        "transcriptome/all-annotation-combined.gtf"
    conda:
        "../envs/gffcompare.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; gffcompare -i {input[0]} {input[1]} -o transcriptome/all-annotation"
     
