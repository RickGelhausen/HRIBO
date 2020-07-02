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

rule mergedtranscripts:
   input:
        lambda wildcards: expand("transcriptome/RNA-{condition}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RNA") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "transcriptome/transcriptome-combined.gtf"
    conda:
        "../envs/gffcompare.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; cd transcriptome; gffcompare -i {input} -o transcriptome"

rule mergedtranscripts:
   input:
        rules.mergedtranscripts.output,
        rules.retrieveAnnotation.output
    output:
        "transcriptome/transcriptome-annotation-combined.gtf"
    conda:
        "../envs/gffcompare.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; cd transcriptome; gffcompare -i transcriptome/{params.prefix}.gtf -o transcriptome"
     
