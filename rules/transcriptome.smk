rule transcripts:
    input:
        "maplink/{method}-{condition}-{replicate}.bam",
    output:
        "transcriptome/{method}-{condition}-{replicate}.gtf"
    conda:
        "../envs/scallop.yaml"
    #wildcard_constraints:
    #    method=["RNA"]
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; scallop -i maplink/{params.prefix}.bam -o {output}"

rule mergedTranscripts:
    input:
        lambda wildcards: expand("transcriptome/RNA-{condition}-{replicate}.gtf", zip, condition=samples["condition"], replicate=samples["replicate"])
    output:
        "transcriptome/all_combined.gtf"
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
        "transcriptome/all_annotation_combined.gtf"
    conda:
        "../envs/gffcompare.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; gffcompare -i {input[0]} {input[1]} -o transcriptome/all_annotation"
     
