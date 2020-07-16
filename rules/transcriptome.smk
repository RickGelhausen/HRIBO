rule transcripts:
    input:
        "maplink/{method}-{condition}-{replicate}.bam",
    output:
        "transcriptome/{method}-{condition}-{replicate}.gtf"
    conda:
        "../envs/scallop.yaml"
    #wildcard_constraints:
    #    method=["RNA"]
    threads: 5
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; stringtie -p {threads} --conservative -l 150 -o {output} maplink/{params.prefix}.bam"

rule mergedTranscripts:
    input:
        lambda wildcards: list(set(expand("transcriptome/{method}-{condition}-{replicate}.gtf", zip, method=samples.loc[(samples["method"] == "RNA"), "method"], condition=samples["condition"], replicate=samples["replicate"]))),
        rules.retrieveAnnotation.output
    output:
        "transcriptome/all.combined.gtf"
    conda:
        "../envs/scallop.yaml"
    threads: 5
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "mkdir -p transcriptome; stringtie -p {threads} --merge {input[0]} -G {input[1]} -C -o {output}"
        #"mkdir -p transcriptome; echo \"{input}\" | tr \" \" \"\\n\" > transcriptome/inputlist;  gffcompare -i transcriptome/inputlist  -o transcriptome/all"

rule mergedAnnotationTranscripts:
    input:
        "transcriptome/all.combined.gtf",
        rules.retrieveAnnotation.output
    output:
        "transcriptome/all_annotation.combined.gtf"
    conda:
        "../envs/gffcompare.yaml"
    threads: 1
    shell:
        "mkdir -p transcriptome; echo \"{input}\" | tr \" \" \"\\n\" > transcriptome/inputlist2;  gffcompare -i transcriptome/inputlist2 -r {input[1]}  -o transcriptome/all_annotation"
