rule processAnnotation:
    input:
        annotation=rules.retrieveAnnotation.output
    output:
        "annotation/processed-annotation.gtf"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p annotation; python3 ribo_benchmark/scripts/processAnnotation.py -a {input.annotation} -o {output}"

rule generateMetageneRoi:
    input:
        annotation=rules.processAnnotation.output
    output:
        "offsets/metagene_rois.txt"
    conda:
        "../envs/auxillary.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; metagene generate -q offsets/metagene --landmark cds_start --annotation_files {input.annotation}"

rule psiteOffsets:
    input:
        rois=rules.generateMetageneRoi.output,
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bamindex="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        psite="offsets/{method}-{condition}-{replicate}_p_offsets.txt"
    conda:
        "../envs/auxillary.yaml"
    threads: 1
    params:
        prefix=lambda wildcards: "{wildcards.method}-{wildcards.condition}-{wildcards.replicate}"
    shell:
        "mkdir -p offsets; psite {input.rois} offsets/{params.prefix} --min_length 22 --max_length 40 --require_upstream --count_files {input.bam}"
