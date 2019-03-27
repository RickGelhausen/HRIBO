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

rule generateMetageneRoiStart:
    input:
        annotation=rules.processAnnotation.output
    output:
        "offsets/metagene_start_rois.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; metagene generate -q offsets/metagene_start --landmark cds_start --annotation_files {input.annotation}"

rule generateMetageneRoiStop:
    input:
        annotation=rules.processAnnotation.output
    output:
        "offsets/metagene_stop_rois.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; metagene generate -q offsets/metagene_stop --landmark cds_stop --annotation_files {input.annotation}"

rule psiteOffsets:
    input:
        rois=rules.generateMetageneRoiStart.output,
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bamindex="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        psite="offsets/{method}-{condition}-{replicate}_p_offsets.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    params:
        prefix=lambda wildcards: "{wildcards.method}-{wildcards.condition}-{wildcards.replicate}"
    shell:
        "mkdir -p offsets; psite {input.rois} offsets/{params.prefix} --min_length 22 --max_length 40 --require_upstream --count_files {input.bam}"
