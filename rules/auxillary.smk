def getGFFtype(filename):
    try:
        with open(filename, "r") as f:
            first_line = f.readline().strip()
        
        if first_line == "##gff-version 3":
            return "GFF3"
        else:
            return "GTF2"
    except FileNotFoundError:
        return "failed"

rule processAnnotation:
    input:
        annotation=rules.retrieveAnnotation.output
    output:
        "offsets/processed-annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p offsets; SPtools/scripts/processAnnotation.py -a {input.annotation} -o {output}"

rule generateMetageneRoiStart:
    input:
        rules.processAnnotation.output
    output:
        "offsets/metagene_start_rois.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    params:
        gfftype=lambda wildcards, input: (getGFFtype(str(input)))
    shell:
        "mkdir -p offsets; echo {params.gfftype}; metagene generate -q offsets/metagene_start --landmark cds_start --annotation_files {input} --annotation_format {params.gfftype}"

rule generateMetageneRoiStop:
    input:
        rules.processAnnotation.output
    output:
        "offsets/metagene_stop_rois.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    params:
        gfftype=lambda wildcards, input: (getGFFtype(str(input)))
    shell:
        "mkdir -p offsets; echo {params.gfftype}; metagene generate -q offsets/metagene_stop --landmark cds_stop --annotation_files {input} --annotation_format {params.gfftype}"

rule psiteOffsetsStart:
    input:
        rois=rules.generateMetageneRoiStart.output,
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bamindex="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        psite="offsets/start/{method}-{condition}-{replicate}_p_offsets.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p offsets/start; psite -q {input.rois} offsets/start/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} --min_length 29 --max_length 35 --require_upstream --count_files {input.bam}"

rule psiteOffsetsStop:
    input:
        rois=rules.generateMetageneRoiStop.output,
        bam="maplink/{method}-{condition}-{replicate}.bam",
        bamindex="maplink/{method}-{condition}-{replicate}.bam.bai"
    output:
        psite="offsets/stop/{method}-{condition}-{replicate}_p_offsets.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p offsets/stop; psite -q {input.rois} offsets/stop/{wildcards.method}-{wildcards.condition}-{wildcards.replicate} --min_length 29 --max_length 35 --require_upstream --count_files {input.bam}"

