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
        rules.ribotishAnnotation.output
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
        rules.ribotishAnnotation.output
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

rule generateReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        combined="tracks/combined_annotated.gff"
    output:
        "auxiliary/read_counts.bed"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary
        cut -f1,4,5,7,9 {input.combined} > tmp.bed
        bedtools multicov -bams {input.bam} -bed tmp.bed > {output}
        sed -i '1i \# {input.bam}\n' {output}
        rm tmp.bed
        """

rule totalMappedReads:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        "auxiliary/total_mapped_reads.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/tmp.py -b {input.bam} -o {output}"

rule createExcelSummary:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/read_counts.bed"
    output:
        report("auxiliary/summary.xlsx", caption="../report/summary.rst", category="Summary table")
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/excel.py -t {input.total} -r {input.reads} -o {output}"