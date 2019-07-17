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

rule generateCombinedReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        combined="tracks/combined_annotated.gff"
    output:
        "auxiliary/combined_read_counts.bed"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary
        awk -F'\\t' '{{ print $1 FS $4 FS $5 FS $9 FS $6 FS $7 }}' {input.combined} > tmp_combined.bed
        bedtools multicov -s -D -bams {input.bam} -bed tmp_combined.bed > {output}
        sed -i '1i \# {input.bam}\n' {output}
        rm tmp_combined.bed
        """

rule generateUniqAnnotation:
    input:
        "annotation/annotation.gtf"
    output:
        "auxiliary/annotation_uniq.gtf"
    threads: 1
    shell:
        "awk -F'\t' '!seen[$1 FS $4 FS $5 FS $7]++' {input} > {output}"

rule generateAnnotationReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/annotation_uniq.gtf"
    output:
        "auxiliary/annotation_read_counts.bed"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary
        awk -F'\t' '{{ print $1 FS $4 FS $5 FS $9 FS $6 FS $7 FS $2 FS $3}}' {input.annotation} > tmp_annotation.bed
        bedtools multicov -s -D -bams {input.bam} -bed tmp_annotation.bed > {output}
        sed -i '1i \# {input.bam}\n' {output}
        rm tmp_annotation.bed
        """

rule annotationToBed:
    input:
        annotation="auxiliary/annotation_uniq.gtf"
    output:
        "auxiliary/annotation_uniq.bed"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary
        cut -f1,4,5,7 {input.annotation} > {output}
        """

rule totalMappedReads:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="auxiliary/total_mapped_reads.txt",
        length="auxiliary/average_read_lengths.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

# rule calculateAverageLengths:
#     input:
#         bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
#         bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
#         annotation="auxiliary/annotation_uniq.bed"
#     output:
#         length="auxiliary/average_read_lengths.bed",
#     conda:
#         "../envs/plastid.yaml"
#     threads: 1
#     shell:
#         "mkdir -p auxiliary; SPtools/scripts/average_read_lengths.py -b {input.bam} -a {input.annotation} -l {output.length}"

rule createExcelAnnotation:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/annotation_read_counts.bed",
        length="auxiliary/average_read_lengths.txt"
    output:
        rpkm= "auxiliary/annotation_rpkm.xlsx",
        tpm= "auxiliary/annotation_tpm.xlsx"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/measures_excel.py -t {input.total} -r {input.reads} -l {input.length} --out_tpm {output.tpm} --out_rpkm {output.rpkm}"

rule createExcelSummary:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/combined_read_counts.bed",
        genome="genomes/genome.fa"
    output:
        report("auxiliary/summary.xlsx", caption="../report/summary.rst", category="Summary table")
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/summary_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"
