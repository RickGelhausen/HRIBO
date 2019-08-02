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

rule enrichAnnotation:
    input:
        "annotation/annotation.gtf"
    output:
        "auxiliary/enriched_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/enrich_annotation.py -a {input} -o {output}"

rule samplesToExcel:
    input:
        "SPtools/samples.tsv"
    output:
        "auxiliary/samples.xlsx"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/samples_to_xlsx.py -i {input} -o {output}"
 
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
        awk -F'\\t' '{{ print $1 FS $4 FS $5 FS $9 FS $6 FS $7 FS $2 FS $3}}' {input.combined} > tmp_combined.bed
        bedtools multicov -s -D -bams {input.bam} -bed tmp_combined.bed > {output}
        sed -i '1i \# {input.bam}\n' {output}
        rm tmp_combined.bed
        """

rule generateAnnotationTotalReadCounts:
    input:
        bam=expand("bammulti/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("bammulti/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/annotation_total_read_counts.bed"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary
        awk -F'\t' '{{ print $1 FS $4 FS $5 FS $9 FS $6 FS $7 FS $2 FS $3}}' {input.annotation} > tmp_annotation_total.bed
        bedtools multicov -s -D -bams {input.bam} -bed tmp_annotation_total.bed > {output}
        sed -i '1i \# {input.bam}\n' {output}
        rm tmp_annotation_total.bed
        """

rule generateAnnotationUniqueReadCounts:
    input:
        bam=expand("rRNAbam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("rRNAbam/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/annotation_unique_read_counts.bed"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary
        awk -F'\t' '{{ print $1 FS $4 FS $5 FS $9 FS $6 FS $7 FS $2 FS $3}}' {input.annotation} > tmp_annotation_unique.bed
        bedtools multicov -s -D -bams {input.bam} -bed tmp_annotation_unique.bed > {output}
        sed -i '1i \# {input.bam}\n' {output}
        rm tmp_annotation_unique.bed
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

rule createExcelTotalAnnotation:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/annotation_total_read_counts.bed",
        genome="genomes/genome.fa"
    output:
        "auxiliary/annotation_total.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelUniqueAnnotation:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/annotation_unique_read_counts.bed",
        genome="genomes/genome.fa"
    output:
        "auxiliary/annotation_unique.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelTotalAnnotationReadCount:
    input:
        reads="auxiliary/annotation_total_read_counts.bed"
    output:
        "auxiliary/total_read_counts.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/generate_read_table.py -r {input.reads} -o {output}"

rule createExcelUniqueAnnotationReadCount:
    input:
        reads="auxiliary/annotation_unique_read_counts.bed"
    output:
        "auxiliary/unique_read_counts.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/generate_read_table.py -r {input.reads} -o {output}"

rule createExcelSummary:
    input:
        total="auxiliary/total_mapped_reads.txt",
        reads="auxiliary/combined_read_counts.bed",
        genome="genomes/genome.fa"
    output:
        report("auxiliary/summary.xlsx", caption="../report/summary.rst", category="Summary table")
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; SPtools/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"
