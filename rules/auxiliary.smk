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
        "mkdir -p offsets; HRIBO/scripts/process_annotation.py -a {input.annotation} -o {output}"

rule enrichAnnotation:
    input:
        "annotation/annotation.gtf"
    output:
        "auxiliary/enriched_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/enrich_annotation.py -a {input} -o {output}"

rule unambigousAnnotation:
    input:
        "auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/unambigous_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary;
        awk -F'\\t' '/^[^#]/ {{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\tID=uid%s;\\n", $1, $2, $3, $4, $5, $6, $7, $8, NR-1}}' {input} > {output}
        """

rule samplesToExcel:
    input:
        "HRIBO/samples.tsv"
    output:
        "auxiliary/samples.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/samples_to_xlsx.py -i {input} -o {output}"

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
        annotation="tracks/combined_annotated.gff"
    output:
        "auxiliary/combined_read_counts.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p auxiliary
        featureCounts -F GTF -s 1 -g ID -O -t CDS -M --fraction -a {input.annotation} {input.bam} -T {threads} -o auxiliary/combined_read_counts.raw.tmp
        cat auxiliary/combined_read_counts.raw.tmp | sed 1,2d | awk '{{print $0"\\tCDS"}}' >> {output}
        rm auxiliary/combined_read_counts.raw.tmp
        """

rule generateAnnotationTotalReadCounts:
    input:
        bam=expand("bammulti/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("bammulti/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/unambigous_annotation.gtf"
    output:
        "auxiliary/annotation_total_reads.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p auxiliary
        UNIQUE="$(cut -f3 {input.annotation} | sort | uniq)"
        IDENTIFIER="ID"
        LINE="$(sed '3q;d' {input.annotation})"
        if [[ $LINE == *"gene_id="* ]]; then IDENTIFIER="gene_id"; fi;
        for f in ${{UNIQUE}}
        do
            featureCounts -F GTF -s 1 -g $IDENTIFIER -O -t $f -M --fraction -a {input.annotation} {input.bam} -T {threads} -o auxiliary/annotation_total_reads.raw.tmp
            cat auxiliary/annotation_total_reads.raw.tmp | sed 1,2d | awk -v var=$f -FS'\\t' '{{print $0"\\t"var}}' >> {output}
            rm auxiliary/annotation_total_reads.raw.tmp
        done
        """

rule generateAnnotationUniqueReadCounts:
    input:
        bam=expand("rRNAbam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("rRNAbam/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/unambigous_annotation.gtf"
    output:
        "auxiliary/annotation_unique_reads.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p auxiliary
        UNIQUE="$(cut -f3 {input.annotation} | sort | uniq)"
        IDENTIFIER="ID"
        LINE="$(sed '3q;d' {input.annotation})"
        if [[ $LINE == *"gene_id="* ]]; then IDENTIFIER="gene_id"; fi;
        for f in ${{UNIQUE}}
        do
            featureCounts -F GTF -s 1 -g $IDENTIFIER -O -t $f -M --fraction -a {input.annotation} {input.bam} -T {threads} -o auxiliary/annotation_unique_reads.raw.tmp
            cat auxiliary/annotation_unique_reads.raw.tmp | sed 1,2d | awk -v var=$f -FS'\\t' '{{print $0"\\t"var}}' >> {output}
            rm auxiliary/annotation_unique_reads.raw.tmp
        done
        """

rule mapPredictionReads:
    input:
        reads="auxiliary/combined_read_counts.raw",
        annotation="tracks/combined_annotated.gff"
    output:
        "auxiliary/prediction_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mapTotalReads:
    input:
        reads="auxiliary/annotation_total_reads.raw",
        annotation="auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/total_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mapUniqueReads:
    input:
        reads="auxiliary/annotation_unique_reads.raw",
        annotation="auxiliary/enriched_annotation.gtf"
    output:
        "auxiliary/unique_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule totalMappedReads:
    input:
        bam=expand("bammulti/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("bammulti/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="auxiliary/total_sum_mapped_reads.txt",
        length="auxiliary/total_average_read_lengths.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

rule uniqueMappedReads:
    input:
        bam=expand("rRNAbam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("rRNAbam/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="auxiliary/unique_sum_mapped_reads.txt",
        length="auxiliary/unique_average_read_lengths.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

rule finalMappedReads:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="auxiliary/final_sum_mapped_reads.txt",
        length="auxiliary/final_average_read_lengths.txt"
    conda:
        "../envs/plastid.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

rule createExcelTotalAnnotation:
    input:
        total="auxiliary/total_sum_mapped_reads.txt",
        reads="auxiliary/total_annotation.gtf",
        genome="genomes/genome.fa"
    output:
        "auxiliary/annotation_total.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelUniqueAnnotation:
    input:
        total="auxiliary/unique_sum_mapped_reads.txt",
        reads="auxiliary/unique_annotation.gtf",
        genome="genomes/genome.fa"
    output:
        "auxiliary/annotation_unique.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelSummary:
    input:
        total="auxiliary/final_sum_mapped_reads.txt",
        reads="auxiliary/prediction_annotation.gtf",
        genome="genomes/genome.fa"
    output:
        report("auxiliary/summary.xlsx", caption="../report/summary.rst", category="Summary table")
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelTotalAnnotationReadCount:
    input:
        reads="auxiliary/total_annotation.gtf",
        total="auxiliary/total_sum_mapped_reads.txt"
    output:
        "auxiliary/total_read_counts.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_read_table.py -r {input.reads} -t {input.total} -o {output}"

rule createExcelUniqueAnnotationReadCount:
    input:
        reads="auxiliary/unique_annotation.gtf",
        total="auxiliary/unique_sum_mapped_reads.txt"
    output:
        "auxiliary/unique_read_counts.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_read_table.py -r {input.reads} -t {input.total} -o {output}"
