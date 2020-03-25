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

rule enrichAnnotation:
    input:
        annotation=rules.retrieveAnnotation.output
    output:
        "auxiliary/enriched_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/enrich_annotation.py -a {input.annotation} -o {output}"

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

rule createExcelTotalAnnotation:
    input:
        total="readcounts/total_mapped_reads.txt",
        reads="readcounts/total_annotation.gtf",
        genome="genomes/genome.fa"
    output:
        "auxiliary/annotation_total.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelUniqueAnnotation:
    input:
        total="readcounts/unique_mapped_reads.txt",
        reads="readcounts/unique_annotation.gtf",
        genome="genomes/genome.fa"
    output:
        "auxiliary/annotation_unique.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelSummary:
    input:
        total="readcounts/bam_mapped_reads.txt",
        reads="readcounts/reparation_annotation.gff",
        genome="genomes/genome.fa"
    output:
        "auxiliary/predictions_reparation.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule createExcelTotalAnnotationReadCount:
    input:
        reads="readcounts/total_annotation.gtf",
        total="readcounts/total_mapped_reads.txt"
    output:
        "auxiliary/total_read_counts.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_read_table.py -r {input.reads} -t {input.total} -o {output}"

rule createExcelUniqueAnnotationReadCount:
    input:
        reads="readcounts/unique_annotation.gtf",
        total="readcounts/unique_mapped_reads.txt"
    output:
        "auxiliary/unique_read_counts.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_read_table.py -r {input.reads} -t {input.total} -o {output}"
