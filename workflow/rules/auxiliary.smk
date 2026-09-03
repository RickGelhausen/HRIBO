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


EXCEL_SCRIPT_DEPS = [
    str(SCRIPTS / "excel_utils.py"),
    str(SCRIPTS / "gff_utils.py"),
]


rule enrichAnnotation:
    input:
        annotation=rules.checkAnnotation.output,
        script=str(SCRIPTS / "enrich_annotation.py")
    output:
        "auxiliary/enriched_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -a {input.annotation:q} -o {output:q}"

rule unambigousAnnotation:
    input:
        "auxiliary/enriched_annotation.gff"
    output:
        "auxiliary/unambigous_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        awk -F'\\t' '/^[^#]/ {{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\tID=uid%s;\\n", $1, $2, $3, $4, $5, $6, $7, $8, NR-1}}' {input:q} > {output:q}
        """

rule samplesToExcel:
    input:
        samples=config["biologySettings"]["samples"],
        script=str(SCRIPTS / "samples_to_xlsx.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/samples.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -i {input.samples:q} -o {output:q}"

rule createExcelTotalAnnotation:
    input:
        total="readcounts/total_mapped_reads.txt",
        reads="readcounts/total_annotation.gtf",
        genome="genomes/genome.fa",
        script=str(SCRIPTS / "generate_excel.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/annotation_total.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -t {input.total:q} -r {input.reads:q} -g {input.genome:q} -o {output:q}"

rule createExcelUniqueAnnotation:
    input:
        total="readcounts/unique_mapped_reads.txt",
        reads="readcounts/unique_annotation.gtf",
        genome="genomes/genome.fa",
        script=str(SCRIPTS / "generate_excel.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/annotation_unique.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -t {input.total:q} -r {input.reads:q} -g {input.genome:q} -o {output:q}"

rule createExcelSummary:
    input:
        total="readcounts/bam_mapped_reads.txt",
        reads="readcounts/reparation_annotation.gff",
        genome="genomes/genome.fa",
        script=str(SCRIPTS / "generate_excel_reparation.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/predictions_reparation.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -t {input.total:q} -r {input.reads:q} -g {input.genome:q} -o {output:q}"

rule createExcelTotalAnnotationReadCount:
    input:
        reads="readcounts/total_annotation.gtf",
        total="readcounts/total_mapped_reads.txt",
        script=str(SCRIPTS / "generate_read_table.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/total_read_counts.xlsx",
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -r {input.reads:q} -t {input.total:q} -o {output:q}"

rule createExcelUniqueAnnotationReadCount:
    input:
        reads="readcounts/unique_annotation.gtf",
        total="readcounts/unique_mapped_reads.txt",
        script=str(SCRIPTS / "generate_read_table.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/unique_read_counts.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -r {input.reads:q} -t {input.total:q} -o {output:q}"

rule createOverviewTable:
    input:
        **overview_sources(),
        script=str(SCRIPTS / "generate_excel_overview.py"),
        script_deps=EXCEL_SCRIPT_DEPS
    output:
        "auxiliary/overview.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=120
    params:
        optional=overview_flags()
    log:
        "logs/overview_table.log"
    shell:
        """
        python3 {input.script:q} \
            -a {input.annotation:q} \
            -g {input.genome:q} \
            -t {input.totalreads:q} \
            {params.optional:q} \
            -o {output:q} 2> {log:q}
        """
