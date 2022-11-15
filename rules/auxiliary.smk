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
        annotation=rules.checkAnnotation.output
    output:
        "auxiliary/enriched_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/enrich_annotation.py -a {input.annotation} -o {output}"

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
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel_reparation.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

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

rule createOverviewTableReparation:
    input:
        annotation="readcounts/independant_annotation.gff",
        genome=rules.retrieveGenome.output,
        totalreads="readcounts/bam_mapped_reads.txt",
        reparation="readcounts/reparation_annotation.gff"
    output:
        "auxiliary/overview.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        contrasts=CONTRASTS
    shell:
        """
        mkdir -p auxiliary;
        if [ -z {params.contrasts} ]
        then
            HRIBO/scripts/generate_excel_overview.py -a {input.annotation} -g {input.genome} -t {input.totalreads} --mapped_reads_reparation {input.reparation} -o {output}
        else
            HRIBO/scripts/generate_excel_overview.py -a {input.annotation} -c {params.contrasts} -g {input.genome} -t {input.totalreads} --mapped_reads_reparation {input.reparation} -o {output}
        fi
        """

rule createOverviewTablePredictions:
    input:
        annotation="readcounts/independant_annotation.gff",
        genome=rules.retrieveGenome.output,
        totalreads="readcounts/bam_mapped_reads.txt",
        reparation="readcounts/reparation_annotation.gff",
        deepribo="readcounts/deepribo_annotation.gff"
    output:
        "auxiliary/overview.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        contrasts=CONTRASTS
    shell:
        """
        mkdir -p auxiliary;
        if [ -z {params.contrasts} ]
        then
            HRIBO/scripts/generate_excel_overview.py -a {input.annotation} -g {input.genome} -t {input.totalreads} --mapped_reads_deepribo {input.deepribo} --mapped_reads_reparation {input.reparation} -o {output}
        else
            HRIBO/scripts/generate_excel_overview.py -a {input.annotation} -c {params.contrasts}  -g {input.genome} -t {input.totalreads} --mapped_reads_deepribo {input.deepribo} --mapped_reads_reparation {input.reparation} -o {output}
        fi
        """

rule createOverviewTableDiffExpr:
    input:
        annotation="readcounts/independant_annotation.gff",
        genome=rules.retrieveGenome.output,
        xtail="xtail/xtail_all.csv",
        riborex="riborex/riborex_all.csv",
        deltate="deltate/deltate_all.csv",
        totalreads="readcounts/bam_mapped_reads.txt",
        reparation="readcounts/reparation_annotation.gff"
    output:
        "auxiliary/overview.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        contrasts=CONTRASTS
    shell:
        """
        if [ -z {params.contrasts} ]
        then
            mkdir -p auxiliary; HRIBO/scripts/generate_excel_overview.py -a {input.annotation} -g {input.genome} --xtail {input.xtail} --deltate {input.deltate} --riborex {input.riborex} -t {input.totalreads} --mapped_reads_reparation {input.reparation} -o {output}
        else
            mkdir -p auxiliary; HRIBO/scripts/generate_excel_overview.py -c {params.contrasts} -a {input.annotation} -g {input.genome} --xtail {input.xtail} --deltate {input.deltate} --riborex {input.riborex} -t {input.totalreads} --mapped_reads_reparation {input.reparation} -o {output}
        fi
        """

rule createOverviewTableAll:
    input:
        annotation="readcounts/independant_annotation.gff",
        genome=rules.retrieveGenome.output,
        xtail="xtail/xtail_all.csv",
        riborex="riborex/riborex_all.csv",
        deltate="deltate/deltate_all.csv",
        totalreads="readcounts/bam_mapped_reads.txt",
        reparation="readcounts/reparation_annotation.gff",
        deepribo="readcounts/deepribo_annotation.gff"
    output:
        "auxiliary/overview.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        contrasts=CONTRASTS
    shell:
        """
        if [ -z {params.contrasts} ]
        then
            mkdir -p auxiliary; HRIBO/scripts/generate_excel_overview.py -a {input.annotation} -g {input.genome} --xtail {input.xtail} --deltate {input.deltate} --riborex {input.riborex} -t {input.totalreads} --mapped_reads_deepribo {input.deepribo} --mapped_reads_reparation {input.reparation} -o {output}
        else
            mkdir -p auxiliary; HRIBO/scripts/generate_excel_overview.py -c {params.contrasts} -a {input.annotation} -g {input.genome} --xtail {input.xtail} --deltate {input.deltate} --riborex {input.riborex} -t {input.totalreads} --mapped_reads_deepribo {input.deepribo} --mapped_reads_reparation {input.reparation} -o {output}
        fi
        """
