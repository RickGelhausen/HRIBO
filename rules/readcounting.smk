rule generateDifferentialExpressionReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="tracks/updated_annotation.gff"
    output:
        "readcounts/differential_expression_read_counts.csv"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p readcounts
        HRIBO/scripts/call_featurecounts.py -b {input.bam} -s 1 --with_O --for_diff_expr -o {output} -t {threads} -a {input.annotation}
        """

rule generateReparationReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="tracks/reparation_annotated.gff"
    output:
        "readcounts/reparation_read_counts.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p readcounts
        HRIBO/scripts/call_featurecounts.py -b {input.bam} -s 1 --with_O -o {output} -t {threads} -a {input.annotation}
        """

rule generateDeepRiboReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="tracks/deepribo_merged.gff"
    output:
        "readcounts/deepribo_read_counts.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p readcounts
        HRIBO/scripts/call_featurecounts.py -b {input.bam} -s 1 --with_O -o {output} -t {threads} -a {input.annotation}
        """

rule generateAnnotationIndependantReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/unambigous_annotation.gtf"
    output:
        "readcounts/annotation_independant_read_counts.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p readcounts
        HRIBO/scripts/call_featurecounts.py -b {input.bam} -s 1 --with_O -o {output} -t {threads} -a {input.annotation}
        """

rule generateAnnotationTotalReadCounts:
    input:
        bam=expand("bammulti/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("bammulti/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/unambigous_annotation.gtf"
    output:
        "readcounts/annotation_total_reads.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p readcounts
        HRIBO/scripts/call_featurecounts.py -b {input.bam} -s 1 --with_O --with_M --fraction -o {output} -t {threads} -a {input.annotation}
        """

rule generateAnnotationUniqueReadCounts:
    input:
        bam=expand("rRNAbam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("rRNAbam/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="auxiliary/unambigous_annotation.gtf"
    output:
        "readcounts/annotation_unique_reads.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p auxiliary
        HRIBO/scripts/call_featurecounts.py -b {input.bam} -s 1 --with_O --fraction -o {output} -t {threads} -a {input.annotation}
        """

rule mapIndependantReads:
    input:
        reads="readcounts/annotation_independant_read_counts.raw",
        annotation=rules.checkAnnotation.output
    output:
        "readcounts/independant_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p readcounts; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mapReparationReads:
    input:
        reads="readcounts/reparation_read_counts.raw",
        annotation="tracks/reparation_annotated.gff"
    output:
        "readcounts/reparation_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p readcounts; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mapDeepRiboReads:
    input:
        reads="readcounts/deepribo_read_counts.raw",
        annotation="tracks/deepribo_merged.gff"
    output:
        "readcounts/deepribo_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p readcounts; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mapTotalReads:
    input:
        reads="readcounts/annotation_total_reads.raw",
        annotation="auxiliary/enriched_annotation.gff"
    output:
        "readcounts/total_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p readcounts; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mapUniqueReads:
    input:
        reads="readcounts/annotation_unique_reads.raw",
        annotation="auxiliary/enriched_annotation.gff"
    output:
        "readcounts/unique_annotation.gtf"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p readcounts; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule totalMappedReads:
    input:
        bam=expand("bammulti/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("bammulti/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="readcounts/total_mapped_reads.txt",
        length="readcounts/total_average_read_lengths.txt"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p readcounts; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

rule uniqueMappedReads:
    input:
        bam=expand("rRNAbam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("rRNAbam/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="readcounts/unique_mapped_reads.txt",
        length="readcounts/unique_average_read_lengths.txt"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p readcounts; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"

rule maplinkMappedReads:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="readcounts/bam_mapped_reads.txt",
        length="readcounts/bam_average_read_lengths.txt"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p readcounts; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"
