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
        featureCounts -F GTF -s 1 -g ID -t CDS -a {input.annotation} {input.bam} -T {threads} -o readcounts/reparation_read_counts.raw.tmp
        cat readcounts/reparation_read_counts.raw.tmp | sed 1,2d | awk '{{print $0"\\tCDS"}}' >> {output}
        rm readcounts/reparation_read_counts.raw.tmp
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
        featureCounts -F GTF -s 1 -g ID -t CDS -a {input.annotation} {input.bam} -T {threads} -o readcounts/deepribo_read_counts.raw.tmp
        cat readcounts/deepribo_read_counts.raw.tmp | sed 1,2d | awk -v var=CDS -FS'\\t' '{{print $0"\\t"var}}' >> {output}
        rm readcounts/deepribo_read_counts.raw.tmp
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
        UNIQUE="$(cut -f3 {input.annotation} | sort | uniq)"
        IDENTIFIER="ID"
        LINE="$(sed '3q;d' {input.annotation})"
        if [[ $LINE == *"gene_id="* ]]; then IDENTIFIER="gene_id"; fi;
        for f in ${{UNIQUE}}
        do
            featureCounts -F GTF -s 1 -g $IDENTIFIER -O -t $f -M --fraction -a {input.annotation} {input.bam} -T {threads} -o readcounts/annotation_total_reads.raw.tmp
            cat readcounts/annotation_total_reads.raw.tmp | sed 1,2d | awk -v var=$f -FS'\\t' '{{print $0"\\t"var}}' >> {output}
            rm readcounts/annotation_total_reads.raw.tmp
        done
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
        UNIQUE="$(cut -f3 {input.annotation} | sort | uniq)"
        IDENTIFIER="ID"
        LINE="$(sed '3q;d' {input.annotation})"
        if [[ $LINE == *"gene_id="* ]]; then IDENTIFIER="gene_id"; fi;
        for f in ${{UNIQUE}}
        do
            featureCounts -F GTF -s 1 -g $IDENTIFIER -O -t $f -M --fraction -a {input.annotation} {input.bam} -T {threads} -o readcounts/annotation_unique_reads.raw.tmp
            cat readcounts/annotation_unique_reads.raw.tmp | sed 1,2d | awk -v var=$f -FS'\\t' '{{print $0"\\t"var}}' >> {output}
            rm readcounts/annotation_unique_reads.raw.tmp
        done
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
        reads="readcounts/annotation_deepribo_reads.raw",
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
        annotation="auxiliary/enriched_annotation.gtf"
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
        annotation="auxiliary/enriched_annotation.gtf"
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

