# Read counting
#
# Every featureCounts run goes through one rule; the runs differ only in which
# alignments they count, which annotation they count against and a few extra
# flags. Those differences live in READ_COUNT_SETS in common.smk, keyed by the
# output basename, which is also the {countset} wildcard.

wildcard_constraints:
    countset="|".join(re.escape(name) for name in READ_COUNT_SETS),
    mapped="|".join(re.escape(name) for name in READ_COUNT_ANNOTATIONS),
    source="|".join(MAPPED_READ_SOURCES)


rule readCounts:
    input:
        bam=readcount_bams,
        bamindex=readcount_bam_indices,
        annotation=readcount_annotation,
        script=str(SCRIPTS / "call_featurecounts.py")
    output:
        "readcounts/{countset}"
    conda:
        "../envs/subread.yaml"
    threads: 5
    resources:
        mem_mb=30000,
        runtime=120
    params:
        extra=readcount_flags
    log:
        "logs/readcounts_{countset}.log"
    shell:
        """
        python3 {input.script:q} \
            -b {input.bam:q} \
            -a {input.annotation:q} \
            -s 1 --with_O {params.extra:q} \
            -t {threads} \
            -o {output:q} 2> {log:q}
        awk '/^WARNING:/' {log:q} >&2
        """


rule mapReadsToAnnotation:
    input:
        reads=mapped_counts_reads,
        annotation=mapped_counts_annotation,
        script=str(SCRIPTS / "map_reads_to_annotation.py")
    output:
        "readcounts/{mapped}"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "logs/map_reads_to_annotation_{mapped}.log"
    shell:
        "python3 {input.script:q} -i {input.reads:q} -a {input.annotation:q} -o {output:q} 2> {log:q}"


rule mappedReadSummary:
    input:
        bam=mapped_read_bams,
        bamindex=mapped_read_bam_indices,
        script=str(SCRIPTS / "total_mapped_reads.py")
    output:
        mapped="readcounts/{source}_mapped_reads.txt",
        length="readcounts/{source}_average_read_lengths.txt"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "logs/mapped_read_summary_{source}.log"
    shell:
        "python3 {input.script:q} -b {input.bam:q} -m {output.mapped:q} -l {output.length:q} 2> {log:q}"
