rule fastqcunique:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/4unique/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/4unique/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/map/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/4unique; fastqc -o qc/4unique -t {threads} -f sam_mapped {input.sam}; mv qc/4unique/{params.prefix}_fastqc.html {output.html}; mv qc/4unique/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcmulti:
    input:
        sam="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/sammulti/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/3mapped; fastqc -o qc/3mapped -t {threads} -f sam_mapped {input.sam}; mv qc/3mapped/{params.prefix}_fastqc.html {output.html}; mv qc/3mapped/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcrrnafilter:
    input:
        reads="bam/{method}-{condition}-{replicate}.bam"
    output:
        html="qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html",
        zip="qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.zip"
        #report("qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA.html", caption="../report/fastqcnorRNA.rst", category="Removing hits mapping to rRNA")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/5removedrRNA; fastqc -o qc/5removedrRNA -t {threads} {input}; mv qc/5removedrRNA/{params.prefix}_fastqc.html {output.html}; mv qc/5removedrRNA/{params.prefix}_fastqc.zip {output.zip}"

rule featurescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/all/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/all; featureCounts -T {threads} -t gene -g ID -a {input.annotation} -o {output.txt} {input.bam}"

rule trnafeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/trnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/trnainall; featureCounts -T {threads} -t tRNA -g ID -a {input.annotation} -o {output.txt} {input.bam}"

rule norrnafeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/rrnainall; featureCounts -T {threads} -t rRNA -g ID -a {input.annotation} -o {output.txt} {input.bam}"

rule rrnatotalfeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bammulti/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainallaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/rrnainallaligned; featureCounts -T {threads} -t rRNA -g ID -a {input.annotation} -o {output.txt} {input.bam}"

rule rrnauniquefeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="rRNAbam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/rrnainuniquelyaligned; featureCounts -T {threads} -t rRNA -g ID -a {input.annotation} -o {output.txt} {input.bam}"

rule coveragedepth:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "coverage/{method}-{condition}-{replicate}.bed"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p coverage; bedtools genomecov -ibam {input} -bg > {output}"
