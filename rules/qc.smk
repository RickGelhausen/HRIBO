rule fastqcunique:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/unique/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/unique/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/map/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/unique; fastqc -o qc/unique -t {threads} -f sam_mapped {input.sam}; mv qc/unique/{params.prefix}_fastqc.html {output.html}; mv qc/unique/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcmulti:
    input:
        sam="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/mapped/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/mapped/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/sammulti/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/mapped; fastqc -o qc/mapped -t {threads} -f sam_mapped {input.sam}; mv qc/mapped/{params.prefix}_fastqc.html {output.html}; mv qc/mapped/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcraw:
    input:
        reads=getfastq,
    output:
        html="qc/raw/{method}-{condition}-{replicate}-raw_fastqc.html",
        zip="qc/raw/{method}-{condition}-{replicate}-raw_fastqc.zip"
        #report("qc/raw/{method}-{condition}-{replicate}-raw.html", caption="../report/fastqcraw.rst", category="Input quality control")
    conda:
        "../envs/fastqc.yaml"
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.reads[0]))[0])[0])
    threads: 8
    shell:
        "mkdir -p qc/raw; fastqc -o qc/raw -t {threads} {input}; mv qc/raw/{params.prefix}_fastqc.html {output.html}; mv qc/raw/{params.prefix}_fastqc.zip {output.zip}"

rule fastqctrimmed:
    input:
        reads="trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        html="qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html",
        zip="qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.zip"
        #report("qc/trimmed/{method}-{condition}-{replicate}-trimmed.html", caption="../report/fastqctrimmed.rst", category="Trimming")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/trimmed; fastqc -o qc/trimmed -t {threads} {input}; mv qc/trimmed/{params.prefix}_fastqc.html {output.html}; mv qc/trimmed/{params.prefix}_fastqc.zip {output.zip}"

rule fastqcrrnafilter:
    input:
        reads="bam/{method}-{condition}-{replicate}.bam"
    output:
        html="qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html",
        zip="qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.zip"
        #report("qc/norRNA/{method}-{condition}-{replicate}-norRNA.html", caption="../report/fastqcnorRNA.rst", category="Removing hits mapping to rRNA")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/norRNA; fastqc -o qc/norRNA -t {threads} {input}; mv qc/norRNA/{params.prefix}_fastqc.html {output.html}; mv qc/norRNA/{params.prefix}_fastqc.zip {output.zip}"

rule gff2gtf:
    input:
        annotation={rules.retrieveAnnotation.output},
    output:
        gtfcds="qc/featurecount/annotation.gtf",
    conda:
        "../envs/cufflinks.yaml"
    threads: 1
    shell:
        "mkdir -p qc/featurecount; gffread -T {input.annotation} -o {output.gtfcds};"

rule featurecountAnnotation:
    input:
        annotation={rules.retrieveAnnotation.output},
    output:
        "qc/featurecount/annotationall.gtf",
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p qc/all; SPtools/scripts/annotation_featurecount.py -a {input.annotation} -o {output};"

rule featurescounts:
    input:
        annotation={rules.featurecountAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/all/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/all; featureCounts -T {threads} -t gene -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"

rule trnafeaturescounts:
    input:
        annotation={rules.featurecountAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/trnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/trnainall; featureCounts -T {threads} -t tRNA -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"

rule norrnafeaturescounts:
    input:
        annotation={rules.featurecountAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/rrnainall; featureCounts -T {threads} -t rRNA -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"

rule rrnatotalfeaturescounts:
    input:
        annotation={rules.featurecountAnnotation.output},
        bam="bammulti/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainallaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/rrnainallaligned; featureCounts -T {threads} -t rRNA -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"

rule rrnauniquefeaturescounts:
    input:
        annotation={rules.featurecountAnnotation.output},
        bam="rRNAbam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/rrnainuniquelyaligned; featureCounts -T {threads} -t rRNA -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"


#rule ncrnafeaturescounts:
#   input:
#       annotation={rules.featurecountAnnotation.output},
#       bam="bam/{method}-{condition}-{replicate}.bam"
#   output:
#       txt="qc/ncrnafeaturecount/{method}-{condition}-{replicate}.txt",
#   conda:
#       "../envs/subread.yaml"
#   threads: 8
#   shell:
#       "mkdir -p qc/ncrnarnafeaturecount; featureCounts -T {threads} -t ncRNA -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"


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

rule multiqc:
    input:
        expand("tracks/{condition}.ribotish.gff", condition=set(samples["condition"])),
        expand("qc/raw/{method}-{condition}-{replicate}-raw_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/unique/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/mapped/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/all/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainallaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
 #       expand("qc/ncrnafeaturecount/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        report("qc/multi/multiqc_report.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        dir="qc/multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "export LC_ALL=en_US.utf8; export LANG=en_US.utf8; multiqc -f -d --exclude picard --exclude gatk -z -o {params.dir} qc/mapped qc/raw qc/trimmed qc/norRNA qc/unique qc/featurecount qc/trnafeaturecount qc/rrnatotalfeaturecount qc/rrnauniquefeaturecount qc/norrnafeaturecount trimmed  2> {log}"
