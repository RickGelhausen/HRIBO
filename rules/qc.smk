rule fastqcmapped:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/map/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/map/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/map/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "mkdir -p qc/map; fastqc -o qc/map -t {threads} -f sam_mapped {input.sam}; mv qc/map/{params.prefix}_fastqc.html {output.html}; mv qc/map/{params.prefix}_fastqc.zip {output.zip}"

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
        reads="norRNA/{method}-{condition}-{replicate}.fastq"
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

rule featurescounts:
    input:
        annotation={rules.retrieveAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/featurecount/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/featurecount; featureCounts -T {threads} -t exon -g gene_id -a {input.annotation} -o {output.txt} {input.bam}"

rule multiqc:
    input:
        expand("tracks/{sample.condition}.ribotish.gff", sample=samples.itertuples()),
        expand("tracks/{sample.method}-{sample.condition}-{sample.replicate}.bw", sample=samples.itertuples()),
        expand("qc/raw/{sample.method}-{sample.condition}-{sample.replicate}-raw_fastqc.html", sample=samples.itertuples()),
        expand("qc/trimmed/{sample.method}-{sample.condition}-{sample.replicate}-trimmed_fastqc.html", sample=samples.itertuples()),
        expand("qc/norRNA/{sample.method}-{sample.condition}-{sample.replicate}-norRNA_fastqc.html", sample=samples.itertuples()),
        expand("qc/map/{sample.method}-{sample.condition}-{sample.replicate}-map_fastqc.html", sample=samples.itertuples()),
        expand("qc/featurecount/{sample.method}-{sample.condition}-{sample.replicate}.txt", sample=samples.itertuples())
        #expand("bam/{sample.method}-{sample.condition}-{sample.replicate}.bam", sample=samples.itertuples())

    output:
        report("qc/multi/multiqc_report.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        dir="qc/multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "export LC_ALL=en_US.utf8; export LANG=en_US.utf8; multiqc -f -d --exclude picard --exclude gatk -z -o {params.dir} qc/raw qc/trimmed qc/norRNA qc/map qc/featurecount trimmed  2> {log}"
