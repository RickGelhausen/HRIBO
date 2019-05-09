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

rule gff2gtf:
    input:
        annotation={rules.retrieveAnnotation.output},
    output:
        gtf="qc/featurecount/annotation.gtf",
    conda:
        "../envs/cufflinks.yaml"
    threads: 1
    shell:
        "mkdir -p qc/featurecount; gffread {input.annotation} -T -o {output.gtf}"

rule featurescounts:
    input:
        annotation={rules.gff2gtf.output.gtf},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/featurecount/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    shell:
        "mkdir -p qc/featurecount; featureCounts -T {threads} -t exon -g transcript_id -a {input.annotation} -o {output.txt} {input.bam}"

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
        expand("tracks/{method}-{condition}-{replicate}.fwd.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("tracks/{method}-{condition}-{replicate}.rev.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/raw/{method}-{condition}-{replicate}-raw_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/norRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/map/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/featurecount/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/min/{method}-{condition}-{replicate}.min.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("globaltracks/min/{method}-{condition}-{replicate}.min.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/min/{method}-{condition}-{replicate}.min.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])

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
