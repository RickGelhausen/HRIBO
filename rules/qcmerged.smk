
rule fastqcraw:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_Q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_P.fastq.gz"
    output:
        html1="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc_Q.html",
        zip1="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc_Q.zip",
        html2="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc_P.html",
        zip2="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc_P.zip"
        #report("qc/1raw/{method}-{condition}-{replicate}-raw.html", caption="../report/fastqcraw.rst", category="Input quality control")
    conda:
        "../envs/fastqc.yaml"
    params:
        prefix1=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.fastq1))[0])[0]),
        prefix2=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.fastq2))[0])[0])
    threads: 8
    shell:
        """
        mkdir -p qc/1raw;
        fastqc -o qc/1raw -t {threads} {input.fastq1}; mv qc/1raw/{params.prefix1}_fastqc.html {output.html1}; mv qc/1raw/{params.prefix1}_fastqc.zip {output.zip1}
        fastqc -o qc/1raw -t {threads} {input.fastq2}; mv qc/1raw/{params.prefix2}_fastqc.html {output.html2}; mv qc/1raw/{params.prefix2}_fastqc.zip {output.zip2}
        """


rule fastqctrimmed:
    input:
        reads="trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        html="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html",
        zip="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.zip"
        #report("qc/2trimmed/{method}-{condition}-{replicate}-trimmed.html", caption="../report/fastqctrimmed.rst", category="Trimming")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/2trimmed; fastqc -o qc/2trimmed -t {threads} {input}; mv qc/2trimmed/{params.prefix}_fastqc.html {output.html}; mv qc/2trimmed/{params.prefix}_fastqc.zip {output.zip}"


rule multiqc:
    input:
        expand("qc/1raw/{method}-{condition}-{replicate}-raw_fastqc_Q.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/1raw/{method}-{condition}-{replicate}-raw_fastqc_P.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/4unique/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/all/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainallaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
 #       expand("qc/ncrnafeaturecount/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        "qc/multi/multiqc_report.html"
    params:
        dir="qc/multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "export LC_ALL=en_US.utf8; export LANG=en_US.utf8; multiqc -f -d --exclude picard --exclude gatk -z -o {params.dir} qc/3mapped qc/1raw qc/2trimmed qc/5removedrRNA qc/4unique qc/all qc/trnainall qc/rrnainallaligned qc/rrnainuniquelyaligned qc/rrnainall trimmed  2> {log}"
