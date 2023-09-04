
rule fastqcraw_single:
    input:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
    output:
        html="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc.html",
        zip="qc/1raw/{method}-{condition}-{replicate}-raw_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.fastq))[0])[0])
    threads: 8
    shell:
        "mkdir -p qc/1raw; fastqc -o qc/1raw -t {threads} {input.fastq}; mv qc/1raw/{params.prefix}_fastqc.html {output.html}; mv qc/1raw/{params.prefix}_fastqc.zip {output.zip}"

rule fastqctrimmed_single:
    input:
        reads="trimmed/{method}-{condition}-{replicate}.fastq"
    output:
        html="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html",
        zip="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "mkdir -p qc/2trimmed; fastqc -o qc/2trimmed -t {threads} {input}; mv qc/2trimmed/{params.prefix}_fastqc.html {output.html}; mv qc/2trimmed/{params.prefix}_fastqc.zip {output.zip}"


rule fastqcraw_paired:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_p.fastq.gz"
    output:
        html1="qc/1raw/{method}-{condition}-{replicate}-raw-q_fastqc.html",
        zip1="qc/1raw/{method}-{condition}-{replicate}-raw-q_fastqc.zip",
        html2="qc/1raw/{method}-{condition}-{replicate}-raw-p_fastqc.html",
        zip2="qc/1raw/{method}-{condition}-{replicate}-raw-p_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    params:
        prefix1=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.fastq1))[0])[0]),
        prefix2=lambda wildcards, input: (os.path.splitext(os.path.splitext(os.path.basename(input.fastq2))[0])[0])
    threads: 8
    shell:
        """
        mkdir -p qc/1raw
        fastqc -o qc/1raw -t {threads} {input.fastq1}; mv qc/1raw/{params.prefix1}_fastqc.html {output.html1}; mv qc/1raw/{params.prefix1}_fastqc.zip {output.zip1}
        fastqc -o qc/1raw -t {threads} {input.fastq2}; mv qc/1raw/{params.prefix2}_fastqc.html {output.html2}; mv qc/1raw/{params.prefix2}_fastqc.zip {output.zip2}
        """

rule fastqctrimmed_paired:
    input:
        reads1="trimmed/{method}-{condition}-{replicate}_q.fastq",
        reads2="trimmed/{method}-{condition}-{replicate}_p.fastq"
    output:
        html1="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_q_fastqc.html",
        zip1="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_q_fastqc.zip",
        html2="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_p_fastqc.html",
        zip2="qc/2trimmed/{method}-{condition}-{replicate}-trimmed_p_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    params:
        prefix1=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads1))[0]),
        prefix2=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads2))[0])
    shell:
        """
        mkdir -p qc/2trimmed;
        fastqc -o qc/2trimmed -t {threads} {input}; mv qc/2trimmed/{params.prefix1}_fastqc.html {output.html1}; mv qc/2trimmed/{params.prefix1}_fastqc.zip {output.zip1}
        fastqc -o qc/2trimmed -t {threads} {input}; mv qc/2trimmed/{params.prefix2}_fastqc.html {output.html2}; mv qc/2trimmed/{params.prefix2}_fastqc.zip {output.zip2}
        """

ruleorder: fastqcraw_paired > fastqcraw_single
ruleorder: fastqctrimmed_paired > fastqctrimmed_single

def get_qc_files():
    qc_files = []
    for index, row in samples.iterrows():
        if pd.isna(row['fastqFile2']):
            qc_files.append("qc/1raw/{method}-{condition}-{replicate}-raw_fastqc.html".format(**row))
            qc_files.append("qc/2trimmed/{method}-{condition}-{replicate}-trimmed_fastqc.html".format(**row))
            qc_files.append("trimmed/{method}-{condition}-{replicate}.fastq".format(**row))
        else:
            qc_files.append("qc/1raw/{method}-{condition}-{replicate}-raw-q_fastqc.html".format(**row))
            qc_files.append("qc/1raw/{method}-{condition}-{replicate}-raw-p_fastqc.html".format(**row))
            qc_files.append("qc/2trimmed/{method}-{condition}-{replicate}-trimmed_q_fastqc.html".format(**row))
            qc_files.append("qc/2trimmed/{method}-{condition}-{replicate}-trimmed_p_fastqc.html".format(**row))
            qc_files.append("trimmed/{method}-{condition}-{replicate}_q.fastq".format(**row))
            qc_files.append("trimmed/{method}-{condition}-{replicate}_p.fastq".format(**row))

    return qc_files

rule multiqc:
    input:
        get_qc_files(),
        expand("qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/4unique/{method}-{condition}-{replicate}-map_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/all/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/trnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainallaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        expand("qc/rrnainall/{method}-{condition}-{replicate}.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    output:
        report("qc/multi/multiqc_report.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        dir="qc/multi"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "export LC_ALL=en_US.utf8; export LANG=en_US.utf8; multiqc -f -d --exclude picard --exclude gatk -z -o {params.dir} qc/1raw qc/2trimmed qc/3mapped qc/4unique qc/5removedrRNA qc/all qc/trnainall qc/rrnainallaligned qc/rrnainuniquelyaligned qc/rrnainall trimmed  2> {log}"
