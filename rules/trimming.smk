def get_inputs(wildcards):
    row = samples[(samples['method'] == wildcards.method) & (samples['condition'] == wildcards.condition) & (samples['replicate'] == wildcards.replicate)]
    if pd.isna(row['fastqFile2'].iloc[0]):
        return row['fastqFile'].iloc[0]
    else:
        return [row['fastqFile'].iloc[0], row['fastqFile2'].iloc[0]]

rule link_single:
    input:
        get_inputs
    output:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
    params:
        prefix=lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.basename(input[0]))[0]),
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output.fastq))
    threads: 1
    shell:
        "mkdir -p trimlink; ln -s {params.inlink} {params.outlink};"

rule link_paired:
    input:
        get_inputs
    output:
        fastq1="trimlink/{method}-{condition}-{replicate}_q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_p.fastq.gz"
    params:
        prefix1=lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.basename(input[0]))[0]),
        prefix2=lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.basename(input[1]))[0]),
        inlink1=lambda wildcards, input:(os.getcwd() + "/" + str(input[0])),
        inlink2=lambda wildcards, input:(os.getcwd() + "/" + str(input[1])),
        outlink1=lambda wildcards, output:(os.getcwd() + "/" + str(output.fastq1)),
        outlink2=lambda wildcards, output:(os.getcwd() + "/" + str(output.fastq2))
    threads: 1
    shell:
        "mkdir -p trimlink; ln -s {params.inlink1} {params.outlink1}; ln -s {params.inlink2} {params.outlink2};"

ruleorder: link_paired > link_single


rule trim_single:
    input:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
    output:
        fastq=temp("trimmed/{method}-{condition}-{replicate}.fastq")
    params:
        adapter3=lambda wildcards, output: ("" if not ADAPTERS_S3 else (" ".join([" -a %s" % adapter for adapter in ADAPTERS_S3.split(",")]))),
        adapter5=lambda wildcards, output: ("" if not ADAPTERS_S5 else (" ".join([" -g %s" % adapter for adapter in ADAPTERS_S5.split(",")]))),
        quality=" -q 20 --trim-n ",
        filtering=" -m 10 "
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    shell:
        "mkdir -p trimmed; cutadapt -j {threads} {params.adapter3} {params.adapter5} {params.quality} {params.filtering} -o {output.fastq} {input.fastq}"

rule trim_paired:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_p.fastq.gz"
    output:
        fastq1=temp("trimmed/{method}-{condition}-{replicate}_q.fastq"),
        fastq2=temp("trimmed/{method}-{condition}-{replicate}_p.fastq")
    params:
        adapter3q=lambda wildcards, output: ("" if not ADAPTERS_P3R1 else (" ".join([" -a %s" % adapter for adapter in ADAPTERS_P3R1.split(",")]))),
        adapter5q=lambda wildcards, output: ("" if not ADAPTERS_P5R1 else (" ".join([" -g %s" % adapter for adapter in ADAPTERS_P5R1.split(",")]))),
        adapter3p=lambda wildcards, output: ("" if not ADAPTERS_P3R2 else (" ".join([" -A %s" % adapter for adapter in ADAPTERS_P3R2.split(",")]))),
        adapter5p=lambda wildcards, output: ("" if not ADAPTERS_P5R2 else (" ".join([" -G %s" % adapter for adapter in ADAPTERS_P5R2.split(",")]))),
        quality=" -q 20 --trim-n ",
        filtering=" -m 10 "
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    shell:
        "mkdir -p trimmed; cutadapt -j {threads} {params.adapter3q} {params.adapter5q} {params.adapter3p} {params.adapter5p} {params.quality} {params.filtering} -o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq2}"

ruleorder: trim_paired > trim_single