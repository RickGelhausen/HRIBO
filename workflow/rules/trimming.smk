def sample_row(wildcards):
    """The single sample sheet row addressed by the wildcards."""
    match = samples[
        (samples["method"] == wildcards.method)
        & (samples["condition"] == wildcards.condition)
        & (samples["replicate"] == wildcards.replicate)
    ]
    if len(match) != 1:
        raise ValueError(
            f"Expected exactly one sample sheet row for "
            f"{wildcards.method}-{wildcards.condition}-{wildcards.replicate}, found {len(match)}."
        )
    return match.iloc[0]


def get_inputs_single(wildcards):
    row = sample_row(wildcards)
    if is_paired_end(row):
        raise ValueError(f"{library_name(wildcards)} is paired-end; the single-end rule does not apply.")
    return row["fastqFile"]


def get_inputs_paired(wildcards):
    row = sample_row(wildcards)
    if not is_paired_end(row):
        raise ValueError(f"{library_name(wildcards)} is single-end; the paired-end rule does not apply.")
    return [row["fastqFile"], row["fastqFile2"]]

rule link_single:
    input:
        fastq=get_inputs_single,
        stager=str(SCRIPTS / "stage_input.py")
    output:
        fastq="trimlink/{method}-{condition}-{replicate}.fastq.gz"
    threads: 1
    shell:
        "python3 {input.stager:q} link {input.fastq:q} {output.fastq:q}"

rule link_paired:
    input:
        fastq1=lambda wildcards: get_inputs_paired(wildcards)[0],
        fastq2=lambda wildcards: get_inputs_paired(wildcards)[1],
        stager=str(SCRIPTS / "stage_input.py")
    output:
        fastq1="trimlink/{method}-{condition}-{replicate}_q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_p.fastq.gz"
    threads: 1
    shell:
        """
        python3 {input.stager:q} link {input.fastq1:q} {output.fastq1:q}
        python3 {input.stager:q} link {input.fastq2:q} {output.fastq2:q}
        """

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
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "cutadapt -j {threads} {params.adapter3} {params.adapter5} {params.quality} {params.filtering} -o {output.fastq} {input.fastq}"

rule trim_paired:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_p.fastq.gz"
    output:
        fastq1=temp("trimmedpaired/{method}-{condition}-{replicate}_q.fastq"),
        fastq2=temp("trimmedpaired/{method}-{condition}-{replicate}_p.fastq")
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
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "cutadapt -j {threads} {params.adapter3q} {params.adapter5q} {params.adapter3p} {params.adapter5p} {params.quality} {params.filtering} -o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq2}"

rule merge_fastq:
    input:
        fastq1="trimmedpaired/{method}-{condition}-{replicate}_q.fastq",
        fastq2="trimmedpaired/{method}-{condition}-{replicate}_p.fastq"
    output:
        fastq="trimmed/{method}-{condition}-{replicate}.fastq"
    conda:
        "../envs/pear.yaml"
    threads: 20
    resources:
        mem_mb=20000,
        runtime=120
    log:
        "logs/{method}-{condition}-{replicate}_pear.log"
    shell:
        """
        mkdir -p pear
        pear -n 10 -f {input.fastq1} -r {input.fastq2} -o pear/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}
        mv pear/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}.assembled.fastq {output.fastq}
        """

ruleorder: trim_single > merge_fastq
