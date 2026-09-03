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
        fastq="trimmed/{method}-{condition}-{replicate}.fastq"
    params:
        adapter3=command_option_values("-a", ADAPTERS_S3),
        adapter5=command_option_values("-g", ADAPTERS_S5),
        quality=["-q", "20", "--trim-n"],
        filtering=["-m", "10"]
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "cutadapt -j {threads} {params.adapter3:q} {params.adapter5:q} {params.quality:q} {params.filtering:q} -o {output.fastq:q} {input.fastq:q}"

rule trim_paired:
    input:
        fastq1="trimlink/{method}-{condition}-{replicate}_q.fastq.gz",
        fastq2="trimlink/{method}-{condition}-{replicate}_p.fastq.gz"
    output:
        fastq1=temp("trimmedpaired/{method}-{condition}-{replicate}_q.fastq"),
        fastq2=temp("trimmedpaired/{method}-{condition}-{replicate}_p.fastq")
    params:
        adapter3q=command_option_values("-a", ADAPTERS_P3R1),
        adapter5q=command_option_values("-g", ADAPTERS_P5R1),
        adapter3p=command_option_values("-A", ADAPTERS_P3R2),
        adapter5p=command_option_values("-G", ADAPTERS_P5R2),
        quality=["-q", "20", "--trim-n"],
        filtering=["-m", "10"]
    conda:
        "../envs/cutadapt.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=120
    shell:
        "cutadapt -j {threads} {params.adapter3q:q} {params.adapter5q:q} {params.adapter3p:q} {params.adapter5p:q} {params.quality:q} {params.filtering:q} -o {output.fastq1:q} -p {output.fastq2:q} {input.fastq1:q} {input.fastq2:q}"

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
    params:
        outdir="pear",
        prefix=lambda wildcards: f"pear/{library_name(wildcards)}",
        assembled=lambda wildcards: f"pear/{library_name(wildcards)}.assembled.fastq"
    log:
        "logs/{method}-{condition}-{replicate}_pear.log"
    shell:
        """
        exec > {log:q} 2>&1
        mkdir -p {params.outdir:q}
        pear -n 10 -f {input.fastq1:q} -r {input.fastq2:q} -o {params.prefix:q}
        mv {params.assembled:q} {output.fastq:q}
        """

ruleorder: trim_single > merge_fastq
