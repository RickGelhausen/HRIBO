from pathlib import Path

rule normalizedmetageneprofilingTIS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.checkAnnotation.output
    output:
        meta=directory("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "mkdir -p metageneprofiling/TIS/norm; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TIS/norm/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --normalization"


rule metageneprofilingTIS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.checkAnnotation.output
    output:
        directory("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "mkdir -p metageneprofiling/TIS/raw; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TIS/raw/{params.prefix} --in_fai_filepath genomes/genome.fa.fai"

rule normalizedmetageneprofilingTTS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.checkAnnotation.output
    output:
        meta=directory("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "mkdir -p metageneprofiling/TTS/norm; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TTS/norm/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --normalization --input_type TTS"


rule metageneprofilingTTS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.checkAnnotation.output
    output:
        directory("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "mkdir -p metageneprofiling/TTS/raw; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TTS/raw/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --input_type TTS"
