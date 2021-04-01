from pathlib import Path

rule normalizedmetageneprofilingTIS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        readlengthstat="metageneprofiling/readlengthstats/{method}-{condition}-{replicate}.bam_read_length_distribution.json",
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
        "mkdir -p metageneprofiling/TIS/norm; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TIS/norm/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --normalization --in_readlengthstat_filepath {input.readlengthstat} --noise_reduction_analysis"

rule readlengthstat:
    input:
        bam=rules.maplink.output,
        bamIndex=rules.bamindex.output
    output:
        "metageneprofiling/readlengthstats/{method}-{condition}-{replicate}.bam_read_length_distribution.json"
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]) + "/")
    shell:
        "mkdir -p metageneprofiling/readlengthstats; HRIBO/scripts/readlengthstat.py --in_bam_filepath {input.bam} --out_folder_filepath {params.prefix}"


rule metageneprofilingTIS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        readlengthstat="metageneprofiling/readlengthstats/{method}-{condition}-{replicate}.bam_read_length_distribution.json",
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
        "mkdir -p metageneprofiling/TIS/raw; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TIS/raw/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --in_readlengthstat_filepath {input.readlengthstat} --noise_reduction_analysis"

rule normalizedmetageneprofilingTTS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        readlengthstat="metageneprofiling/readlengthstats/{method}-{condition}-{replicate}.bam_read_length_distribution.json",
        bamIndex=rules.bamindex.output,
        annotation=rules.checkAnnotation.output
    output:
        meta=directory("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem),
    shell:
        "mkdir -p metageneprofiling/TTS/norm; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TTS/norm/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --normalization --input_type TTS --in_readlengthstat_filepath {input.readlengthstat} --noise_reduction_analysis"


rule metageneprofilingTTS:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        readlengthstat="metageneprofiling/readlengthstats/{method}-{condition}-{replicate}.bam_read_length_distribution.json",
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
        "mkdir -p metageneprofiling/TTS/raw; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/TTS/raw/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --input_type TTS --in_readlengthstat_filepath {input.readlengthstat} --noise_reduction_analysis"

rule merged_offsets:
    input:
        expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
        expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
        expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
        expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
    output:
        "metageneprofiling/merged_offsets.json"
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0]) + "/")
    shell:
        "mkdir -p metageneprofiling/; HRIBO/scripts/merged_offsets.py --in_metagene_directorypath metageneprofiling --out_filepath metageneprofiling"

