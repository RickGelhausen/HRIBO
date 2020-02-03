from pathlib import Path

rule normalizedmetageneprofiling:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.retrieveAnnotation.output
    output: 
        meta=directory("metageneprofiling/norm/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "mkdir -p metageneprofiling/norm; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/norm/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --normalization"


rule metageneprofiling:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.retrieveAnnotation.output
    output:
        directory("metageneprofiling/raw/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "mkdir -p metageneprofiling/raw; HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metageneprofiling/raw/{params.prefix} --in_fai_filepath genomes/genome.fa.fai"

