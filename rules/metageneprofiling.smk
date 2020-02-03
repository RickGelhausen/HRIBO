from pathlib import Path

rule metageneprofiling:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.retrieveAnnotation.output
    output: 
        meta=directory("metageneprofiling/raw/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem)
    shell:
        "HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metagene/norm/{params.prefix} --in_fai_filepath genomes/genome.fa.fai --normalization"


rule normalizedmetageneprofiling:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        annotation=rules.retrieveAnnotation.output
    output:
        directory("metageneprofiling/norm/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig')
    shell:
        "HRIBO/scripts/metageneprofiling.py --in_bam_filepath {input.bam} --in_gff_filepath {input.annotation} --out_plot_filepath metagene/norm/raw/{params.prefix} --in_fai_filepath genomes/genome.fa.fai"

