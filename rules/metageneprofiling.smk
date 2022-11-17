from pathlib import Path

rule metageneProfiling:
    input:
        bam=rules.maplink.output,
        bamIndex=rules.bamindex.output,
        genome=rules.retrieveGenome.output,
        annotation=rules.checkAnnotation.output
    output:
        meta=directory("metageneprofiling/{method}-{condition}-{replicate}")
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    params:
        readlengths=config["metageneSettings"]["readLengths"],
        positionsInORF=config["metageneSettings"]["positionsInORF"],
        positionsOutORF=config["metageneSettings"]["positionsOutsideORF"],
        filteringMethod=" ".join(config["metageneSettings"]["filteringMethod"]),
        neighboringGenesDistance=config["metageneSettings"]["neighboringGenesDistance"],
        rpkmThreshold=config["metageneSettings"]["rpkmThreshold"],
        lengthCutoff=config["metageneSettings"]["lengthCutoff"],
        mappingMethods=" ".join(config["metageneSettings"]["mappingMethods"]),
        normalizationMethod=config["metageneSettings"]["normalizationMethod"],
        outputFormats=" ".join(config["metageneSettings"]["outputFormats"]),
        includePlotlyJS=config["metageneSettings"]["includePlotlyJS"],
        colorList= [] if len(config["metageneSettings"]["colorList"]) == 0 else " ".join(config["metageneSettings"]["colorList"])
    shell:
        """
        mkdir -p metageneprofiling;
        HRIBO/scripts/metagene_profiling.py -b {input.bam} -g {input.genome} -a {input.annotation} -o {output.meta} \
            --read_lengths {params.readlengths} \
            --normalization_methods {params.normalizationMethod} \
            --mapping_methods {params.mappingMethods} \
            --positions_in_ORF {params.positionsInORF} \
            --positions_out_ORF {params.positionsOutORF} \
            --filtering_method {params.filteringMethod} \
            --neighboring_genes_distance {params.neighboringGenesDistance} \
            --rpkm_threshold {params.rpkmThreshold} \
            --length_cutoff {params.lengthCutoff} \
            --color_list {params.colorList} \
            --output_formats {params.outputFormats} \
            --include_plotly_js {params.includePlotlyJS}
        """

