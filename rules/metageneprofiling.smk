from pathlib import Path

rule readLengthStatistics:
    input:
        bamfiles=expand(maplink/{method}-{condition}-{replicate}.bam, method=config["methods"], condition=config["conditions"], replicate=config["replicates"]),
        bamIndex=expand(bamindex/{method}-{condition}-{replicate}.bam.bai, method=config["methods"], condition=config["conditions"], replicate=config["replicates"])
    output:
        plot="metageneprofiling/read_length_fractions.html",
        table="metageneprofiling/read_length_fractions.tsv",
        raw="metageneprofiling/read_length_counts.tsv"
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    log: "logs/read_length_statistics.log"
    params:
        readlengths=config["readstatSettings"]["readLengths"]
    shell:
        """
        mkdir -p metageneprofiling;
        HRIBO/scripts/read_length_statistics.py -b {input.bamfiles} -r {params.readlengths} -o metageneprofiling/ > {log}
        """

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
        filteringMethods=config["metageneSettings"]["filteringMethods"],
        neighboringGenesDistance=config["metageneSettings"]["neighboringGenesDistance"],
        rpkmThreshold=config["metageneSettings"]["rpkmThreshold"],
        lengthCutoff=config["metageneSettings"]["lengthCutoff"],
        mappingMethods=config["metageneSettings"]["mappingMethods"],
        normalizationMethods=config["metageneSettings"]["normalizationMethods"],
        outputFormats=config["metageneSettings"]["outputFormats"],
        includePlotlyJS=config["metageneSettings"]["includePlotlyJS"],
        colorList= "nocolor" if len(config["metageneSettings"]["colorList"]) == 0 else config["metageneSettings"]["colorList"]
    log: "logs/{method}-{condition}-{replicate}.log"
    shell:
        """
        mkdir -p metageneprofiling;
        if [ {params.colorList} == nocolor ]; then
            colorList="";
        else
            colorList="--color_list {params.colorList}";
        fi;
        HRIBO/scripts/metagene_profiling.py -b {input.bam} -g {input.genome} -a {input.annotation} -o {output.meta} \
            --read_lengths {params.readlengths} \
            --normalization_methods {params.normalizationMethods} \
            --mapping_methods {params.mappingMethods} \
            --positions_in_ORF {params.positionsInORF} \
            --positions_out_ORF {params.positionsOutORF} \
            --filtering_method {params.filteringMethods} \
            --neighboring_genes_distance {params.neighboringGenesDistance} \
            --rpkm_threshold {params.rpkmThreshold} \
            --length_cutoff {params.lengthCutoff} \
            --output_formats {params.outputFormats} \
            --include_plotly_js {params.includePlotlyJS} \
            ${{colorList}}; > {log}
        """

