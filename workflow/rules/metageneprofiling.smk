from pathlib import Path

rule readLengthStatistics:
    input:
        bamfiles=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples_metagene["method"], condition=samples_metagene["condition"], replicate=samples_metagene["replicate"]),
        bamIndex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples_metagene["method"], condition=samples_metagene["condition"], replicate=samples_metagene["replicate"])
    output:
        plot="metageneprofiling/read_length_fractions.html",
        table="metageneprofiling/read_length_fractions.xlsx",
        raw="metageneprofiling/read_length_counts.xlsx"
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    resources:
        mem_mb=20000,
        runtime=60
    log: "logs/read_length_statistics.log"
    params:
        readlengths=config["readstatSettings"]["readLengths"]
    shell:
        """
        {SCRIPTS}/read_length_statistics.py -a {input.bamfiles} -r {params.readlengths} -o metageneprofiling/ > {log}
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
    resources:
        mem_mb=20000,
        runtime=120
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
    log: "logs/{method}-{condition}-{replicate}_metageneprofiling.log"
    shell:
        """
        if [ {params.colorList} == nocolor ]; then
            colorList="";
        else
            colorList="--color_list {params.colorList}";
        fi;
        {SCRIPTS}/metagene_profiling.py -b {input.bam} -g {input.genome} -a {input.annotation} -o {output.meta} \
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


rule tisAdvisor:
    input:
        bam=rules.maplink.output,
        bamIndex=rules.bamindex.output,
        genome=rules.retrieveGenome.output,
        annotation=rules.checkAnnotation.output
    output:
        report_html=report(
            "tis_advice/{method}-{condition}-{replicate}/tis_recommendation.html",
            caption="../report/tisadvice.rst",
            category="TIS caller advice",
            labels={"library": "{method}-{condition}-{replicate}"}
        ),
        recommendation="tis_advice/{method}-{condition}-{replicate}/tis_recommendation.json",
        evidence="tis_advice/{method}-{condition}-{replicate}/read_length_evidence.tsv"
    conda:
        "../envs/metageneprofiling.yaml"
    threads: 1
    resources:
        mem_mb=20000,
        runtime=120
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.report_html),
        readlengths=config["tisAdvisorSettings"]["readLengths"],
        mappingMethod=config["tisAdvisorSettings"]["mappingMethod"],
        positionsInORF=config["metageneSettings"]["positionsInORF"],
        positionsOutORF=config["metageneSettings"]["positionsOutsideORF"],
        filteringMethods=config["metageneSettings"]["filteringMethods"],
        neighboringGenesDistance=config["metageneSettings"]["neighboringGenesDistance"],
        rpkmThreshold=config["metageneSettings"]["rpkmThreshold"],
        includePlotlyJS=config["metageneSettings"]["includePlotlyJS"]
    log:
        "logs/{method}-{condition}-{replicate}_tis_advisor.log"
    shell:
        """
        {SCRIPTS}/tis_advisor.py \
            -b {input.bam} \
            -a {input.annotation} \
            -g {input.genome} \
            -o {params.outdir} \
            -r {params.readlengths} \
            --mapping_method {params.mappingMethod} \
            --positions_in_ORF {params.positionsInORF} \
            --positions_out_ORF {params.positionsOutORF} \
            --filtering_methods {params.filteringMethods} \
            --neighboring_genes_distance {params.neighboringGenesDistance} \
            --rpkm_threshold {params.rpkmThreshold} \
            --include_plotly_js {params.includePlotlyJS} > {log} 2>&1
        """
