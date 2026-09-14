from pathlib import Path

rule readLengthStatistics:
    input:
        bamfiles=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples_metagene["method"], condition=samples_metagene["condition"], replicate=samples_metagene["replicate"]),
        bamIndex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples_metagene["method"], condition=samples_metagene["condition"], replicate=samples_metagene["replicate"]),
        script=str(SCRIPTS / "read_length_statistics.py"),
        script_deps=[
            str(SCRIPTS / "lib" / "__init__.py"),
            str(SCRIPTS / "lib" / "alignment.py"),
            str(SCRIPTS / "lib" / "io.py"),
            str(SCRIPTS / "lib" / "theme.py"),
        ]
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
        python3 {input.script:q} \
            -a {input.bamfiles:q} \
            -r {params.readlengths:q} \
            -o "metageneprofiling/" \
            > {log:q} 2>&1
        """

rule metageneProfiling:
    input:
        bam=rules.maplink.output,
        bamIndex=rules.bamindex.output,
        genome=rules.retrieveGenome.output,
        annotation=rules.checkAnnotation.output,
        script=str(SCRIPTS / "metagene_profiling.py"),
        script_deps=[
            str(SCRIPTS / "lib" / "__init__.py"),
            str(SCRIPTS / "lib" / "alignment.py"),
            str(SCRIPTS / "lib" / "annotation.py"),
            str(SCRIPTS / "lib" / "io.py"),
            str(SCRIPTS / "lib" / "metagene.py"),
            str(SCRIPTS / "lib" / "misc.py"),
            str(SCRIPTS / "lib" / "plotting.py"),
            str(SCRIPTS / "lib" / "psite.py"),
            str(SCRIPTS / "lib" / "theme.py"),
        ]
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
        colorArgs=(
            ["--color_list", *config["metageneSettings"]["colorList"]]
            if config["metageneSettings"]["colorList"]
            else []
        )
    log: "logs/{method}-{condition}-{replicate}_metageneprofiling.log"
    shell:
        """
        python3 {input.script:q} \
            -b {input.bam:q} \
            -g {input.genome:q} \
            -a {input.annotation:q} \
            -o {output.meta:q} \
            --read_lengths {params.readlengths:q} \
            --normalization_methods {params.normalizationMethods:q} \
            --mapping_methods {params.mappingMethods:q} \
            --positions_in_ORF {params.positionsInORF:q} \
            --positions_out_ORF {params.positionsOutORF:q} \
            --filtering_methods {params.filteringMethods:q} \
            --neighboring_genes_distance {params.neighboringGenesDistance:q} \
            --rpkm_threshold {params.rpkmThreshold:q} \
            --length_cutoff {params.lengthCutoff:q} \
            --output_formats {params.outputFormats:q} \
            --include_plotly_js {params.includePlotlyJS:q} \
            {params.colorArgs:q} \
            > {log:q} 2>&1
        """


rule tisAdvisor:
    input:
        bam=rules.maplink.output,
        bamIndex=rules.bamindex.output,
        genome=rules.retrieveGenome.output,
        annotation=rules.checkAnnotation.output,
        script=str(SCRIPTS / "tis_advisor.py"),
        script_deps=[
            str(SCRIPTS / "lib" / "__init__.py"),
            str(SCRIPTS / "lib" / "alignment.py"),
            str(SCRIPTS / "lib" / "annotation.py"),
            str(SCRIPTS / "lib" / "io.py"),
            str(SCRIPTS / "lib" / "metagene.py"),
            str(SCRIPTS / "lib" / "misc.py"),
            str(SCRIPTS / "lib" / "plotting.py"),
            str(SCRIPTS / "lib" / "psite.py"),
            str(SCRIPTS / "lib" / "theme.py"),
        ]
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
        mappingMethods=config["tisAdvisorSettings"]["mappingMethods"],
        positionsInORF=config["metageneSettings"]["positionsInORF"],
        positionsOutORF=config["metageneSettings"]["positionsOutsideORF"],
        filteringMethods=config["metageneSettings"]["filteringMethods"],
        neighboringGenesDistance=config["metageneSettings"]["neighboringGenesDistance"],
        rpkmThreshold=config["metageneSettings"]["rpkmThreshold"],
        lengthCutoff=config["metageneSettings"]["lengthCutoff"],
        includePlotlyJS=config["metageneSettings"]["includePlotlyJS"],
        deepriboASiteOffset=config["predictionSettings"]["deepriboASiteOffset"]
    log:
        "logs/{method}-{condition}-{replicate}_tis_advisor.log"
    shell:
        """
        python3 {input.script:q} \
            -b {input.bam:q} \
            -a {input.annotation:q} \
            -g {input.genome:q} \
            -o {params.outdir:q} \
            -r {params.readlengths:q} \
            --mapping_methods {params.mappingMethods:q} \
            --positions_in_ORF {params.positionsInORF:q} \
            --positions_out_ORF {params.positionsOutORF:q} \
            --filtering_methods {params.filteringMethods:q} \
            --neighboring_genes_distance {params.neighboringGenesDistance:q} \
            --rpkm_threshold {params.rpkmThreshold:q} \
            --length_cutoff {params.lengthCutoff:q} \
            --include_plotly_js {params.includePlotlyJS:q} \
            --deepribo_asite_offset {params.deepriboASiteOffset:q} \
            > {log:q} 2>&1
        """
