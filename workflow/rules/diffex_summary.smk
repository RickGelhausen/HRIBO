"""Cross-condition detection and differential-expression overview."""


rule conditionOverview:
    input:
        counts="readcounts/differential_expression_read_counts.csv",
        annotation="readcounts/independant_annotation.gff",
        xtail="xtail/xtail_all.csv",
        riborex="riborex/riborex_all.csv",
        deltate="deltate/deltate_all.csv",
        script=str(SCRIPTS / "generate_diffex_summary.py"),
        script_deps=[
            str(SCRIPTS / "lib" / "__init__.py"),
            str(SCRIPTS / "lib" / "diffex_browser.py"),
        ]
    output:
        html="diffex_summary/condition_overview.html",
        xlsx="diffex_summary/condition_overview.xlsx",
        condition_matrix="diffex_summary/condition_matrix.tsv",
        contrast_matrix="diffex_summary/contrast_matrix.tsv",
        browser_manifest="diffex_summary/browser_tracks.tsv",
        browser=directory("diffex_summary/browser")
    conda:
        "../envs/excel.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    params:
        contrasts=CONTRASTS,
        min_cpm=config["differentialExpressionSettings"].get("detectionMinCPM", 1),
        min_count=config["differentialExpressionSettings"].get("detectionMinCount", 10),
        min_replicates=config["differentialExpressionSettings"].get("detectionMinReplicates", 2),
        padj_cutoff=config["differentialExpressionSettings"]["padjCutoff"],
        log2fc_cutoff=config["differentialExpressionSettings"]["log2fcCutoff"]
    log:
        "logs/diffex_summary.log"
    shell:
        """
        python3 {input.script:q} \
            --counts {input.counts:q} \
            --annotation {input.annotation:q} \
            --xtail {input.xtail:q} \
            --riborex {input.riborex:q} \
            --deltate {input.deltate:q} \
            --contrasts {params.contrasts:q} \
            --output_dir diffex_summary \
            --min_cpm {params.min_cpm:q} \
            --min_count {params.min_count:q} \
            --min_replicates {params.min_replicates:q} \
            --padj_cutoff {params.padj_cutoff:q} \
            --log2fc_cutoff {params.log2fc_cutoff:q} \
            > {log:q} 2>&1
        """
