DELTATE_CONTAINER = (
    "docker://gelhausr/deltate@sha256:"
    "f7611403ff14417b495ce1806ec88416dd5f6b7cfeb0a8ad4fa6300889265492"
)


rule prepareDeltaTEScript:
    input:
        patcher=str(SCRIPTS / "patch_deltate.py")
    output:
        script="deltate/DTEG.R"
    container:
        DELTATE_CONTAINER
    threads: 1
    resources:
        mem_mb=256,
        runtime=1
    shell:
        """
        python3 {input.patcher:q} /usr/local/bin/DTEG.R {output.script:q}
        """


rule deltatePrepareInput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}",
        bam=expand("bam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        script=str(SCRIPTS / "prepare_deltate_input.py")
    output:
        ribo="deltate/{contrast}/ribo_counts.txt",
        rna="deltate/{contrast}/rna_counts.txt",
        samples="deltate/{contrast}/samples_info.txt",
        replicates="deltate/{contrast}/has_replicates.txt"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        contrast = lambda wildcards, input: input[1].split("/")[1],
        out_dir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        python3 {input.script:q} -c {params.contrast:q} -r {input.rawreads:q} --bam_files {input.bam:q} -o {params.out_dir:q}
        """

rule deltate:
    input:
        contrastfile="contrasts/{contrast}",
        ribo="deltate/{contrast}/ribo_counts.txt",
        rna="deltate/{contrast}/rna_counts.txt",
        samples="deltate/{contrast}/samples_info.txt",
        replicates="deltate/{contrast}/has_replicates.txt",
        runner=str(SCRIPTS / "run_deltate.sh"),
        engine=rules.prepareDeltaTEScript.output.script
    output:
        fcribo=ensure("deltate/{contrast}/fold_changes/deltaRibo.txt", non_empty=True),
        fcrna=ensure("deltate/{contrast}/fold_changes/deltaRNA.txt", non_empty=True),
        fcte=ensure("deltate/{contrast}/fold_changes/deltaTE.txt", non_empty=True),
        fig=ensure("deltate/{contrast}_figures.pdf", non_empty=True)
    container:
        DELTATE_CONTAINER
    threads: 1
    params:
        result_dir=lambda wildcards: f"deltate/{wildcards.contrast}",
        result_fig=lambda wildcards: f"deltate/{wildcards.contrast}/Result_figures.pdf"
    log:
        "logs/deltate/{contrast}.log"
    shell:
        """
        bash {input.runner:q} \
            {input.ribo:q} \
            {input.rna:q} \
            {input.samples:q} \
            {params.result_dir:q} \
            {output.fcribo:q} \
            {output.fcrna:q} \
            {output.fcte:q} \
            {params.result_fig:q} \
            {output.fig:q} \
            {input.engine:q} \
            > {log:q} 2>&1
        """

rule deltatexlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        deltate_ribo="deltate/{contrast}/fold_changes/deltaRibo.txt",
        deltate_rna="deltate/{contrast}/fold_changes/deltaRNA.txt",
        deltate_te="deltate/{contrast}/fold_changes/deltaTE.txt",
        script=str(SCRIPTS / "generate_excel_deltate.py"),
        script_deps=[
            str(SCRIPTS / "excel_utils.py"),
            str(SCRIPTS / "gff_utils.py"),
        ]
    output:
        xlsx_sorted="deltate/{contrast}_sorted.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    params:
        padj_cutoff=config["differentialExpressionSettings"]["padjCutoff"],
        log2fc_cutoff=config["differentialExpressionSettings"]["log2fcCutoff"]
    shell:
        """
        python3 {input.script:q} -a {input.annotation:q} -g {input.genome:q} -i {input.deltate_ribo:q} -r {input.deltate_rna:q} -t {input.deltate_te:q} -o {output.xlsx_sorted:q} --padj_cutoff {params.padj_cutoff:q} --log2fc_cutoff {params.log2fc_cutoff:q}
        """

rule pooldeltate:
    input:
        deltate=expand("deltate/{contr}_sorted.xlsx", contr=CONTRASTS),
        script=str(SCRIPTS / "merge_differential_expression.py")
    output:
        "deltate/deltate_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 {input.script:q} {input.deltate:q} -o {output:q} -t deltate
        """
