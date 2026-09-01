rule deltatePrepareInput:
    input:
        rawreads="readcounts/differential_expression_read_counts.csv",
        contrastfile="contrasts/{contrast}",
        bam=expand("bam/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
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
        {SCRIPTS}/prepare_deltate_input.py -c {params.contrast} -r {input.rawreads} -b bam/ -o {params.out_dir}
        """

rule deltate:
    input:
        contrastfile="contrasts/{contrast}",
        ribo="deltate/{contrast}/ribo_counts.txt",
        rna="deltate/{contrast}/rna_counts.txt",
        samples="deltate/{contrast}/samples_info.txt",
        replicates="deltate/{contrast}/has_replicates.txt",
        runner=str(SCRIPTS / "run_deltate.sh")
    output:
        fcribo=ensure("deltate/{contrast}/fold_changes/deltaRibo.txt", non_empty=True),
        fcrna=ensure("deltate/{contrast}/fold_changes/deltaRNA.txt", non_empty=True),
        fcte=ensure("deltate/{contrast}/fold_changes/deltaTE.txt", non_empty=True),
        fig=ensure("deltate/{contrast}_figures.pdf", non_empty=True)
    container:
        "docker://gelhausr/deltate:latest"
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
            > {log:q} 2>&1
        """

rule deltatexlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        deltate_ribo="deltate/{contrast}/fold_changes/deltaRibo.txt",
        deltate_rna="deltate/{contrast}/fold_changes/deltaRNA.txt",
        deltate_te="deltate/{contrast}/fold_changes/deltaTE.txt"
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
        python3 {SCRIPTS}/generate_excel_deltate.py -a {input.annotation} -g {input.genome} -i {input.deltate_ribo} -r {input.deltate_rna} -t {input.deltate_te} -o {output.xlsx_sorted} --padj_cutoff {params.padj_cutoff} --log2fc_cutoff {params.log2fc_cutoff}
        """

rule pooldeltate:
    input:
        deltate=expand("deltate/{contr}_sorted.xlsx", contr=CONTRASTS)
    output:
        "deltate/deltate_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 {SCRIPTS}/merge_differential_expression.py {input.deltate} -o {output} -t deltate
        """
