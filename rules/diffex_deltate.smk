def read_has_replicates(filename):
    try:
        line = False
        with open(filename, "r") as f:
            line = bool(f.readline().strip())
        return line
    except FileNotFoundError:
        return "failed"

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
        "../envs/pytools.yaml"
    threads: 1
    shell:
        """
        mkdir -p deltate;
        HRIBO/scripts/prepare_deltate_input.py -c {contrast} -r {input.rawreads} -b bam/ -o deltate/{contrast}
        """

rule deltate:
    input:
        contrastfile="contrasts/{contrast}",
        ribo="deltate/{contrast}/ribo_counts.txt",
        rna="deltate/{contrast}/rna_counts.txt",
        samples="deltate/{contrast}/samples_info.txt",
        replicates="deltate/{contrast}/has_replicates.txt"
    output:
        fcribo="deltate/{contrast}/fold_changes/deltaRibo.txt",
        fcrna="deltate/{contrast}/fold_changes/deltaRNA.txt",
        fcte="deltate/{contrast}/fold_changes/deltaTE.txt",
        fig="deltate/{contrast}/Result_figures.pdf"
    singularity:
        "docker://gelhausr/deltate:latest"
    threads: 1
    params:
        has_replicates=lambda wildcards, input: read_has_replicates(input[2])
    shell:
        """
        mkdir -p deltate;
        if [ params.has_replicates ]
        then
            DTEG.R -c {input.ribo} {input.rna} {input.samples} 1 deltate/{contrast}/
        else
            touch {output.fcribo}
            touch {output.fcrna}
            touch {output.fcte}
            touch {output.fig}
        fi
        """

rule deltatexlsx:
    input:
        annotation=rules.checkAnnotation.output,
        genome=rules.retrieveGenome.output,
        deltate_ribo="deltate/{contrast}/fold_changes/deltaRibo.txt",
        deltate_rna="deltate/{contrast}/fold_changes/deltaRNA.txt,
        deltate_te="deltate/{contrast}/fold_changes/deltaTE.txt
    output:
        xlsx_sorted="deltate/{contrast}_sorted.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/generate_excel_deltate.py -a {input.annotation} -g {input.genome} -i {input.deltate_ribo} -r {input.deltate_rna} -t {input.deltate_te} -o {output.xlsx_sorted}
        """


cur_contrast=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(samples["condition"].unique(),2)))] for item in sublist]

rule pooldeltate:
    input:
        deltate=expand("deltate/{contr}_sorted.xlsx", contr=cur_contrast)
    output:
        "deltate/deltate_all.csv"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        """
        python3 HRIBO/scripts/merge_differential_expression.py {input.deltate} -o {output} -t deltate
        """