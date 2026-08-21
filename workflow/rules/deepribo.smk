
def read_parameters(filename, idx):
    """One of the two cutoffs written by parameter_estimation.R.

    Raises rather than returning a placeholder: passing something like "failed"
    on to DeepRibo as -r/-c produces a far more confusing failure than stopping
    here does.
    """
    with open(filename) as handle:
        fields = handle.readline().strip().split(",")
    if len(fields) < 2 or not all(field.strip() for field in fields[:2]):
        raise ValueError(
            f"{filename} does not contain the expected 'min_RPKM,min_coverage' pair "
            f"(found {fields!r}). The S-curve estimation step probably failed."
        )
    return fields[idx].strip()

rule deepriboGetModel:
    output:
        "deepribo/DeepRibo_model_v1.pt"
    params:
        url="https://github.com/Biobix/DeepRibo/raw/master/models/DeepRibo_model_v1.pt"
    conda:
        "../envs/download.yaml"
    threads: 1
    retries: 3
    resources:
        mem_mb=1000,
        runtime=30
    log:
        "logs/deepriboGetModel.log"
    shell:
        "curl -sSL --fail --retry 3 --retry-delay 5 {params.url} -o {output} 2> {log}"


rule asiteOccupancy:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bai="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        asitefwd="coverage_deepribo/{condition}-{replicate}_asite_fwd.bedgraph",
        asiterev="coverage_deepribo/{condition}-{replicate}_asite_rev.bedgraph"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    params:
        prefix=lambda wildcards, output: output.asitefwd[: -len("_asite_fwd.bedgraph")]
    log:
        "logs/{condition}-{replicate}_asite_occupancy.log"
    shell:
        "{SCRIPTS}/coverage_deepribo.py --alignment_file {input.bam} --output_file_prefix {params.prefix} 2> {log}"

rule coverage:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bai="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        covfwd="coverage_deepribo/{condition}-{replicate}_cov_fwd.bedgraph",
        covrev="coverage_deepribo/{condition}-{replicate}_cov_rev.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    log:
        "logs/{condition}-{replicate}_deepribo_coverage.log"
    shell:
        """
        bedtools genomecov -bg -ibam {input.bam} -strand + > {output.covfwd} 2> {log}
        bedtools genomecov -bg -ibam {input.bam} -strand - > {output.covrev} 2>> {log}
        """

rule parseDeepRibo:
    input:
        covS= "coverage_deepribo/{condition}-{replicate}_cov_fwd.bedgraph",
        covAS= "coverage_deepribo/{condition}-{replicate}_cov_rev.bedgraph",
        asiteS= "coverage_deepribo/{condition}-{replicate}_asite_fwd.bedgraph",
        asiteAS= "coverage_deepribo/{condition}-{replicate}_asite_rev.bedgraph",
        genome= rules.retrieveGenome.output,
        annotation= rules.checkAnnotation.output
    output:
        "deepribo/{condition}-{replicate}/data_list.csv"
    container:
        "docker://gelhausr/deepribo:latest"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/{condition}-{replicate}_parse_deepribo.log"
    shell:
        """
        mkdir -p deepribo/{wildcards.condition}-{wildcards.replicate}/0/;
        mkdir -p deepribo/{wildcards.condition}-{wildcards.replicate}/1/;
        DataParser.py {input.covS} {input.covAS} {input.asiteS} {input.asiteAS} {input.genome} deepribo/{wildcards.condition}-{wildcards.replicate} -g {input.annotation} > {log} 2>&1
        """

rule parameterEstimation:
    input:
        "deepribo/{condition}-{replicate}/data_list.csv"
    output:
        "deepribo/{condition}-{replicate}/parameters.txt"
    container:
        "docker://gelhausr/deepribo:latest"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    params:
        # A per-library prefix: the previous constant "figure" meant every
        # library wrote its S-curve diagnostic to the same path.
        dest=lambda wildcards, output: os.path.join(os.path.dirname(output[0]), "s_curve")
    log:
        "logs/{condition}-{replicate}_parameter_estimation.log"
    shell:
        "Rscript {SCRIPTS}/parameter_estimation.R -f {input} -o {output} -d {params.dest} > {log} 2>&1"

rule predictDeepRibo:
    input:
        model= "deepribo/DeepRibo_model_v1.pt",
        data= "deepribo/{condition}-{replicate}/data_list.csv",
        parameter= "deepribo/{condition}-{replicate}/parameters.txt"
    output:
        "deepribo/{condition}-{replicate}/predictions.csv"
    container:
        "docker://gelhausr/deepribo:latest"
    threads: 10
    resources:
        mem_mb=20000,
        runtime=240
    params:
        rpkm=lambda wildcards, input: read_parameters(input.parameter, 0),
        cov=lambda wildcards, input: read_parameters(input.parameter, 1)
    log:
        "logs/{condition}-{replicate}_predict_deepribo.log"
    shell:
        """
        DeepRibo.py predict deepribo/ --pred_data {wildcards.condition}-{wildcards.replicate}/ -r {params.rpkm} -c {params.cov} --model {input.model} --dest {output} --num_workers {threads} > {log} 2>&1
        """

rule deepriboGFF:
    input:
        "deepribo/{condition}-{replicate}/predictions.csv"
    output:
        r"deepribo/{condition}-{replicate,\d+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/create_deepribo_gff.py -c {wildcards.condition} -r {wildcards.replicate} -i {input} -o {output}"

rule concatDeepRibo:
    input:
        lambda wildcards: expand("deepribo/{{condition}}-{replicate}.deepribo.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/concatenate_gff.py {input} -o {output}"

rule allDeepRibo:
    input:
        merged_gff=expand("tracks/{condition}.deepribo.gff", zip, condition=set(samples["condition"]))
    output:
        "tracks/deepribo_all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/concatenate_gff.py {input.merged_gff} -o {output}"

rule filterDeepRibo:
    input:
        ingff="tracks/deepribo_all.gff",
        annotation=rules.checkAnnotation.output
    output:
        merged="tracks/deepribo_merged.gff",
        plus="tracks/deepribo_merged_plus.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/merge_duplicates_deepribo.py -i {input.ingff} -o {output.merged} -a {input.annotation}"


rule createExcelSummaryDeepRibo:
    input:
        total="readcounts/bam_mapped_reads.txt",
        reads="readcounts/deepribo_annotation.gff",
        genome="genomes/genome.fa"
    output:
        "auxiliary/predictions_deepribo.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/generate_excel_deepribo.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule newAnnotationDeepRibo:
    input:
        reparation_orfs="tracks/reparation_annotated.gff",
        deepribo_orfs="tracks/deepribo_merged.gff",
        currentAnnotation=rules.checkAnnotation.output
    output:
        "tracks/totalAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        {SCRIPTS}/concatenate_gff.py {input.deepribo_orfs} {input.reparation_orfs} {input.currentAnnotation} -o {output}
        """
