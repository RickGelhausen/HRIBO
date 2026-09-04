
import math
import os


DEEPRIBO_CONTAINER = (
    "docker://gelhausr/deepribo@sha256:"
    "4c4aef7d7d8eea790c6fcb360f69eda864888ea5effe24492ec88b8ef807a780"
)
DEEPRIBO_MODEL_URL = (
    "https://raw.githubusercontent.com/Biobix/DeepRibo/"
    "322dbd85e60c45466bb48ff428f01e1e7e7f8a3a/models/DeepRibo_model_v1.pt"
)
DEEPRIBO_MODEL_SHA256 = (
    "e3742ed7666f07eb72a28c35c70f2c62ac38f1d201e082b69b4f8cf35fbe4aa3"
)
DEEPRIBO_MODEL_SIZE = 8292882


def read_parameters(filename, idx, receipt=None):
    """Validate and return one cutoff written by parameter_estimation.R."""
    try:
        with open(filename) as handle:
            lines = handle.read().strip().splitlines()
    except FileNotFoundError as error:
        if receipt is not None and not os.path.exists(receipt):
            # Params functions are evaluated while constructing dry-run DAGs,
            # before the receipt-producing job has run. Point the exception at
            # that declared input so Snakemake represents the value as <TBD>.
            raise FileNotFoundError(
                2, "cutoff publication is not ready", receipt
            ) from error
        raise ValueError(
            f"{filename} is missing despite its completed cutoff receipt"
        ) from error
    fields = lines[0].split(",") if len(lines) == 1 else []
    if len(fields) != 2 or not all(field.strip() for field in fields):
        raise ValueError(
            f"{filename} must contain exactly one 'min_RPKM,min_coverage' pair "
            f"(found {fields!r}). The S-curve estimation step probably failed."
        )
    try:
        values = [float(field) for field in fields]
    except ValueError as exc:
        raise ValueError(
            f"{filename} contains a non-numeric DeepRibo cutoff: {fields!r}. "
            "The S-curve estimation step probably failed."
        ) from exc
    if not all(math.isfinite(value) for value in values):
        raise ValueError(
            f"{filename} contains a non-finite DeepRibo cutoff: {fields!r}. "
            "The S-curve estimation step probably failed."
        )
    if values[0] <= 0:
        raise ValueError(
            f"{filename} has a non-positive min_RPKM cutoff: {fields[0]!r}"
        )
    if not 0 <= values[1] <= 0.60:
        raise ValueError(
            f"{filename} has a min_coverage cutoff outside [0, 0.60]: "
            f"{fields[1]!r}"
        )
    return fields[idx].strip()

rule deepriboGetModel:
    input:
        fetcher=str(SCRIPTS / "fetch_verified.py")
    output:
        "deepribo/DeepRibo_model_v1.pt"
    params:
        url=DEEPRIBO_MODEL_URL,
        sha256=DEEPRIBO_MODEL_SHA256,
        size=DEEPRIBO_MODEL_SIZE
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
        """
        python3 {input.fetcher:q} \
            --url {params.url:q} \
            --sha256 {params.sha256:q} \
            --size {params.size} \
            --output {output:q} > {log:q} 2>&1
        """


rule prepareDeepRiboSCurveScript:
    input:
        patcher=workflow.source_path("../scripts/patch_deepribo_scurve.py")
    output:
        script="deepribo/s_curve_cutoff_estimation.R"
    container:
        DEEPRIBO_CONTAINER
    threads: 1
    resources:
        mem_mb=256,
        runtime=1
    shell:
        """
        python3 {input.patcher:q} \
            /usr/local/bin/s_curve_cutoff_estimation.R \
            {output.script:q}
        """


rule asiteOccupancy:
    input:
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bai="maplink/RIBO-{condition}-{replicate}.bam.bai",
        script=str(SCRIPTS / "coverage_deepribo.py")
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
        prefix=lambda wildcards, output: output.asitefwd[: -len("_asite_fwd.bedgraph")],
        offset=config["predictionSettings"]["deepriboASiteOffset"]
    log:
        "logs/{condition}-{replicate}_asite_occupancy.log"
    shell:
        """
        python3 {input.script:q} \
            --alignment_file {input.bam:q} \
            --output_file_prefix {params.prefix:q} \
            --offset {params.offset:q} 2> {log:q}
        """

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
        bedtools genomecov -bg -ibam {input.bam:q} -strand + > {output.covfwd:q} 2> {log:q}
        bedtools genomecov -bg -ibam {input.bam:q} -strand - > {output.covrev:q} 2>> {log:q}
        """

rule parseDeepRibo:
    input:
        covS= "coverage_deepribo/{condition}-{replicate}_cov_fwd.bedgraph",
        covAS= "coverage_deepribo/{condition}-{replicate}_cov_rev.bedgraph",
        asiteS= "coverage_deepribo/{condition}-{replicate}_asite_fwd.bedgraph",
        asiteAS= "coverage_deepribo/{condition}-{replicate}_asite_rev.bedgraph",
        genome= rules.retrieveGenome.output,
        annotation= rules.checkAnnotation.output,
        parser=workflow.source_path("../scripts/deepribo_data_parser.py")
    output:
        # DataParser writes the CSV plus paired *_seq.pt and *_reads.pt tensors
        # below 0/ and 1/.  They are one indivisible predictor input, so expose
        # the complete directory rather than leaving the tensors as untracked
        # side effects.  Keeping it outside the parameter/prediction directory
        # also lets Snakemake replace a stale parse without deleting later
        # results before their dependent jobs are rerun.  The tree is temporary:
        # Snakemake retains it through both consumers, removes it afterwards,
        # and recreates the complete tree if either downstream result is needed.
        parsed=temp(directory("deepribo/parsed/{condition}-{replicate}"))
    container:
        DEEPRIBO_CONTAINER
    threads: 1
    resources:
        mem_mb=16000,
        runtime=120
    params:
        class_zero=lambda wildcards, output: os.path.join(output.parsed, "0"),
        class_one=lambda wildcards, output: os.path.join(output.parsed, "1")
    log:
        "logs/{condition}-{replicate}_parse_deepribo.log"
    shell:
        """
        mkdir -p {params.class_zero:q} {params.class_one:q}
        python3 {input.parser:q} {input.covS:q} {input.covAS:q} {input.asiteS:q} {input.asiteAS:q} {input.genome:q} {output.parsed:q} -g {input.annotation:q} > {log:q} 2>&1
        """

rule parameterEstimation:
    input:
        launcher=workflow.source_path("../scripts/run_parameter_estimation.py"),
        parsed=rules.parseDeepRibo.output.parsed,
        wrapper_deps=[workflow.source_path("../scripts/parameter_estimation.R")],
        engine=rules.prepareDeepRiboSCurveScript.output.script
    output:
        # Snakemake removes ordinary outputs before a job and its update()
        # backup can itself survive a killed scheduler. Declare only a receipt;
        # the checked runner owns the sibling pair and publishes the receipt
        # only after both artifacts are durably committed.
        receipt=ensure(
            "deepribo/cutoffs/{condition}-{replicate}/.complete",
            non_empty=True
        )
    container:
        DEEPRIBO_CONTAINER
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    params:
        data=lambda wildcards, input: os.path.join(input.parsed, "data_list.csv"),
        parameters=lambda wildcards, output: os.path.join(
            os.path.dirname(output.receipt), "parameters.txt"
        ),
        plot=lambda wildcards, output: os.path.join(
            os.path.dirname(output.receipt), "s_curve.png"
        )
    log:
        "logs/{condition}-{replicate}_parameter_estimation.log"
    shell:
        """
        python3 {input.launcher:q} {input.wrapper_deps:q} \
            --file {params.data:q} \
            --out {params.parameters:q} \
            --plot {params.plot:q} \
            --receipt {output.receipt:q} \
            --engine {input.engine:q} \
            > {log:q} 2>&1
        """

rule predictDeepRibo:
    input:
        model= "deepribo/DeepRibo_model_v1.pt",
        parsed=rules.parseDeepRibo.output.parsed,
        cutoff_receipt=rules.parameterEstimation.output.receipt
    output:
        "deepribo/{condition}-{replicate}/predictions.csv"
    container:
        DEEPRIBO_CONTAINER
    threads: 10
    resources:
        mem_mb=20000,
        runtime=240
    params:
        rpkm=lambda wildcards, input: read_parameters(
            os.path.join(os.path.dirname(input.cutoff_receipt), "parameters.txt"),
            0,
            input.cutoff_receipt
        ),
        cov=lambda wildcards, input: read_parameters(
            os.path.join(os.path.dirname(input.cutoff_receipt), "parameters.txt"),
            1,
            input.cutoff_receipt
        ),
        # The pinned DeepRibo loader inserts separators between data_path and
        # pred_data itself, so pass the two path components without trailing `/`.
        prediction_data=lambda wildcards, input: os.path.basename(input.parsed),
        root=lambda wildcards, input: os.path.dirname(input.parsed)
    log:
        "logs/{condition}-{replicate}_predict_deepribo.log"
    shell:
        """
        DeepRibo.py predict {params.root:q} --pred_data {params.prediction_data:q} -r {params.rpkm:q} -c {params.cov:q} --model {input.model:q} --dest {output:q} --num_workers {threads} > {log:q} 2>&1
        """

rule deepriboGFF:
    input:
        predictions="deepribo/{condition}-{replicate}/predictions.csv",
        genome=rules.retrieveGenome.output,
        script=str(SCRIPTS / "create_deepribo_gff.py"),
        script_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        r"deepribo/{condition}-{replicate,\d+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -c {wildcards.condition:q} -r {wildcards.replicate:q} -i {input.predictions:q} -g {input.genome:q} -o {output:q}"

rule concatDeepRibo:
    input:
        gffs=lambda wildcards: expand("deepribo/{{condition}}-{replicate}.deepribo.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        script=str(SCRIPTS / "concatenate_gff.py")
    output:
        "tracks/{condition}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} {input.gffs:q} -o {output:q}"

rule allDeepRibo:
    input:
        merged_gff=expand("tracks/{condition}.deepribo.gff", condition=conditions),
        script=str(SCRIPTS / "concatenate_gff.py")
    output:
        "tracks/deepribo_all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} {input.merged_gff:q} -o {output:q}"

rule filterDeepRibo:
    input:
        ingff="tracks/deepribo_all.gff",
        annotation=rules.checkAnnotation.output,
        script=str(SCRIPTS / "merge_duplicates_deepribo.py"),
        script_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        merged="tracks/deepribo_merged.gff",
        plus="tracks/deepribo_merged_plus.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        python3 {input.script:q} \
            -i {input.ingff:q} \
            -o {output.merged:q} \
            --plus-output {output.plus:q} \
            -a {input.annotation:q}
        """


rule createExcelSummaryDeepRibo:
    input:
        total="readcounts/bam_mapped_reads.txt",
        reads="readcounts/deepribo_annotation.gff",
        genome="genomes/genome.fa",
        script=str(SCRIPTS / "generate_excel_deepribo.py"),
        script_deps=[
            str(SCRIPTS / "excel_utils.py"),
            str(SCRIPTS / "gff_utils.py"),
        ]
    output:
        "auxiliary/predictions_deepribo.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -t {input.total:q} -r {input.reads:q} -g {input.genome:q} -o {output:q}"
