from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

def read_parameters(filename, idx):
    try:
        line = ""
        with open(filename, "r") as f:
            line = f.readline()[:-1].split(",")
        return line[idx]
    except FileNotFoundError:
        return "failed"

rule deepriboGetModel:
    input:
        HTTP.remote("github.com/Biobix/DeepRibo/blob/master/models/DeepRibo_model_v1.pt?raw=true", keep_local=True)
    output:
        "deepribo/DeepRibo_model_v1.pt"
    run:
        shell("mkdir -p deepribo; mv {input} deepribo/DeepRibo_model_v1.pt")


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
    shell:
        "mkdir -p coverage_deepribo; HRIBO/scripts/coverage_deepribo.py --alignment_file {input.bam} --output_file_prefix coverage_deepribo/{wildcards.condition}-{wildcards.replicate}"

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
    shell:
        """
        mkdir -p coverage_deepribo
        bedtools genomecov -bg -ibam {input.bam} -strand + > {output.covfwd}
        bedtools genomecov -bg -ibam {input.bam} -strand - > {output.covrev}
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
    singularity:
        "docker://gelhausr/deepribo:minimal"
    threads: 1
    shell:
        """
        mkdir -p deepribo/{wildcards.condition}-{wildcards.replicate}/0/;
        mkdir -p deepribo/{wildcards.condition}-{wildcards.replicate}/1/;
        DataParser.py {input.covS} {input.covAS} {input.asiteS} {input.asiteAS} {input.genome} deepribo/{wildcards.condition}-{wildcards.replicate} -g {input.annotation}
        """

rule parameterEstimation:
    input:
        "deepribo/{condition}-{replicate}/data_list.csv"
    output:
        "deepribo/{condition}-{replicate}/parameters.txt"
    singularity:
        "docker://gelhausr/deepribo:minimal"
    threads: 1
    shell:
        "mkdir -p deepribo; Rscript HRIBO/scripts/parameter_estimation.R -f {input} -o {output}"

rule predictDeepRibo:
    input:
        model= "deepribo/DeepRibo_model_v1.pt",
        data= "deepribo/{condition}-{replicate}/data_list.csv",
        parameter= "deepribo/{condition}-{replicate}/parameters.txt"
    output:
        "deepribo/{condition}-{replicate}/predictions.csv"
    singularity:
        "docker://gelhausr/deepribo:minimal"
    threads: 10
    params:
        rpkm= lambda wildcards, input: read_parameters(input[2], 0),
        cov= lambda wildcards, input: read_parameters(input[2], 1)
    shell:
        """
        mkdir -p deepribo;
        DeepRibo.py predict deepribo/ --pred_data {wildcards.condition}-{wildcards.replicate}/ -r {params.rpkm} -c {params.cov} --model {input.model} --dest {output} --num_workers {threads}
        """

rule deepriboGFF:
    input:
        "deepribo/{condition}-{replicate}/predictions.csv"
    output:
        "deepribo/{condition}-{replicate,\d+}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/create_deepribo_gff.py -c {wildcards.condition} -r {wildcards.replicate} -i {input} -o {output}"

rule concatDeepRibo:
    input:
        lambda wildcards: expand("deepribo/{{condition}}-{replicate}.deepribo.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition}.deepribo.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/concatenate_gff.py {input} -o {output}"

rule allDeepRibo:
    input:
        merged_gff=expand("tracks/{condition}.deepribo.gff", zip, condition=set(samples["condition"]))
    output:
        "tracks/deepribo_all.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/concatenate_gff.py {input.merged_gff} -o {output}"

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
        "mkdir -p tracks; HRIBO/scripts/merge_duplicates_deepribo.py -i {input.ingff} -o {output.merged} -a {input.annotation}"


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
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel_deepribo.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"

rule newAnnotationDeepRibo:
    input:
        reparation_orfs="tracks/reparation_annotated.gff",
        deepribo_orfs="tracks/deepribo_merged_plus.gff",
        currentAnnotation=rules.checkAnnotation.output
    output:
        "tracks/totalAnnotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p tracks;
        HRIBO/scripts/concatenate_gff.py {input.deepribo_orfs} {input.reparation_orfs} {input.currentAnnotation} -o {output}
        """
