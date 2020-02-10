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
        annotation= rules.retrieveAnnotation.output
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
        annotation=rules.retrieveAnnotation.output
    output:
        "tracks/deepribo_merged.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/merge_duplicates_deepribo.py -i {input.ingff} -o {output} -a {input.annotation}"

rule generateAnnotationDeepRiboReadCounts:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        annotation="tracks/deepribo_merged.gff"
    output:
        "auxiliary/annotation_deepribo_reads.raw"
    conda:
        "../envs/subread.yaml"
    threads: 5
    shell:
        """
        mkdir -p auxiliary
        featureCounts -F GTF -s 1 -g ID -O -t CDS -M --fraction -a {input.annotation} {input.bam} -T {threads} -o auxiliary/annotation_deepribo_reads.raw.tmp
        cat auxiliary/annotation_deepribo_reads.raw.tmp | sed 1,2d | awk -v var=CDS -FS'\\t' '{{print $0"\\t"var}}' >> {output}
        rm auxiliary/annotation_deepribo_reads.raw.tmp
        """

rule mapDeepRiboReads:
    input:
        reads="auxiliary/annotation_deepribo_reads.raw",
        annotation="tracks/deepribo_merged.gff"
    output:
        "auxiliary/deepribo_annotation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        """
        mkdir -p auxiliary; HRIBO/scripts/map_reads_to_annotation.py -i {input.reads} -a {input.annotation} -o {output}
        """

rule mappedReadsDeepRibo:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        bamindex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        mapped="auxiliary/deepribo_sum_mapped_reads.txt",
        length="auxiliary/deepribo_average_read_lengths.txt"
    conda:
        "../envs/pytools.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/total_mapped_reads.py -b {input.bam} -m {output.mapped} -l {output.length}"


rule createExcelSummaryDeepRibo:
    input:
        total="auxiliary/deepribo_sum_mapped_reads.txt",
        reads="auxiliary/deepribo_annotation.gff",
        genome="genomes/genome.fa"
    output:
        "auxiliary/predictions_deepribo.xlsx"
    conda:
        "../envs/excel.yaml"
    threads: 1
    shell:
        "mkdir -p auxiliary; HRIBO/scripts/generate_excel_deepribo.py -t {input.total} -r {input.reads} -g {input.genome} -o {output}"
