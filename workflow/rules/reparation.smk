
rule uniprotDBRetrieve:
    output:
        "uniprotDB/uniprot_sprot.fasta"
    params:
        url="https://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
    conda:
        "../envs/download.yaml"
    threads: 1
    retries: 3
    resources:
        mem_mb=1000,
        runtime=60
    log:
        "logs/uniprotDBRetrieve.log"
    shell:
        """
        curl -sSL --fail --retry 3 --retry-delay 5 {params.url} -o {output}.gz 2> {log}
        gunzip -f {output}.gz 2>> {log}
        """


rule reparation:
    input:
        genome=rules.retrieveGenome.output,
        gtf=rules.checkAnnotation.output,
        db="uniprotDB/uniprot_sprot.fasta",
        bam="bam/RIBO-{condition}-{replicate}.bam",
        bamindex="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        orfs="reparation/{condition}-{replicate}/Predicted_ORFs.txt",
        metagene="reparation/{condition}-{replicate}/metagene_profile.pdf",
        roc="reparation/{condition}-{replicate}/PR_and_ROC_curve.pdf",
        psite="reparation/{condition}-{replicate}/p_site_offset.png",
        scurve="reparation/{condition}-{replicate}/S_Curve.pdf"
    conda:
        "../envs/reparation.yaml"
    threads: 12
    resources:
        reparation_instances=1,
        mem_mb=30000,
        runtime=240
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output.orfs))
    log:
        r"logs/{condition, [a-zA-Z]+}-{replicate,\d+}_reparation.log"
    shell:
        """
        mkdir -p {params.prefix}
        mkdir -p {params.prefix}/tmp
        reparation.pl -bam {input.bam} -g {input.genome} -gtf {input.gtf} -db {input.db} -out {params.prefix} -threads {threads}
        """

rule reparationGFF:
    input:
        "reparation/{condition}-{replicate}/Predicted_ORFs.txt"
    output:
        r"reparation/{condition, [a-zA-Z0-9]+}-{replicate,\d+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/create_reparation_gff.py -c {wildcards.condition} -r {wildcards.replicate} -i {input} -o {output}"

rule concatReparation:
    input:
        lambda wildcards: expand("reparation/{{condition}}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z0-9]+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/concatenate_gff.py {input} -o {output}"
