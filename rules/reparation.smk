
rule uniprotDBRetrieve:
    input:
        storage.ftp("ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz")
    output:
        "uniprotDB/uniprot_sprot.fasta"
    threads: 1
    shell:
        """
        mkdir -p uniprotDB
        mv {input} {output}.gz
        gunzip {output}.gz
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
        reparation_instances=1
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
        "mkdir -p tracks; HRIBO/scripts/create_reparation_gff.py -c {wildcards.condition} -r {wildcards.replicate} -i {input} -o {output}"

rule concatReparation:
    input:
        lambda wildcards: expand("reparation/{{condition}}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z0-9]+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; HRIBO/scripts/concatenate_gff.py {input} -o {output}"
