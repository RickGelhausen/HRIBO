from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule uniprotDBRetrieve:
    input:
        FTP.remote("ftp://biftp.informatik.uni-freiburg.de/pub/HRIBO/uniprot_sprot.fasta.gz",keep_local=True,allow_redirects=True)
    output:
        "uniprotDB/uniprot_sprot.fasta"
    threads: 1
    run:
        outputName = os.path.basename(input[0])
        shell("mkdir -p uniprotDB; mv {input} uniprotDB/{outputName}; gunzip uniprotDB/{outputName}")

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
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output.orfs))
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_reparation.log"
    shell:
        "mkdir -p reparation; if [ uniprotDB/uniprot_sprot.fasta.bak does not exist ]; then cp -p uniprotDB/uniprot_sprot.fasta uniprotDB/uniprot_sprot.fasta.bak; fi; mkdir -p {params.prefix}/tmp; reparation.pl -bam {input.bam} -g {input.genome} -gtf {input.gtf} -db {input.db} -out {params.prefix} -threads {threads}; if [ uniprotDB/uniprot_sprot.fasta does not exist ]; then cp -p uniprotDB/uniprot_sprot.fasta.bak uniprotDB/uniprot_sprot.fasta; fi;"

rule reparationGFF:
    input:
        "reparation/{condition}-{replicate}/Predicted_ORFs.txt"
    output:
        "reparation/{condition, [a-zA-Z0-9]+}-{replicate,\d+}.reparation.gff"
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
