from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule uniprotDBRetrieve:
    input:
        FTP.remote("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",keep_local=True,allow_redirects=True)
    output:
        "uniprotDB/uniprot_sprot.fasta"
    threads: 1
    run:
        outputName = os.path.basename(input[0])
        shell("mkdir -p uniprotDB; mv {input} uniprotDB/{outputName}; gunzip uniprotDB/{outputName}")

rule reparation:
    input:
        sam="sam/RIBO-{condition}-{replicate}.sam",
        genome=rules.retrieveGenome.output,
        gtf=rules.retrieveAnnotation.output,
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
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output.orfs))
    log:
        "logs/{condition, [a-zA-Z]+}-{replicate,\d+}_reparation.log"
    shell:
        "mkdir -p reparation; if [ uniprotDB/uniprot_sprot.fasta.bak does not exist ]; then cp -p uniprotDB/uniprot_sprot.fasta uniprotDB/uniprot_sprot.fasta.bak; fi; mkdir -p {params.prefix}/tmp; reparation.pl -sam {input.sam} -g {input.genome} -gtf {input.gtf} -db {input.db} -out {params.prefix} -bam {input.bam} -threads {threads}; if [ uniprotDB/uniprot_sprot.fasta does not exist ]; then cp -p uniprotDB/uniprot_sprot.fasta.bak uniprotDB/uniprot_sprot.fasta; fi;"

rule reparationGFF:
    input:
        "reparation/{condition}-{replicate}/Predicted_ORFs.txt"
    output:
        "reparation/{condition, [a-zA-Z]+}-{replicate,\d+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/reparationGFF.py -i {input} -o {output}"

rule concatReparation:
    input:
        lambda wildcards: expand("reparation/{{condition}}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"])
    output:
        "tracks/{condition, [a-zA-Z]+}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input} -o {output}"

rule reparationreport:
    input:
        metagene="reparation/{condition}-{replicate}/metagene_profile.pdf",
        roc="reparation/{condition}-{replicate}/PR_and_ROC_curve.pdf",
        psite="reparation/{condition}-{replicate}/p_site_offset.png",
        scurve="reparation/{condition}-{replicate}/S_Curve.pdf"
    output:
        metagene="figures/{condition}-{replicate}_metagene.jpg",
        roc="figures/{condition}-{replicate}_roc.jpg",
        psite="figures/{condition}-{replicate}_psite.png",
        scurve="figures/{condition}-{replicate}_scurve.jpg"
        #metagene=report("figures/{condition}-{replicate}_metagene.jpg", caption="../report/reparation_metagene.rst", category="Novel ORFs - Reparation"),
        #roc=report("figures/{condition}-{replicate}_roc.jpg", caption="../report/reparation_roc.rst", category="Novel ORFs - Reparation"),
        #psite=report("figures/{condition}-{replicate}_psite.png", caption="../report/reparation_psite.rst", category="Novel ORFs - Reparation"),
        #scurve=report("figures/{condition}-{replicate}_scurve.jpg", caption="../report/reparation_scurve.rst", category="Novel ORFs - Reparation")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input.metagene} -quality 100 -flatten -sharpen 0x1.0 {output.metagene}; convert -density 150 -trim {input.roc} -quality 100 -flatten -sharpen 0x1.0 {output.roc}; cp {input.psite} {output.psite}; convert -density 150 -trim {input.scurve} -quality 100 -flatten -sharpen 0x1.0 {output.scurve};")
