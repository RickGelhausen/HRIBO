
UNIPROT_SPROT_RELEASE = "2026_02"
UNIPROT_SPROT_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/"
    "knowledgebase/complete/uniprot_sprot.fasta.gz"
)
UNIPROT_SPROT_SHA256 = (
    "b774748a050fd3de0bf2ad49b359ed59a2c2b02c89df3fe0679fafa446751794"
)
UNIPROT_SPROT_SIZE = 93706469
UNIPROT_SPROT_ARCHIVE_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/"
    "release-2026_02/knowledgebase/uniprot_sprot-only2026_02.tar.gz"
)
UNIPROT_SPROT_ARCHIVE_MD5 = "e2ddaeb6739a3771e9748e7db5573c6a"
UNIPROT_SPROT_ARCHIVE_SIZE = 1737439496
UNIPROT_SPROT_ARCHIVE_MEMBER = "uniprot_sprot.fasta.gz"
# BioContainers publishes this historical dependency stack for Linux/amd64.
REPARATION_CONTAINER = (
    "docker://quay.io/biocontainers/reparation_blast@sha256:"
    "6852b3b69b532039a5d674115b9cbfd2953f6479cf187bcc769e2d899fcbc288"
)

# UniProtKB/Swiss-Prot release 2026_02 is 94 MB in the normal download area,
# while its immutable previous-release bundle is 1.7 GB. The checksum-locked
# small file is preferred; after UniProt advances current_release, the official
# archived bundle is verified and used as the reproducible fallback.


rule prepareReparationAnnotation:
    input:
        annotation=rules.checkAnnotation.output,
        adapter=str(SCRIPTS / "prepare_reparation_annotation.py"),
        adapter_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        "reparation/annotation.gtf"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        "python3 {input.adapter:q} -a {input.annotation:q} -o {output:q}"


rule uniprotDBRetrieve:
    input:
        fetcher=str(SCRIPTS / "fetch_verified.py")
    output:
        "uniprotDB/uniprot_sprot.fasta"
    params:
        url=UNIPROT_SPROT_URL,
        sha256=UNIPROT_SPROT_SHA256,
        size=UNIPROT_SPROT_SIZE,
        release=UNIPROT_SPROT_RELEASE,
        archive_url=UNIPROT_SPROT_ARCHIVE_URL,
        archive_md5=UNIPROT_SPROT_ARCHIVE_MD5,
        archive_size=UNIPROT_SPROT_ARCHIVE_SIZE,
        archive_member=UNIPROT_SPROT_ARCHIVE_MEMBER
    conda:
        "../envs/download.yaml"
    threads: 1
    retries: 3
    resources:
        mem_mb=1000,
        disk_mb=3000,
        runtime=60
    log:
        "logs/uniprotDBRetrieve.log"
    shell:
        """
        compact_status=0
        python3 {input.fetcher:q} \
            --url {params.url:q} \
            --sha256 {params.sha256:q} \
            --size {params.size} \
            --decompress-gzip \
            --output {output:q} > {log:q} 2>&1 || compact_status=$?
        if [ "$compact_status" -eq 42 ]
        then
            echo "The compact mirror no longer contains UniProt release 2026_02; using its archived bundle." >> {log:q}
            python3 {input.fetcher:q} \
                --url {params.archive_url:q} \
                --md5 {params.archive_md5:q} \
                --size {params.archive_size:q} \
                --tar-member {params.archive_member:q} \
                --decompress-gzip \
                --output {output:q} >> {log:q} 2>&1
        elif [ "$compact_status" -ne 0 ]
        then
            exit "$compact_status"
        fi
        """


rule reparation:
    input:
        genome=rules.retrieveGenome.output,
        gtf=rules.prepareReparationAnnotation.output,
        db="uniprotDB/uniprot_sprot.fasta",
        bam="maplink/RIBO-{condition}-{replicate}.bam",
        bamindex="maplink/RIBO-{condition}-{replicate}.bam.bai"
    output:
        orfs="reparation/{condition}-{replicate}/Predicted_ORFs.txt",
        metagene="reparation/{condition}-{replicate}/metagene_profile.pdf",
        roc="reparation/{condition}-{replicate}/PR_and_ROC_curve.pdf",
        psite="reparation/{condition}-{replicate}/p_site_offset.png",
        scurve="reparation/{condition}-{replicate}/S_Curve.pdf"
    container:
        REPARATION_CONTAINER
    threads: 12
    resources:
        reparation_instances=1,
        mem_mb=30000,
        runtime=240
    params:
        prefix=lambda wildcards, output: os.path.dirname(output.orfs),
        temporary=lambda wildcards, output: os.path.join(
            os.path.dirname(output.orfs), "tmp"
        )
    log:
        "logs/{condition}-{replicate}_reparation.log"
    shell:
        """
        exec > {log:q} 2>&1
        mkdir -p {params.prefix:q} {params.temporary:q}
        reparation.pl -bam {input.bam:q} -g {input.genome:q} -gtf {input.gtf:q} -db {input.db:q} -wdir {params.prefix:q} -threads {threads}
        """

rule reparationGFF:
    input:
        orfs="reparation/{condition}-{replicate}/Predicted_ORFs.txt",
        script=str(SCRIPTS / "create_reparation_gff.py"),
        script_deps=[str(SCRIPTS / "gff_utils.py")]
    output:
        "reparation/{condition}-{replicate}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} -c {wildcards.condition:q} -r {wildcards.replicate:q} -i {input.orfs:q} -o {output:q}"

rule concatReparation:
    input:
        gffs=lambda wildcards: expand("reparation/{{condition}}-{replicate}.reparation.gff", zip, replicate=samples.loc[(samples["method"] == "RIBO") & (samples["condition"] == wildcards.condition), "replicate"]),
        script=str(SCRIPTS / "concatenate_gff.py")
    output:
        "tracks/{condition}.reparation.gff"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} {input.gffs:q} -o {output:q}"
