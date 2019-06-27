from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

def get_dbs():
    group = config["taxonomy"]
    dbs=[]
    if group == "Eukarya":
        dbs=["rfam-5.8s-database-id98","rfam-5s-database-id98","silva-euk-18s-id95","silva-euk-28s-id98"]
    elif group == "Bacteria":
        dbs=["rfam-5s-database-id98","silva-bac-23s-id98","silva-bac-16s-id90"]
    elif group == "Archea":
        dbs=["rfam-5s-database-id98","silva-arc-16s-id95","silva-arc-23s-id98"]
    else:
        dbs=[]
    return dbs

def get_indexfiles ():
    dbs=get_dbs()
    indexstring=""
    for rrnadb in dbs:
        dbstring = "./rRNA_databases/" + rrnadb + ".fasta" + ",./index/rRNA/" + rrnadb + ":"
        indexstring = dbstring + indexstring
    return str(indexstring)

rule rrnaretrieve:
    input:
        HTTP.remote("https://github.com/biocore/sortmerna/raw/master/rRNA_databases/{rrnadb}.fasta",keep_local=True,allow_redirects=True)
    output:
        "rRNA_databases/{rrnadb}.fasta"
    threads: 1
    run:
        outputName = os.path.basename(input[0])
        shell("mkdir -p rRNA_databases; mv {input} rRNA_databases/{outputName}")

rule rrnaindex:
    input:
        ["rRNA_databases/{rrnadb}.fasta".format(rrnadb=rrnadb) for rrnadb in get_dbs()]
    output:
        ["index/rRNA/{rrnadb}.bursttrie_0.dat".format(rrnadb=rrnadb) for rrnadb in get_dbs()]
    conda:
        "../envs/sortmerna.yaml"
    params:
        dbstring = get_indexfiles()
    threads: 1
    shell:
        "mkdir -p index/rRNA; indexdb_rna --ref {params.dbstring}"

rule rrnafilter:
    input:
        "trimmed/{method}-{condition}-{replicate}.fastq",
        rules.rrnaindex.output
    output:
        norrna="norRNA/{method}-{condition}-{replicate}.fastq",
        rrna="norRNA/rRNA/reject/{method}-{condition}-{replicate}.fastq"
    conda:
        "../envs/sortmerna.yaml"
    params:
        prefix=lambda wildcards, output: (os.path.splitext(output.norrna)[0]),
        rejectprefix=lambda wildcards, output: (os.path.splitext(output.rrna)[0]),
        dbstring = get_indexfiles()
    threads: 20
    shell:
        "mkdir -p norRNA; mkdir -p norRNA/rRNA; sortmerna -a {threads} --ref {params.dbstring} --reads {input[0]} --num_alignments 1 --fastx --log --aligned {params.rejectprefix} --other {params.prefix} 2> /dev/null"

