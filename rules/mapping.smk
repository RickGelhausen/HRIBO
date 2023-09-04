rule genomeSegemehlIndex:
    input:
        genome=rules.retrieveGenome.output
    output:
        index=temp("genomeSegemehlIndex/genome.idx")
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    log:
        "logs/genomeIndex.log"
    shell:
        "mkdir -p genomeSegemehlIndex; echo \"Computing Segemehl index\"; segemehl.x --threads {threads} -x {output.index} -d {input.genome} 2> {log}"


def get_fastq_files(wildcards):
    output = dict()
    row = samples[(samples['method'] == wildcards.method) & (samples['condition'] == wildcards.condition) & (samples['replicate'] == wildcards.replicate)]
    if pd.isna(row['fastqFile2'].iloc[0]):
        output["fastq"] = f"trimmed/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}.fastq"
    else:
        output["fastq1"] = f"trimmed/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}_q.fastq"
        output["fastq2"] = f"trimmed/{wildcards.method}-{wildcards.condition}-{wildcards.replicate}_p.fastq"
    return output

rule map:
    input:
        unpack(get_fastq_files),
        genome=rules.retrieveGenome.output,
        genomeSegemehlIndex="genomeSegemehlIndex/genome.idx"
    output:
        sammulti=temp("sammulti/{method}-{condition}-{replicate}.sam")
    conda:
        "../envs/segemehl.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.dirname(output[0])),
        fastq=lambda wildcards, input: "-q %s" % (input.fastq) if len(input) == 3 else "-q %s -p %s" % (input.fastq1, input.fastq2)
    log:
        "logs/{method}-{condition}-{replicate}_segemehl.log"
    shell:
        """
        mkdir -p sammulti; segemehl.x -e -d {input.genome} -i {input.genomeSegemehlIndex} {params.fastq} --threads {threads} -o {output.sammulti} 2> {log}
        """

rule samuniq:
    input:
        sammulti="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        sam=temp("sam/{method}-{condition}-{replicate}.rawsam"),
        unmapped=temp("sammulti/{method}-{condition}-{replicate}.sam.unmapped")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    shell:
        """
        set +e
        mkdir -p sam
        awk '$2 == "4"' {input.sammulti} > {input.sammulti}.unmapped
        gawk -i inplace '$2 != "4"' {input.sammulti}
        samtools view -H <(cat {input.sammulti}) | grep '@HD' > {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V >> {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@RG' >> {output.sam}
        samtools view -H <(cat {input.sammulti}) | grep '@PG' >> {output.sam}
        cat {input.sammulti} |grep -v '^@' | grep -w 'NH:i:1' >> {output.sam}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

rule samstrandswap:
    input:
        sam="sam/{method}-{condition}-{replicate}.rawsam"
    output:
        sam=temp("sam/{method}-{condition}-{replicate}.sam")
    threads: 1
    params:
         method=lambda wildcards: wildcards.method
    shell: "if [ \"{params.method}\" == \"NOTSET\" ]; then HRIBO/scripts/sam_strand_inverter.py --sam_in_filepath={input.sam} --sam_out_filepath={output.sam}; else cp {input.sam} {output.sam}; fi"

rule sammultitobam:
    input:
        sam="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        temp("bammulti/{method}-{condition}-{replicate}.bam")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    shell:
        "mkdir -p bammulti; samtools view -@ {threads} -bh {input.sam} | samtools sort -@ {threads} -o {output} -O bam"

rule samtobam:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        temp("rRNAbam/{method}-{condition}-{replicate}.bam")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    shell:
        "mkdir -p rRNAbam; samtools view -@ {threads} -bh {input.sam} | samtools sort -@ {threads} -o {output} -O bam"

rule maplink:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "maplink/{method}-{condition}-{replicate}.bam"
    params:
        inlink=lambda wildcards, input:(os.getcwd() + "/" + str(input)),
        outlink=lambda wildcards, output:(os.getcwd() + "/" + str(output))
    threads: 1
    shell:
        "mkdir -p maplink; ln -s {params.inlink} {params.outlink}"
