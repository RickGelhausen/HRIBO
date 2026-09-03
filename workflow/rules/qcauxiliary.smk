rule fastqcunique:
    input:
        sam="sam/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/4unique/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/4unique/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/map/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    resources:
        mem_mb=30000,
        runtime=60
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "fastqc -o qc/4unique -t {threads} -f sam_mapped {input.sam:q}; mv qc/4unique/{params.prefix:q}_fastqc.html {output.html:q}; mv qc/4unique/{params.prefix:q}_fastqc.zip {output.zip:q}"

rule fastqcmulti:
    input:
        sam="sammulti/{method}-{condition}-{replicate}.sam"
    output:
        html="qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.html",
        zip="qc/3mapped/{method}-{condition}-{replicate}-map_fastqc.zip",
        #report("qc/sammulti/{method}-{condition}-{replicate}-map.html", caption="../report/fastqcmapped.rst", category="Mapped reads")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    resources:
        mem_mb=40000,
        runtime=60
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.sam))[0])
    shell:
        "fastqc -o qc/3mapped -t {threads} -f sam_mapped {input.sam:q}; mv qc/3mapped/{params.prefix:q}_fastqc.html {output.html:q}; mv qc/3mapped/{params.prefix:q}_fastqc.zip {output.zip:q}"

rule fastqcrrnafilter:
    input:
        reads="bam/{method}-{condition}-{replicate}.bam"
    output:
        html="qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.html",
        zip="qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA_fastqc.zip"
        #report("qc/5removedrRNA/{method}-{condition}-{replicate}-norRNA.html", caption="../report/fastqcnorRNA.rst", category="Removing hits mapping to rRNA")
    conda:
        "../envs/fastqc.yaml"
    threads: 8
    resources:
        mem_mb=30000,
        runtime=60
    params:
        prefix=lambda wildcards, input: (os.path.splitext(os.path.basename(input.reads))[0])
    shell:
        "fastqc -o qc/5removedrRNA -t {threads} {input:q}; mv qc/5removedrRNA/{params.prefix:q}_fastqc.html {output.html:q}; mv qc/5removedrRNA/{params.prefix:q}_fastqc.zip {output.zip:q}"

rule featurescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/all/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    resources:
        mem_mb=40000,
        runtime=60
    shell:
        """
        column3=$(cut -f3 {input.annotation:q} | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "gene" ]];
        then
            featureCounts -T {threads} -t gene -g ID -a {input.annotation:q} -o {output.txt:q} {input.bam:q};
        else
            touch {output.txt:q};
        fi
        """

rule trnafeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/trnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    resources:
        mem_mb=30000,
        runtime=60
    shell:
        """
        column3=$(cut -f3 {input.annotation:q} | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "tRNA" ]];
        then
            featureCounts -T {threads} -t tRNA -g ID -a {input.annotation:q} -o {output.txt:q} {input.bam:q};
        else
            touch {output.txt:q};
        fi
        """

rule norrnafeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainall/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    resources:
        mem_mb=30000,
        runtime=60
    shell:
        """
        column3=$(cut -f3 {input.annotation:q} | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "rRNA" ]];
        then
            featureCounts -T {threads} -t rRNA -g ID -a {input.annotation:q} -o {output.txt:q} {input.bam:q};
        else
            touch {output.txt:q};
        fi
        """

rule rrnatotalfeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="bammulti/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainallaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    resources:
        mem_mb=30000,
        runtime=60
    shell:
        """
        column3=$(cut -f3 {input.annotation:q} | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "rRNA" ]];
        then
            featureCounts -T {threads} -t rRNA -g ID -a {input.annotation:q} -o {output.txt:q} {input.bam:q};
        else
            touch {output.txt:q};
        fi
        """
        
rule rrnauniquefeaturescounts:
    input:
        annotation={rules.unambigousAnnotation.output},
        bam="rRNAbam/{method}-{condition}-{replicate}.bam"
    output:
        txt="qc/rrnainuniquelyaligned/{method}-{condition}-{replicate}.txt",
    conda:
        "../envs/subread.yaml"
    threads: 8
    resources:
        mem_mb=30000,
        runtime=60
    shell:
        """
        column3=$(cut -f3 {input.annotation:q} | sort | uniq)
        if [[ " ${{column3[@]}} " =~ "rRNA" ]];
        then
            featureCounts -T {threads} -t rRNA -g ID -a {input.annotation:q} -o {output.txt:q} {input.bam:q};
        else
            touch {output.txt:q};
        fi
        """

rule coveragedepth:
    input:
        "bam/{method}-{condition}-{replicate}.bam"
    output:
        "coverage/{method}-{condition}-{replicate}.bed"
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "bedtools genomecov -ibam {input:q} -bg > {output:q}"
