rule xtailclassicnormalize:
    input:
        bam="maplink/{method}-{condition}-{replicate}.bam",
        annotation="annotation/annotation.gtf"
    output:
        xtailclassic="xtailclassic/{method}-{condition}-{replicate}.csv"
    conda:
        "../envs/xtailcounts.yaml"
    threads: 1
    shell: ("mkdir -p xtailclassic; python RPF_count_CDS.py {input.bam} {input.annotation} > {output}")


#rule xtailclassic:
#    input:
#        "uORFs/norm_CDS_reads.csv"
#    output:
#        table=report("xtail/xtail_cds.csv", caption="../report/xtail_cds_fc.rst", category="CDS"),
#        fcplot="xtailclassic/xtail_cds_fc.pdf",
#        rplot="xtailclassic/xtail_cds_r.pdf"
#    conda:
#        "../envs/xtail.yaml"
#    threads: 1
#    shell: ("mkdir -p xtailclassic; SPtools/scripts/xtail_classic.R -t SPtools/samples.tsv -r {input} -x {output.table} -f {output.fcplot} -p {output.rplot};")
