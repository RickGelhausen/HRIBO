rule mergeConditions:
    input:
        ribotish="tracks/{condition, [a-zA-Z]+}.ribotish.gff",
        reparation="tracks/{condition, [a-zA-Z]+}.reparation.gff"
    output:
        "tracks/{condition, [a-zA-Z]+}.merged.gff"
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.ribotish} > {output}.unsorted; cat {input.reparation} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};" 

#rule mergeAll:
#    input:
#        orfs=expand("tracks/{sample.condition}.merged.gff", sample=samples.itertuples())
#    output:
#        "tracks/novelORFs.gff"
#    conda:
#        "envs/merge.py"
#    threads: 1
#    shell:
#        "mkdir -p tracks; merge_orfs.py -i {input.orfs} -o {output};"

#rule newAnnotation:
#    input:
#        newOrfs="tracks/novelORFs.gff"
#        currentAnnoation="xtail/longest_protein_coding_transcripts.gtf" 
#    output:
#        "tracks/newAnnotation.gff"
#    conda:
#        "envs/merge.py"
#    threads: 1
#    shell:
#        "mkdir -p tracks; merge_orfs.py -i {input.currentAnnotation} {input.newOrfs} -o {output};"


