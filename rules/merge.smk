rule mergeConditions:
    input:
        ribotish="tracks/{condition}.ribotish.gff",
        reparation="tracks/{condition}.reparation.gff"
    output:
        report("tracks/{condition, [a-zA-Z]+}.merged.gff", caption="../report/novelmerged.rst", category="Novel ORFs")
    conda:
        "../envs/bedtools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input.ribotish} > {output}.unsorted; cat {input.reparation} >> {output}.unsorted; bedtools sort -i {output}.unsorted > {output};"

rule noverlap:
    input:
        mergedGff="tracks/{condition}.merged.gff"
    output:
         report("tracks/{condition}.filtered.gff", caption="../report/novelfiltered.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/noverlapper.py -i {input.mergedGff} -o {output}"

rule mergeAll:
    input:
        mergedGff=expand("tracks/{condition}.filtered.gff", zip, condition=samples["condition"])
    output:
        report("tracks/all.gff", caption="../report/novelall.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input.mergedGff} -o {output}"

rule filterAll:
    input:
        "tracks/all.gff"
    output:
        report("tracks/filtered.gff", caption="../report/novelallfiltered.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/noverlapper.py -i {input} -o {output}"

rule newAnnotation:
    input:
        newOrfs="tracks/filtered.gff",
        currentAnnotation="xtail/longest_protein_coding_transcripts.gtf"
    output:
        report("xtail/newAnnotation.gff", caption="../report/novelannotation.rst", category="Novel ORFs")
    conda:
        "../envs/mergetools.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/concatGFF.py {input.newOrfs} {input.currentAnnotation} -o xtail/tmp.gff; SPtools/scripts/noverlapper.py -i xtail/tmp.gff -o {output};"
