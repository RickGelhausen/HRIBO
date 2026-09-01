from pathlib import Path

rule genomeSamToolsIndex:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.fa.fai"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    shell:
        "samtools faidx {rules.retrieveGenome.output}"

rule genomeSize:
    input:
        rules.genomeSamToolsIndex.output
    output:
        "genomes/sizes.genome"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    log: "logs/genomeSamToolsIndex.log"
    shell:
        "cut -f1,2 {input[0]} > genomes/sizes.genome"

rule reversecomplementGenome:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.rev.fa"
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/reverse_complement.py --input_fasta_filepath genomes/genome.fa --output_fasta_filepath genomes/genome.rev.fa"

rule startCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        report("tracks/potentialStartCodons.gff", caption="../report/startCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/motif_to_gff.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string ATG --output_gff3_filepath {output}"

rule alternativeStartCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output,
        script=str(SCRIPTS / "motif_to_gff.py")
    output:
        report("tracks/potentialAlternativeStartCodons.gff", caption="../report/startCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    params:
        motifs=",".join(str(codon).upper() for codon in CODONS)
    shell:
        "python3 {input.script:q} --input_genome_fasta_filepath {input.fwd:q} --input_reverse_genome_fasta_filepath {input.rev:q} --motif_string {params.motifs:q} --output_gff3_filepath {output:q}"


rule stopCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        report("tracks/potentialStopCodons.gff", caption="../report/stopCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/motif_to_gff.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string TAG,TGA,TAA --output_gff3_filepath {output}"

rule rbsTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        report("tracks/potentialRibosomeBindingSite.gff", caption="../report/rbsTrack.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "{SCRIPTS}/motif_to_gff.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string AAGG --output_gff3_filepath {output}"


rule bamindex:
    input:
        rules.maplink.output,
        rules.genomeSize.output
    output:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=60
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "samtools index -@ {threads} maplink/{params.prefix}"

rule totalmappedbamindex:
    input:
        rules.sammultitobam.output,
        rules.genomeSize.output
    output:
        temp("bammulti/{method}-{condition}-{replicate}.bam.bai")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=60
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "samtools index -@ {threads} bammulti/{params.prefix}"

rule uniquemappedbamindex:
    input:
        rules.samtobam.output,
        rules.genomeSize.output
    output:
        temp("rRNAbam/{method}-{condition}-{replicate}.bam.bai")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=60
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "samtools index -@ {threads} rRNAbam/{params.prefix}"

# Coverage tracks
#
# All track flavours are produced by the same rule; they differ only in which
# alignments they read and which mapping style scripts/mapping.py applies. The
# {mapping} wildcard is the directory prefix and indexes TRACK_SOURCES, defined
# in common.smk.

wildcard_constraints:
    mapping="|".join(TRACK_SOURCES),
    norm="|".join(TRACK_NORMALIZATIONS),
    strand="|".join(TRACK_STRANDS)


rule coverageTracks:
    input:
        bam=track_bam,
        bamIndex=track_bam_index,
        stats=track_stats,
        genomeSize=rules.genomeSize.output
    output:
        # mapping.py writes all three normalisations and both strands in one pass.
        wig=expand(
            "{{mapping}}tracks/{norm}/{{method}}-{{condition}}-{{replicate}}.{norm}.{strand}.wig",
            norm=TRACK_NORMALIZATIONS,
            strand=TRACK_STRANDS
        )
    conda:
        "../envs/coverage.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    params:
        mapping_style=track_mapping_style,
        library=library_name
    log:
        "logs/{mapping}tracks_{method}-{condition}-{replicate}.log"
    shell:
        """
        {SCRIPTS}/mapping.py \
            --mapping_style {params.mapping_style} \
            --bam_path {input.bam} \
            --wiggle_file_path {wildcards.mapping}tracks/ \
            --no_of_aligned_reads_file_path {input.stats} \
            --library_name {params.library} 2> {log}
        """


rule wigToBigWig:
    input:
        wig="{mapping}tracks/{norm}/{method}-{condition}-{replicate}.{norm}.{strand}.wig",
        genomeSize=rules.genomeSize.output
    output:
        bw=report(
            "{mapping}tracks/{norm}/{method}-{condition}-{replicate}.{norm}.{strand}.{mapping}.bw",
            caption="../report/coveragetracks.rst",
            category="Coverage tracks",
            subcategory="{mapping}",
            labels={
                "library": "{method}-{condition}-{replicate}",
                "normalization": "{norm}",
                "strand": "{strand}"
            }
        )
    conda:
        "../envs/wig.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    log:
        "logs/wigToBigWig_{mapping}_{norm}_{method}-{condition}-{replicate}_{strand}.log"
    shell:
        "wigToBigWig {input.wig} {input.genomeSize} {output.bw} 2> {log}"

rule bamcompare:
    input:
        bam=expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
        genomeSize=rules.genomeSize.output,
        bamIndex=expand("maplink/{method}-{condition}-{replicate}.bam.bai", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
    output:
        "figures/results.npz"
    conda:
        "../envs/wig.yaml"
    threads: 5
    resources:
        mem_mb=20000,
        runtime=120
    shell:
        "multiBamSummary bins --smartLabels --bamfiles {input.bam} -o {output} -p {threads};"

rule plotCorrelation:
    input:
        npz="figures/results.npz"
    output:
        correlation=report("figures/heatmap_SpearmanCorr_readCounts.pdf", caption="../report/correlation.rst", category="Quality control")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "plotCorrelation -in {input.npz} --corMethod spearman --skipZeros --plotTitle \"Spearman Correlation of Read Counts\" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.correlation} --outFileCorMatrix SpearmanCorr_readCounts.tab"

rule annotationBed:
    input:
        rules.checkAnnotation.output
    output:
        "tracks/annotation.bed"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "cat {input[0]} | grep -v '\tgene\t' > tracks/annotation-woGenes.gtf; gtf2bed < tracks/annotation-woGenes.gtf > tracks/annotation.bed"

rule annotationBigBed:
    input:
        rules.annotationBed.output,
        rules.genomeSize.output
    output:
        "tracks/annotation.bb"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "cut -f1-6 {input[0]} > tracks/annotationNScore.bed6;  awk '{{$5=1 ; print ;}}' tracks/annotation.bed6 > tracks/annotation.bed6; bedToBigBed -type=bed6 -tab tracks/annotation.bed6 {input[1]} tracks/annotation.bb"

rule colorBigWig:
    input:
        infwd= "tracks/{method}-{condition}-{replicate}.fwd.bw",
        inrev= "tracks/{method}-{condition}-{replicate}.rev.bw"
    output:
        outfwd= "tracks/color/{method}-{condition}-{replicate}.fwd.bedgraph.gz",
        outrev= "tracks/color/{method}-{condition}-{replicate}.rev.bedgraph.gz"
    conda:
        "../envs/color.yaml"
    threads: 1
    params:
        unzippedfwd=lambda wildcards, output: (os.path.splitext(output.outfwd)[0]),
        unzippedrev=lambda wildcards, output: (os.path.splitext(output.outrev)[0])
    shell:
        """
        set +e
        bigWigToWig {input.infwd} {params.unzippedfwd}
        bigWigToWig {input.inrev} {params.unzippedrev}
        sed -i '2s/^/track type=wiggle_0 visibility=full color=0,0,128 autoscale=on\\n/' {params.unzippedfwd}
        sed -i '2s/^/track type=wiggle_0 visibility=full color=0,130,200 autoscale=on\\n/' {params.unzippedrev}
        gzip -f {params.unzippedfwd}
        gzip -f {params.unzippedrev}
        """

rule colorGFF:
    input:
        rbs="tracks/potentialRibosomeBindingSite.gff",
        start="tracks/potentialStartCodons.gff",
        stop="tracks/potentialStopCodons.gff"
    output:
        outrbs="tracks/color/potentialRibosomeBindingSite.gff",
        outstart="tracks/color/potentialStartCodons.gff",
        outstop="tracks/color/potentialStopCodons.gff"
    threads: 1
    shell:
        """
        set +e
        cp {input.rbs} ./tracks/color/
        cp {input.start} ./tracks/color/
        cp {input.stop} ./tracks/color/
        sed -i '1s/^/##track type=wiggle_0 visibility=full color=145,30,180 autoscale=on\\n/' {output.outrbs}
        sed -i '1s/^/##track type=wiggle_0 visibility=full color=210,245,60 autoscale=on\\n/' {output.outstart}
        sed -i '1s/^/##track type=wiggle_0 visibility=full color=230,25,75 autoscale=on\\n/' {output.outstop}
        """
