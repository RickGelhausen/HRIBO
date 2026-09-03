rule genomeSamToolsIndex:
    input:
        genome=rules.retrieveGenome.output
    output:
        "genomes/genome.fa.fai"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    log:
        "logs/genomeSamToolsIndex.log"
    shell:
        "samtools faidx {input.genome:q} 2> {log:q}"

rule genomeSize:
    input:
        index=rules.genomeSamToolsIndex.output
    output:
        "genomes/sizes.genome"
    conda:
        "../envs/samtools.yaml"
    threads: 1
    shell:
        "cut -f1,2 {input.index:q} > {output:q}"

rule reversecomplementGenome:
    input:
        genome=rules.retrieveGenome.output,
        script=str(SCRIPTS / "reverse_complement.py")
    output:
        "genomes/genome.rev.fa"
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} --input_fasta_filepath {input.genome:q} --output_fasta_filepath {output:q}"

rule startCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output,
        script=str(SCRIPTS / "motif_to_gff.py")
    output:
        report("tracks/potentialStartCodons.gff", caption="../report/startCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} --input_genome_fasta_filepath {input.fwd:q} --input_reverse_genome_fasta_filepath {input.rev:q} --motif_string ATG --output_gff3_filepath {output:q}"

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
        rev=rules.reversecomplementGenome.output,
        script=str(SCRIPTS / "motif_to_gff.py")
    output:
        report("tracks/potentialStopCodons.gff", caption="../report/stopCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} --input_genome_fasta_filepath {input.fwd:q} --input_reverse_genome_fasta_filepath {input.rev:q} --motif_string TAG,TGA,TAA --output_gff3_filepath {output:q}"

rule rbsTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output,
        script=str(SCRIPTS / "motif_to_gff.py")
    output:
        report("tracks/potentialRibosomeBindingSite.gff", caption="../report/rbsTrack.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "python3 {input.script:q} --input_genome_fasta_filepath {input.fwd:q} --input_reverse_genome_fasta_filepath {input.rev:q} --motif_string AAGG --output_gff3_filepath {output:q}"


rule bamindex:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output
    output:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=60
    shell:
        "samtools index -@ {threads} {input.bam:q}"

rule totalmappedbamindex:
    input:
        bam=rules.sammultitobam.output,
        genomeSize=rules.genomeSize.output
    output:
        temp("bammulti/{method}-{condition}-{replicate}.bam.bai")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=60
    shell:
        "samtools index -@ {threads} {input.bam:q}"

rule uniquemappedbamindex:
    input:
        bam=rules.samtobam.output,
        genomeSize=rules.genomeSize.output
    output:
        temp("rRNAbam/{method}-{condition}-{replicate}.bam.bai")
    conda:
        "../envs/samtools.yaml"
    threads: 20
    resources:
        mem_mb=40000,
        runtime=60
    shell:
        "samtools index -@ {threads} {input.bam:q}"

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
        genomeSize=rules.genomeSize.output,
        script=str(SCRIPTS / "mapping.py"),
        script_deps=[
            str(SCRIPTS / "lib" / "__init__.py"),
            str(SCRIPTS / "lib" / "library.py"),
        ]
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
        library=library_name,
        wiggle_dir=lambda wildcards: f"{wildcards.mapping}tracks/"
    log:
        "logs/{mapping}tracks_{method}-{condition}-{replicate}.log"
    shell:
        """
        python3 {input.script:q} \
            --mapping_style {params.mapping_style:q} \
            --bam_path {input.bam:q} \
            --wiggle_file_path {params.wiggle_dir:q} \
            --no_of_aligned_reads_file_path {input.stats:q} \
            --library_name {params.library:q} 2> {log:q}
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
        "wigToBigWig {input.wig:q} {input.genomeSize:q} {output.bw:q} 2> {log:q}"

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
        "multiBamSummary bins --smartLabels --bamfiles {input.bam:q} -o {output:q} -p {threads}"

rule plotCorrelation:
    input:
        npz="figures/results.npz"
    output:
        correlation=report("figures/heatmap_SpearmanCorr_readCounts.pdf", caption="../report/correlation.rst", category="Quality control"),
        matrix="figures/SpearmanCorr_readCounts.tab"
    conda:
        "../envs/wig.yaml"
    threads: 1
    log:
        "logs/plotCorrelation.log"
    shell:
        "plotCorrelation -in {input.npz:q} --corMethod spearman --skipZeros --plotTitle \"Spearman Correlation of Read Counts\" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.correlation:q} --outFileCorMatrix {output.matrix:q} > {log:q} 2>&1"
