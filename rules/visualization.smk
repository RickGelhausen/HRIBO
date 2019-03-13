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
        "mkdir -p genomes; cut -f1,2 {input[0]} > genomes/sizes.genome"

rule reversecomplementGenome:
    input:
        rules.retrieveGenome.output
    output:
        "genomes/genome.rev.fa"
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "mkdir -p genomes; SPtools/scripts/reverseComplement.py --input_fasta_filepath genomes/genome.fa --output_fasta_filepath genomes/genome.rev.fa"

rule startCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        "tracks/potentialStartCodons.gff"
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string ATG,GTG,TTG --output_gff3_filepath {output}"

rule stopCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        "tracks/potentialStopCodons.gff"
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string ATG,GTG,TTG --output_gff3_filepath {output}"

rule rbsTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        "tracks/potentialRibosomeBindingSite.gff"
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string TAG,TGA,TAA --output_gff3_filepath {output}"


rule bamindex:
    input:
        rules.maplink.output,
        rules.genomeSize.output
    output:
        "maplink/{method}-{condition}-{replicate}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "samtools index -@ {threads} maplink/{params.prefix}"

rule wigrev:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output
    output:
        rev=report("tracks/{method}-{condition}-{replicate}.rev.bw", caption="../report/wig.rst", category="Mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 5
    params:
        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
    shell:
        "mkdir -p tracks; bamCoverage --normalizeUsing BPM -p {threads} --filterRNAstrand forward -b {input.bam} -o {output.rev};"

rule wigfwd:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output
    output:
        fwd=report("tracks/{method}-{condition}-{replicate}.fwd.bw", caption="../report/wig.rst", category="Mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 5
    params:
        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
    shell:
        "mkdir -p tracks; bamCoverage --normalizeUsing BPM -p {threads} --filterRNAstrand reverse -b {input.bam} -o {output.fwd};"

rule bamcompare:
    input:
        "qc/multi/multiqc_report.html"
    output:
        "figures/results.npz"
    conda:
        "../envs/wig.yaml"
    threads: 5
    shell:
        "mkdir -p tracks; multiBamSummary bins --smartLabels --bamfiles maplink/*.bam -o {output} -p {threads};"

rule plotCorrelation:
    input:
        npz="figures/results.npz"
    output:
        correlation=report("figures/heatmap_SpearmanCorr_readCounts.pdf", caption="../report/correlation.rst", category="Quality control")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "mkdir -p figures; plotCorrelation -in {input.npz} --corMethod spearman --skipZeros --plotTitle \"Spearman Correlation of Read Counts\" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o {output.correlation} --outFileCorMatrix SpearmanCorr_readCounts.tab"

rule annotationBed:
    input:
        rules.retrieveAnnotation.output
    output:
        "tracks/annotation.bed"
    conda:
        "../envs/bed.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; cat {input[0]} | grep -v '\tgene\t' > tracks/annotation-woGenes.gtf; gtf2bed < tracks/annotation-woGenes.gtf > tracks/annotation.bed"

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
        "mkdir -p tracks; cut -f1-6 {input[0]} > tracks/annotationNScore.bed6;  awk '{{$5=1 ; print ;}}' tracks/annotation.bed6 > tracks/annotation.bed6; bedToBigBed -type=bed6 -tab tracks/annotation.bed6 {input[1]} tracks/annotation.bb"
