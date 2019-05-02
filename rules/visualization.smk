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
        report("tracks/potentialStartCodons.gff", caption="../report/startCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string ATG --output_gff3_filepath {output}"

rule alternativeStartCodonTrack:
    input:
        fwd=rules.retrieveGenome.output,
        rev=rules.reversecomplementGenome.output
    output:
        report("tracks/potentialAlternativeStartCodons.gff", caption="../report/startCodons.rst", category="Annotation")
    conda:
        "../envs/biopython.yaml"
    threads: 1
    shell:
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string GTG,TTG,CTG --output_gff3_filepath {output}"


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
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string TAG,TGA,TAA --output_gff3_filepath {output}"

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
        "mkdir -p tracks; SPtools/scripts/motif2GFF3.py --input_genome_fasta_filepath {input.fwd} --input_reverse_genome_fasta_filepath {input.rev} --motif_string AAGG --output_gff3_filepath {output}"


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
        "mkdir -p tracks; bamCoverage --normalizeUsing BPM -p {threads} --scaleFactor=-1 --binSize=1 --filterRNAstrand forward -b {input.bam} -o {output.rev};"

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
        "mkdir -p tracks; bamCoverage --normalizeUsing BPM -p {threads} --binSize=1 --filterRNAstrand reverse -b {input.bam} -o {output.fwd};"

rule readcountstats:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output
    output:
        stat="maplink/{method}-{condition}-{replicate}.readstats"
    threads: 1
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p tracks; readstats.py --bam_path {input.bam} > {output.stat}; source deactivate;"

rule minreadcounts:
    input:
        stats=expand("maplink/{method}-{condition}-{replicate}.readstats", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        minreads="maplink/minreads.txt"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p tracks; minreads.py {input.stats} > {output.minreads}; source deactivate;"

rule centeredwig:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        stats="maplink/{method}-{condition}-{replicate}.readstats",
        min="maplink/minreads.txt"
    output:
        fwd="centeredtracks/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="centeredtracks/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="centeredtracks/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="centeredtracks/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="centeredtracks/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="centeredtracks/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p centeredtracks; coverage.py --coverage_style centered --bam_path {input.bam} --wiggle_file_path centeredtracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule wigtobigwigrawforward:
    input:
        fwd="centeredtracks/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("centeredtracks/{method}-{condition}-{replicate}.raw.forward.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule wigtobigwigminrawreverse:
    input:
        rev="centeredtracks/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/{method}-{condition}-{replicate}.raw.reverse.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule wigtobigwigminforward:
    input:
        fwd="centeredtracks/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("centeredtracks/{method}-{condition}-{replicate}.min.forward.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule wigtobigwigminreverse:
    input:
        rev="centeredtracks/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/{method}-{condition}-{replicate}.min.reverse.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule wigtobigwigmilforward:
    input:
        fwd="centeredtracks/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("centeredtracks/{method}-{condition}-{replicate}.mil.forward.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule wigtobigwigmilreverse:
    input:
        rev="centeredtracks/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/{method}-{condition}-{replicate}.mil.reverse.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"
	
rule fiveprimewig:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        stats="maplink/{method}-{condition}-{replicate}.readstats",
        min="maplink/minreads.txt"
    output:
        fwd="fiveprimetracks/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="fiveprimetracks/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="fiveprimetracks/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="fiveprimetracks/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="fiveprimetracks/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="fiveprimetracks/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p fiveprimetracks; coverage.py --coverage_style first_base_only --bam_path {input.bam} --wiggle_file_path fiveprimetracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule fiveprimewigtobigwigrawforward:
    input:
        fwd="fiveprimetracks/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("fiveprimetracks/{method}-{condition}-{replicate}.raw.forward.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule fiveprimewigtobigwigrawreverse:
    input:
        rev="centeredtracks/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("fiveprimetracks/{method}-{condition}-{replicate}.raw.reverse.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule fiveprimewigtobigwigminforward:
    input:
        fwd="centeredtracks/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("fiveprimetracks/{method}-{condition}-{replicate}.min.forward.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule fiveprimewigtobigwigminreverse:
    input:
        rev="fiveprimetracks/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("fiveprimetracks/{method}-{condition}-{replicate}.min.reverse.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule fiveprimewigtobigwigmilforward:
    input:
        fwd="fiveprimetracks/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("fiveprimetracks/{method}-{condition}-{replicate}.mil.forward.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule fiveprimewigtobigwimilgreverse:
    input:
        rev="fiveprimetracks/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/{method}-{condition}-{replicate}.mil.reverse.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule threeprimewig:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        stats="maplink/{method}-{condition}-{replicate}.readstats",
        min="maplink/minreads.txt"
    output:
        fwd="threeprimetracks/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="threeprimetracks/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="threeprimetracks/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="threeprimetracks/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="threeprimetracks/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="threeprimetracks/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p threeprimetracks; coverage.py --coverage_style last_base_only --bam_path {input.bam} --wiggle_file_path threeprimetracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule threeprimewigtobigwigrawforward:
    input:
        fwd="centeredtracks/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("threeprime/{method}-{condition}-{replicate}.raw.forward.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule threeprimewigtobigwigrawreverse:
    input:
        rev="centeredtracks/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("threeprimetracks/{method}-{condition}-{replicate}.raw.reverse.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule threeprimewigtobigwigminforward:
    input:
        fwd="centeredtracks/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("threeprimetracks/{method}-{condition}-{replicate}.min.forward.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule threeprimewigtobigwigminreverse:
    input:
        rev="threetracks/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("threeprimetracks/{method}-{condition}-{replicate}.min.reverse.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule threeprimewigtobigwigmilforward:
    input:
        fwd="threetracks/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("threeprimetracks/{method}-{condition}-{replicate}.mil.forward.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule threeprimewigtobigwigmilreverse:
    input:
        rev="threetracks/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("threeprimetracks/{method}-{condition}-{replicate}.mil.reverse.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mappping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"


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
        mkdir -p tracks/color
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
        mkdir -p tracks/color
        cp {input.rbs} ./tracks/color/
        cp {input.start} ./tracks/color/
        cp {input.stop} ./tracks/color/
        sed -i '1s/^/##track type=wiggle_0 visibility=full color=145,30,180 autoscale=on\\n/' {output.outrbs}
        sed -i '1s/^/##track type=wiggle_0 visibility=full color=210,245,60 autoscale=on\\n/' {output.outstart}
        sed -i '1s/^/##track type=wiggle_0 visibility=full color=230,25,75 autoscale=on\\n/' {output.outstop}
        """
