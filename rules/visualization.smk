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

#rule wigrev:
#    input:
#        bam=rules.maplink.output,
#        genomeSize=rules.genomeSize.output,
#        bamIndex=rules.bamindex.output
#    output:
#        rev=report("tracks/{method}-{condition}-{replicate}.rev.bw", caption="../report/wig.rst", category="Mapped tracks")
#    conda:
#        "../envs/wig.yaml"
#    threads: 5
#    params:
#        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
#    shell:
#        "mkdir -p tracks; bamCoverage --normalizeUsing BPM -p {threads} --scaleFactor=-1 --binSize=1 --smoothLength=0 --filterRNAstrand forward -b {input.bam} -o {output.rev};"

#rule wigfwd:
#    input:
#        bam=rules.maplink.output,
#        genomeSize=rules.genomeSize.output,
#        bamIndex=rules.bamindex.output
#    output:
#        fwd=report("tracks/{method}-{condition}-{replicate}.fwd.bw", caption="../report/wig.rst", category="Mapped tracks")
#    conda:
#        "../envs/wig.yaml"
#    threads: 5
#    params:
#        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
#    shell:
#        "mkdir -p tracks; bamCoverage --normalizeUsing BPM -p {threads} --binSize=1 --smoothLength=0 --filterRNAstrand reverse -b {input.bam} -o {output.fwd};"

rule totalmappedbamindex:
    input:
        rules.sammultitobam.output,
        rules.genomeSize.output
    output:
        "bammulti/{method}-{condition}-{replicate}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "samtools index -@ {threads} bammulti/{params.prefix}"

rule uniquemappedbamindex:
    input:
        rules.samtobam.output,
        rules.genomeSize.output
    output:
        "rRNAbam/{method}-{condition}-{replicate}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads: 20
    params:
        prefix=lambda wildcards, output: (os.path.splitext(os.path.basename(output[0]))[0])
    shell:
        "samtools index -@ {threads} rRNAbam/{params.prefix}"

rule totalmappedreadcountstats:
    input:
        bam=rules.totalmappedsamtobam.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.totalmappedbamindex.output
    output:
        stat="bammulti/{method}-{condition}-{replicate}.readstats"
    threads: 1
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p rRNAbam; readstats.py --bam_path {input.bam} > {output.stat}; source deactivate;"

rule uniquemappedreadcountstats:
    input:
        bam=rules.samtobam.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.uniquemappedbamindex.output
    output:
        stat="rRNAbam/{method}-{condition}-{replicate}.readstats"
    threads: 1
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p rRNAbam; readstats.py --bam_path {input.bam} > {output.stat}; source deactivate;"

rule totalmappedminreadcounts:
    input:
        stats=expand("bammulti/{method}-{condition}-{replicate}.readstats", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        minreads="bammulti/minreads.txt"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p rRNAbam; minreads.py {input.stats} > {output.minreads}; source deactivate;"

rule uniquemappedminreadcounts:
    input:
        stats=expand("rRNAbam/{method}-{condition}-{replicate}.readstats", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"])
    output:
        minreads="rRNAbam/minreads.txt"
    threads: 1
    params:
        prefix=lambda wildcards, output: (os.path.splitext(output[0])[0])
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p rRNAbam; minreads.py {input.stats} > {output.minreads}; source deactivate;"

rule totalmappedwig:
    input:
        bam=rules.sammultitobam.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.totalmappedbamindex.output,
        stats="bammulti/{method}-{condition}-{replicate}.readstats",
        min="bammulti/minreads.txt"
    output:
        fwd="totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="totalmappedtracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="totalmappedtracks/min/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p totalmappedtracks; mkdir -p totalmappedtracks/raw; mkdir -p totalmappedtracks/mil; mkdir -p totalmappedtracks/min; coverage.py --coverage_style global --bam_path {input.bam} --wiggle_file_path totalmappedtracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule uniquemappedwig:
    input:
        bam=rules.samtobam.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.uniquemappedbamindex.output,
        stats="rRNAbam/{method}-{condition}-{replicate}.readstats",
        min="rRNAbam/minreads.txt"
    output:
        fwd="uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="uniquemappedtracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="uniquemappedtracks/min/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p uniquemappedtracks; mkdir -p uniquemappedtracks/raw; mkdir -p uniquemappedtracks/mil; mkdir -p uniquemappedtracks/min; coverage.py --coverage_style global --bam_path {input.bam} --wiggle_file_path uniquemappedtracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule totalmappedwigtobigwigrawforward:
    input:
        fwd="totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.totalmapped.bw", caption="../report/totalmappedwig.rst", category="Total mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule totalmappedwigtobigwigminrawreverse:
    input:
        rev="totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.totalmapped.bw", caption="../report/totalmappedwig.rst", category="Total mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule totalmappedwigtobigwigminforward:
    input:
        fwd="totalmappedtracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("totalmappedtracks/min/{method}-{condition}-{replicate}.min.forward.totalmapped.bw", caption="../report/totalmappedwig.rst", category="Total mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule totalmappedwigtobigwigminreverse:
    input:
        rev="totalmappedtracks/min/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("totalmappedtracks/min/{method}-{condition}-{replicate}.min.reverse.totalmapped.bw", caption="../report/totalmappedwig.rst", category="Total mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule totalmappedwigtobigwigmilforward:
    input:
        fwd="totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.totalmapped.bw", caption="../report/totalmappedwig.rst", category="Total mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule totalmappedwigtobigwigmilreverse:
    input:
        rev="totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.totalmapped.bw", caption="../report/totalmappedwig.rst", category="Total mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule uniquemappedwigtobigwigrawforward:
    input:
        fwd="uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.uniquemapped.bw", caption="../report/uniquemappedwig.rst", category="Unique mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule uniquemappedwigtobigwigminrawreverse:
    input:
        rev="uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.uniquemapped.bw", caption="../report/uniquemappedwig.rst", category="Unique mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule uniquemappedwigtobigwigminforward:
    input:
        fwd="uniquemappedtracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("uniquemappedtracks/min/{method}-{condition}-{replicate}.min.forward.uniquemapped.bw", caption="../report/uniquemappedwig.rst", category="Unique mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule uniquemappedwigtobigwigminreverse:
    input:
        rev="uniquemappedtracks/min/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("uniquemappedtracks/min/{method}-{condition}-{replicate}.min.reverse.uniquemapped.bw", caption="../report/uniquemappedwig.rst", category="Unique mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule uniquemappedwigtobigwigmilforward:
    input:
        fwd="uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.uniquemapped.bw", caption="../report/uniquemappedwig.rst", category="Unique mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule uniquemappedwigtobigwigmilreverse:
    input:
        rev="uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.uniquemapped.bw", caption="../report/uniquemappedwig.rst", category="Unique mapped tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

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

rule globalwig:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        stats="maplink/{method}-{condition}-{replicate}.readstats",
        min="maplink/minreads.txt"
    output:
        fwd="globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="globaltracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="globaltracks/min/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p globaltracks; mkdir -p globaltracks/raw; mkdir -p globaltracks/mil; mkdir -p globaltracks/min; coverage.py --coverage_style global --bam_path {input.bam} --wiggle_file_path globaltracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule globalwigtobigwigrawforward:
    input:
        fwd="globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.global.bw", caption="../report/globalwig.rst", category="Global tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule globalwigtobigwigminrawreverse:
    input:
        rev="globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.global.bw", caption="../report/globalwig.rst", category="Global tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule globalwigtobigwigminforward:
    input:
        fwd="globaltracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("globaltracks/min/{method}-{condition}-{replicate}.min.forward.global.bw", caption="../report/globalwig.rst", category="Global tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule globalwigtobigwigminreverse:
    input:
        rev="globaltracks/min/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("globaltracks/min/{method}-{condition}-{replicate}.min.reverse.global.bw", caption="../report/globalwig.rst", category="Global tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule globalwigtobigwigmilforward:
    input:
        fwd="globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.global.bw", caption="../report/globalwig.rst", category="Global tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule globalwigtobigwigmilreverse:
    input:
        rev="globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.global.bw", caption="../report/globalwig.rst", category="Global tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"


rule centeredwig:
    input:
        bam=rules.maplink.output,
        genomeSize=rules.genomeSize.output,
        bamIndex=rules.bamindex.output,
        stats="maplink/{method}-{condition}-{replicate}.readstats",
        min="maplink/minreads.txt"
    output:
        fwd="centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="centeredtracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p centeredtracks; mkdir -p centeredtracks/raw; mkdir -p centeredtracks/mil; mkdir -p centeredtracks/min; coverage.py --coverage_style centered --bam_path {input.bam} --wiggle_file_path centeredtracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule centeredwigtobigwigrawforward:
    input:
        fwd="centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.centered.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule centeredwigtobigwigminrawreverse:
    input:
        rev="centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.centered.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule centeredwigtobigwigminforward:
    input:
        fwd="centeredtracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("centeredtracks/min/{method}-{condition}-{replicate}.min.forward.centered.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule centeredwigtobigwigminreverse:
    input:
        rev="centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.centered.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule centeredwigtobigwigmilforward:
    input:
        fwd="centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.centered.bw", caption="../report/centeredwig.rst", category="Centered tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule centeredwigtobigwigmilreverse:
    input:
        rev="centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.centered.bw", caption="../report/centeredwig.rst", category="Centered tracks")
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
        fwd="fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p fiveprimetracks; mkdir -p fiveprimetracks/raw; mkdir -p fiveprimetracks/mil; mkdir -p fiveprimetracks/min; coverage.py --coverage_style first_base_only --bam_path {input.bam} --wiggle_file_path fiveprimetracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule fiveprimewigtobigwigrawforward:
    input:
        fwd="fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.fiveprime.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule fiveprimewigtobigwigrawreverse:
    input:
        rev="fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.fiveprime.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule fiveprimewigtobigwigminforward:
    input:
        fwd="fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.fiveprime.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule fiveprimewigtobigwigminreverse:
    input:
        rev="fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.fiveprime.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule fiveprimewigtobigwigmilforward:
    input:
        fwd="fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.fiveprime.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule fiveprimewigtobigwimilgreverse:
    input:
        rev="fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.fiveprime.bw", caption="../report/fiveprimewig.rst", category="5' single nucleotide mapping tracks")
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
        fwd="threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        rev="threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        fmil="threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        rmil="threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        fmin="threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        rmin="threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.wig"
    threads: 1
    params:
        prefix=lambda wildcards, output: (Path(output[0]).stem).strip('.raw.forward.wig'),
        prefixpath=lambda wildcards, output: (os.path.dirname(output.fwd))
    shell:
        "source activate /scratch/bi03/egg/miniconda3/envs/coverage; mkdir -p threeprimetracks; mkdir -p threeprimetracks/raw; mkdir -p threeprimetracks/mil; mkdir -p threeprimetracks/min; coverage.py --coverage_style last_base_only --bam_path {input.bam} --wiggle_file_path threeprimetracks/ --no_of_aligned_reads_file_path {input.stats} --library_name {params.prefix} --min_no_of_aligned_reads_file_path {input.min}; source deactivate;"

rule threeprimewigtobigwigrawforward:
    input:
        fwd="threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.threeprime.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule threeprimewigtobigwigrawreverse:
    input:
        rev="threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.threeprime.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule threeprimewigtobigwigminforward:
    input:
        fwd="threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.threeprime.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule threeprimewigtobigwigminreverse:
    input:
        rev="threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.threeprime.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.rev} {input.genomeSize} {output.rev}"

rule threeprimewigtobigwigmilforward:
    input:
        fwd="threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.wig",
        genomeSize=rules.genomeSize.output
    output:
        fwd=report("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.threeprime.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mapping tracks")
    conda:
        "../envs/wig.yaml"
    threads: 1
    shell:
        "wigToBigWig {input.fwd} {input.genomeSize} {output.fwd}"

rule threeprimewigtobigwigmilreverse:
    input:
        rev="threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.wig",
        genomeSize=rules.genomeSize.output
    output:
        rev=report("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.threeprime.bw", caption="../report/threeprimewig.rst", category="3' single nucleotide mapping tracks")
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
