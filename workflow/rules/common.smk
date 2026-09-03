"""
Helper functions and lookup tables shared by the rule files.

Kept separate from the rules themselves, as recommended by the Snakemake best
practices, so that the rule files stay declarative.
"""

from dataclasses import dataclass

from lib.cli import command_option_values


# --------------------------------------------------------------------------
# Coverage tracks
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class TrackSource:
    """How one flavour of coverage track is produced.

    bam_dir      directory holding the alignments the track is computed from
    mapping_style value understood by scripts/mapping.py
    stats        per-library mapped read counts used for the normalisations
    description  shown as the report subcategory
    """

    bam_dir: str
    mapping_style: str
    stats: str
    description: str


# The directory prefix of a track ("globaltracks", "centeredtracks", ...) is the
# {mapping} wildcard, so this table is also the wildcard's domain.
TRACK_SOURCES = {
    "global": TrackSource(
        "maplink", "global", "readcounts/bam_mapped_reads.txt", "Global tracks"
    ),
    "centered": TrackSource(
        "maplink", "centered", "readcounts/bam_mapped_reads.txt", "Centered tracks"
    ),
    "fiveprime": TrackSource(
        "maplink", "first_base_only", "readcounts/bam_mapped_reads.txt", "5' single nucleotide tracks"
    ),
    "threeprime": TrackSource(
        "maplink", "last_base_only", "readcounts/bam_mapped_reads.txt", "3' single nucleotide tracks"
    ),
    # Computed from all alignments including multi-mappers, and from the
    # uniquely mapping alignments before rRNA/tRNA removal. Not requested by the
    # default targets; ask for the files explicitly to build them.
    "totalmapped": TrackSource(
        "bammulti", "global", "readcounts/total_mapped_reads.txt", "Total mapped tracks"
    ),
    "uniquemapped": TrackSource(
        "rRNAbam", "global", "readcounts/unique_mapped_reads.txt", "Uniquely mapped tracks"
    ),
}

# Normalisations written by scripts/mapping.py in a single pass.
#   raw  unnormalised counts
#   mil  counts per million mapped reads
#   min  counts scaled to the smallest library
TRACK_NORMALIZATIONS = ["raw", "mil", "min"]
TRACK_STRANDS = ["forward", "reverse"]

# Track flavours built by the 'tracks' stage.
DEFAULT_TRACK_MAPPINGS = ["global", "centered", "fiveprime", "threeprime"]


def library_name(wildcards):
    """The <method>-<condition>-<replicate> stem shared by every per-library file."""
    return f"{wildcards.method}-{wildcards.condition}-{wildcards.replicate}"


def _library_files(directory, suffix):
    """One path per library, as `<directory>/<library name><suffix>`."""
    return [
        f"{directory}/{method}-{condition}-{replicate}{suffix}"
        for method, condition, replicate in zip(
            samples["method"], samples["condition"], samples["replicate"]
        )
    ]


def get_mapping_files():
    """Targets of the 'mapping' stage: the final alignments and their indices.

    maplink/ holds the uniquely mapping reads left after rRNA and tRNA removal,
    which is what every downstream analysis is computed from. The intermediate
    bammulti/ and rRNAbam/ alignments are not requested here; ask for those
    files explicitly if you need them.
    """
    return _library_files("maplink", ".bam") + _library_files("maplink", ".bam.bai")


def track_bam(wildcards):
    return f"{TRACK_SOURCES[wildcards.mapping].bam_dir}/{library_name(wildcards)}.bam"


def track_bam_index(wildcards):
    return f"{TRACK_SOURCES[wildcards.mapping].bam_dir}/{library_name(wildcards)}.bam.bai"


def track_stats(wildcards):
    return TRACK_SOURCES[wildcards.mapping].stats


def track_mapping_style(wildcards):
    return TRACK_SOURCES[wildcards.mapping].mapping_style


def get_wigfiles():
    """Targets of the 'tracks' stage: one bigwig per library, flavour and strand."""
    return [
        f"{mapping}tracks/{norm}/{method}-{condition}-{replicate}.{norm}.{strand}.{mapping}.bw"
        for mapping in DEFAULT_TRACK_MAPPINGS
        for norm in TRACK_NORMALIZATIONS
        for strand in TRACK_STRANDS
        for method, condition, replicate in zip(
            samples["method"], samples["condition"], samples["replicate"]
        )
    ]


# --------------------------------------------------------------------------
# Quality control
# --------------------------------------------------------------------------


def is_paired_end(row):
    return isinstance(row["fastqFile2"], str) and row["fastqFile2"].strip() != ""


def get_raw_qc_files():
    """FastQC reports of the raw reads, which differ between single and paired end."""
    qc_files = []
    for _, row in samples.iterrows():
        stem = "{method}-{condition}-{replicate}".format(**row)
        if is_paired_end(row):
            qc_files.append(f"qc/1raw/{stem}-raw-q_fastqc.html")
            qc_files.append(f"qc/1raw/{stem}-raw-p_fastqc.html")
        else:
            qc_files.append(f"qc/1raw/{stem}-raw_fastqc.html")
    return qc_files


def get_trimmed_qc_files():
    """FastQC reports of the trimmed reads, before paired reads are merged."""
    qc_files = []
    for _, row in samples.iterrows():
        stem = "{method}-{condition}-{replicate}".format(**row)
        if is_paired_end(row):
            qc_files.append(f"qc/2trimmed/{stem}-trimmed_q_fastqc.html")
            qc_files.append(f"qc/2trimmed/{stem}-trimmed_p_fastqc.html")
        else:
            qc_files.append(f"qc/2trimmed/{stem}-trimmed_fastqc.html")
    return qc_files


def get_processed_read_files():
    """Mapping-ready reads produced by trimming for every library.

    Single-end libraries use Cutadapt's output directly. Paired-end libraries
    use the PEAR assembly consumed by the mapping rules. Both layouts therefore
    publish the same stable path contract.
    """
    return _library_files("trimmed", ".fastq")


def get_qc_files():
    """Inputs of the MultiQC summary that depend on the library layout."""
    qc_files = get_raw_qc_files()
    for _, row in samples.iterrows():
        stem = "{method}-{condition}-{replicate}".format(**row)
        qc_files.append(f"qc/2trimmed/{stem}-trimmed_fastqc.html")
        qc_files.append(f"trimmed/{stem}.fastq")
    return qc_files


def get_trimming_files():
    """Targets of 'trimming': processed reads and their FastQC reports."""
    return get_processed_read_files() + get_raw_qc_files() + get_trimmed_qc_files()


# --------------------------------------------------------------------------
# Read counting
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ReadCountSet:
    """One featureCounts run.

    bam_dir     alignments to count
    annotation  annotation to count against
    extra       additional flags for scripts/call_featurecounts.py
    """

    bam_dir: str
    annotation: str
    extra: tuple[str, ...] = ()


# Keyed by output basename under readcounts/, which is the {countset} wildcard.
READ_COUNT_SETS = {
    "differential_expression_read_counts.csv": ReadCountSet(
        "maplink", "auxiliary/unambigous_annotation.gff", ("--for_diff_expr",)
    ),
    "annotation_independant_read_counts.raw": ReadCountSet(
        "maplink", "auxiliary/unambigous_annotation.gff"
    ),
    "reparation_read_counts.raw": ReadCountSet(
        "maplink", "tracks/reparation_annotated.gff"
    ),
    "deepribo_read_counts.raw": ReadCountSet(
        "maplink", "tracks/deepribo_merged.gff"
    ),
    # Multi-mappers are counted fractionally against the enriched annotation.
    "annotation_total_reads.raw": ReadCountSet(
        "bammulti", "auxiliary/unambigous_annotation.gff", ("--with_M", "--fraction")
    ),
    "annotation_unique_reads.raw": ReadCountSet(
        "rRNAbam", "auxiliary/unambigous_annotation.gff", ("--fraction",)
    ),
}


def readcount_bams(wildcards):
    return _library_files(READ_COUNT_SETS[wildcards.countset].bam_dir, ".bam")


def readcount_bam_indices(wildcards):
    return _library_files(READ_COUNT_SETS[wildcards.countset].bam_dir, ".bam.bai")


def readcount_annotation(wildcards):
    return READ_COUNT_SETS[wildcards.countset].annotation


def readcount_flags(wildcards):
    """Flags beyond the ones every run shares."""
    extra = list(READ_COUNT_SETS[wildcards.countset].extra)
    if wildcards.countset == "differential_expression_read_counts.csv":
        features = config["differentialExpressionSettings"]["features"]
        if features:
            extra.extend(["--use_features", *features])
    return extra


# Annotations that read counts are mapped back onto, keyed by output basename.
READ_COUNT_ANNOTATIONS = {
    "independant_annotation.gff": (
        "readcounts/annotation_independant_read_counts.raw",
        "annotation/annotation_processed.gff",
    ),
    "reparation_annotation.gff": (
        "readcounts/reparation_read_counts.raw",
        "tracks/reparation_annotated.gff",
    ),
    "deepribo_annotation.gff": (
        "readcounts/deepribo_read_counts.raw",
        "tracks/deepribo_merged.gff",
    ),
    "total_annotation.gtf": (
        "readcounts/annotation_total_reads.raw",
        "auxiliary/enriched_annotation.gff",
    ),
    "unique_annotation.gtf": (
        "readcounts/annotation_unique_reads.raw",
        "auxiliary/enriched_annotation.gff",
    ),
}


def mapped_counts_reads(wildcards):
    return READ_COUNT_ANNOTATIONS[wildcards.mapped][0]


def mapped_counts_annotation(wildcards):
    return READ_COUNT_ANNOTATIONS[wildcards.mapped][1]


# Alignment directories summarised by scripts/total_mapped_reads.py, keyed by the
# {source} wildcard of the summary filenames.
MAPPED_READ_SOURCES = {
    "total": "bammulti",
    "unique": "rRNAbam",
    "bam": "maplink",
}


def mapped_read_bams(wildcards):
    return _library_files(MAPPED_READ_SOURCES[wildcards.source], ".bam")


def mapped_read_bam_indices(wildcards):
    return _library_files(MAPPED_READ_SOURCES[wildcards.source], ".bam.bai")


# --------------------------------------------------------------------------
# Overview table
# --------------------------------------------------------------------------


def overview_sources():
    """Inputs of the overview table, which depend on the enabled analyses.

    Returns a mapping of input name to path. The four previous
    createOverviewTable* rules differed only in this set.
    """
    sources = {
        "annotation": "readcounts/independant_annotation.gff",
        "genome": "genomes/genome.fa",
        "totalreads": "readcounts/bam_mapped_reads.txt",
        "reparation": "readcounts/reparation_annotation.gff",
    }
    if DEEPRIBO:
        sources["deepribo"] = "readcounts/deepribo_annotation.gff"
    if DIFFEXPRESS:
        sources["xtail"] = "xtail/xtail_all.csv"
        sources["riborex"] = "riborex/riborex_all.csv"
        sources["deltate"] = "deltate/deltate_all.csv"
    return sources


def overview_flags():
    """Command-line argv matching ``overview_sources`` without shell fragments."""
    sources = overview_sources()
    flags = ["--mapped_reads_reparation", sources["reparation"]]
    if "deepribo" in sources:
        flags.extend(["--mapped_reads_deepribo", sources["deepribo"]])
    for tool in ("xtail", "riborex", "deltate"):
        if tool in sources:
            flags.extend([f"--{tool}", sources[tool]])
    if CONTRASTS:
        flags.extend(["-c", *CONTRASTS])
    return flags
