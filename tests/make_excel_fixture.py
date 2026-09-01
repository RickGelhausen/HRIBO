"""Build deterministic inputs for the excel-generating scripts.

These scripts have no tests and produce spreadsheets that people read directly,
so any refactor has to be checked against the exact output of the version before
it. This module builds the inputs; tests/test_excel_outputs.py captures the
outputs and compares them.
"""

import random
from pathlib import Path

CONTIGS = {"NC_000913.3": 40000, "pPlasmid1": 14000}

# Deliberately a mix: RIBO with matching RNA (so translational efficiency
# columns appear), two conditions, two replicates.
WILDCARDS = [
    "RIBO-A-1", "RIBO-A-2", "RIBO-B-1", "RIBO-B-2",
    "RNA-A-1", "RNA-A-2", "RNA-B-1", "RNA-B-2",
]

FEATURES = ["CDS", "gene", "rRNA", "tRNA", "sRNA", "pseudogene", "transcript", "five_prime_UTR"]


def write_genome(path):
    rng = random.Random(11)
    with open(path, "w") as handle:
        for name, length in CONTIGS.items():
            handle.write(f">{name} simulated sequence\n")
            sequence = "".join(rng.choice("ACGT") for _ in range(length))
            for start in range(0, length, 70):
                handle.write(sequence[start : start + 70] + "\n")


def _features():
    """(contig, feature, start, stop, strand, index), 1-based inclusive."""
    entries = []
    index = 0
    for name, length in CONTIGS.items():
        position = 100
        while position + 350 < length:
            index += 1
            feature = FEATURES[index % len(FEATURES)]
            strand = "+" if index % 2 else "-"
            # Lengths divisible by three, so amino acid translation is exercised.
            entries.append((name, feature, position, position + 299, strand, index))
            position += 400
    return entries


def write_total_mapped(path):
    """<wildcard>\\t<contig>\\t<count>, as scripts/total_mapped_reads.py writes."""
    rng = random.Random(5)
    with open(path, "w") as handle:
        for wildcard in WILDCARDS:
            for contig in CONTIGS:
                handle.write(f"{wildcard}\t{contig}\t{rng.randint(200_000, 900_000)}\n")


def _attributes(index, feature, with_prediction=None):
    parts = [
        f"ID=feature{index}",
        f"locus_tag=b{index:04d}",
        f"old_locus_tag=old{index:04d}",
        f"Name=gene{index}",
        f"product=hypothetical protein {index}",
        f"Note=an annotation note {index}",
    ]
    if with_prediction == "reparation":
        parts.append(f"prob={0.5 + (index % 50) / 100:.2f}")
        parts.append("evidence=reparation")
    elif with_prediction == "deepribo":
        parts.append(f"pred_value={1.0 + (index % 40) / 10:.2f}")
        parts.append("evidence=deepribo")
        parts.append(f"novel_rank={index % 5}")
    return ";".join(parts) + ";"


def write_read_counts(path, source, with_prediction=None, features=None):
    """A GFF whose extra columns hold one read count per library.

    This is the shape map_reads_to_annotation.py produces and that the excel
    scripts consume: nine GFF columns followed by one column per wildcard.
    """
    rng = random.Random(7)
    with open(path, "w") as handle:
        for contig, feature, start, stop, strand, index in _features():
            if features is not None and feature not in features:
                continue
            score = f"{index}" if with_prediction == "deepribo" else "."
            phase = "0"
            counts = [str(rng.randint(0, 5000)) for _ in WILDCARDS]
            row = [
                contig, source, feature, str(start), str(stop), score, strand, phase,
                _attributes(index, feature, with_prediction),
            ] + counts
            handle.write("\t".join(row) + "\n")


def write_annotation(path):
    """A plain GFF3 without read count columns, for the diffex scripts."""
    with open(path, "w") as handle:
        handle.write("##gff-version 3\n")
        for contig, feature, start, stop, strand, index in _features():
            if feature in ("gene", "pseudogene"):
                attributes = (
                    f"ID=gene{index};locus_tag=b{index:04d};"
                    f"old_locus_tag=old{index:04d};gene=name{index};"
                )
            else:
                attributes = (
                    f"ID=feature{index};Parent=gene{index};"
                    f"locus_tag=b{index:04d};Name=gene{index};"
                )
            handle.write("\t".join([
                contig, "sim", feature, str(start), str(stop), ".", strand, "0", attributes,
            ]) + "\n")


def _diffex_ids():
    return [
        (f"{contig}:{start}-{stop}:{strand}", index)
        for contig, feature, start, stop, strand, index in _features()
        if feature == "CDS"
    ]


def write_riborex(path):
    """As riborex.R writes it: write.csv with row names, so column one is unnamed."""
    rng = random.Random(21)
    with open(path, "w") as handle:
        handle.write('"",baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n')
        for identifier, index in _diffex_ids():
            handle.write(
                f"{identifier},{rng.uniform(1, 900):.4f},{rng.uniform(-4, 4):.4f},"
                f"{rng.uniform(0, 1):.4f},{rng.uniform(-3, 3):.4f},"
                f"{rng.uniform(0, 1):.5f},{rng.uniform(0, 1):.5f}\n"
            )


def write_xtail(path):
    """As xtail.R writes it. The last column is R-named 'pvalue.adjust'."""
    rng = random.Random(22)
    columns = [
        "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1",
        "log2FC_TE_v2", "pvalue_v2", "log2FC_TE_final", "pvalue_final",
        "pvalue.adjust",
    ]
    with open(path, "w") as handle:
        handle.write('"",' + ",".join(columns) + "\n")
        for identifier, index in _diffex_ids():
            values = [f"{rng.uniform(-4, 4):.4f}" for _ in range(8)]
            handle.write(f"{identifier}," + ",".join(values) + f",{rng.uniform(0, 1):.5f}\n")


def write_deltate(path, prefix, seed=23):
    """deltaTE fold-change tables, as DESeq2's write.table produces them.

    Tab separated with row names, and a header one field shorter than the data
    rows, which is what makes pandas take the first column as the index. The
    consuming script relies on exactly that inference.
    """
    rng = random.Random(seed)
    # Only the TE table carries a Wald statistic; the RIBO and RNA tables do not.
    # The consuming script's own empty-input fallback documents this asymmetry.
    columns = ["baseMean", "log2FoldChange", "lfcSE"]
    if prefix == "TE":
        columns.append("stat")
    columns += ["pvalue", "padj"]

    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for identifier, index in _diffex_ids():
            values = [f"{rng.uniform(-4, 900):.4f}" for _ in columns]
            handle.write(identifier + "\t" + "\t".join(values) + "\n")



def write_pooled_diffex(path, tool, contrasts=("B-A",)):
    """The pooled <tool>_all.csv the overview table reads.

    Different schema from the per-contrast files: keyed by an explicit gene_id
    column with a contrast column written as "contrast_<name>".
    """
    columns = {
        "riborex": ["log2FoldChange", "pvalue", "padj"],
        "xtail": ["log2FC_TE_final", "pvalue_final", "pvalue_adjust"],
        "deltate": [
            "RIBO_log2FoldChange", "RIBO_pvalue", "RIBO_padj",
            "RNA_log2FoldChange", "RNA_pvalue", "RNA_padj",
            "TE_log2FoldChange", "TE_pvalue", "TE_padj",
        ],
    }[tool]

    rng = random.Random(31 + len(tool))
    with open(path, "w") as handle:
        handle.write("gene_id," + ",".join(columns) + ",contrast\n")
        for contrast in contrasts:
            for identifier, index in _diffex_ids():
                values = [f"{rng.uniform(-4, 4):.4f}" for _ in columns]
                handle.write(f"{identifier}," + ",".join(values) + f",contrast_{contrast}\n")


def build(directory):
    """Write every input the excel scripts need. Returns the directory."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    write_genome(directory / "genome.fa")
    write_annotation(directory / "annotation.gff")
    write_total_mapped(directory / "total_mapped_reads.txt")

    write_read_counts(directory / "total_annotation.gtf", "HRIBO")
    write_read_counts(
        directory / "reparation_annotation.gff", "reparation",
        with_prediction="reparation", features=["CDS"],
    )
    write_read_counts(
        directory / "deepribo_annotation.gff", "deepribo",
        with_prediction="deepribo", features=["CDS"],
    )

    write_riborex(directory / "riborex_all.csv")
    write_xtail(directory / "xtail_all.csv")
    write_pooled_diffex(directory / "riborex_pooled.csv", "riborex")
    write_pooled_diffex(directory / "xtail_pooled.csv", "xtail")
    write_pooled_diffex(directory / "deltate_pooled.csv", "deltate")

    write_deltate(directory / "deltaRibo.txt", "RIBO", seed=23)
    write_deltate(directory / "deltaRNA.txt", "RNA", seed=24)
    write_deltate(directory / "deltaTE.txt", "TE", seed=25)

    return directory


if __name__ == "__main__":
    import sys

    print(build(sys.argv[1]))
