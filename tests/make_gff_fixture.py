"""Build deterministic inputs for the annotation-transforming scripts.

These scripts rewrite GFF files in the middle of the workflow and had no tests,
so their outputs are snapshotted the same way the excel workbooks are.
"""

from pathlib import Path

CONTIGS = {"NC_000913.3": 20000, "pPlasmid1": 8000}


def _features():
    """(contig, feature, start, stop, strand, index), 1-based inclusive."""
    entries = []
    index = 0
    for name, length in CONTIGS.items():
        position = 200
        while position + 400 < length:
            index += 1
            feature = "CDS" if index % 3 else ("rRNA" if index % 2 else "tRNA")
            strand = "+" if index % 2 else "-"
            entries.append((name, feature, position, position + 299, strand, index))
            position += 450
    return entries


def write_gff3_annotation(path):
    """GFF3 with gene parents, the shape checkAnnotation emits.

    Every fourth feature deliberately deviates, to exercise the paths that
    otherwise never run:
      - no gene feature at all, so the CDS is an orphan whose Parent cannot be
        resolved
      - a gene without an old_locus_tag while the CDS carries one, so the
        fallback to the feature's own attribute is taken
    """
    lines = ["##gff-version 3"]
    for contig, feature, start, stop, strand, index in _features():
        orphan = index % 4 == 0
        gene_has_old_locus_tag = index % 4 != 1

        if not orphan:
            gene_attributes = f"ID=gene{index};locus_tag=b{index:04d};"
            if gene_has_old_locus_tag:
                gene_attributes += f"old_locus_tag=old{index:04d};"
            gene_attributes += f"Name=gene{index};gene=g{index};"
            lines.append("\t".join([
                contig, "RefSeq", "gene", str(start), str(stop), ".", strand, ".",
                gene_attributes,
            ]))

        lines.append("\t".join([
            contig, "RefSeq", feature, str(start), str(stop), ".", strand, "0",
            f"ID=cds{index};Parent=gene{index};locus_tag=b{index:04d};"
            f"old_locus_tag=cdsold{index:04d};"
            f"Name=gene{index};product=hypothetical protein {index};",
        ]))
    Path(path).write_text("\n".join(lines) + "\n")


def write_gtf2_annotation(path):
    """GTF2 with the gene_id "x"; attribute style, for the gff3 conversion."""
    lines = []
    for contig, feature, start, stop, strand, index in _features():
        for kind in ("gene", feature):
            lines.append("\t".join([
                contig, "Ensembl", kind, str(start), str(stop), ".", strand, "0",
                f'gene_id "gene{index}"; gene_name "name{index}"; '
                f'locus_tag "b{index:04d}"; gene_biotype "protein_coding";',
            ]))
    Path(path).write_text("\n".join(lines) + "\n")


def write_reparation_tracks(path, orf_type=""):
    """Reparation predictions, as create_reparation_gff.py writes them.

    Deliberately includes the same ORF reported by two replicates, so that the
    duplicate merging is exercised.

    `orf_type` defaults to empty, which is what Reparation reports when it has no
    type for an ORF and which writes "orf_type=;" into the attributes. That empty
    value used to crash reannotate_orfs.py.
    """
    lines = ["##gff-version 3"]
    for contig, feature, start, stop, strand, index in _features():
        if feature != "CDS":
            continue
        identifier = f"{contig}:{start}-{stop}:{strand}"
        for replicate in ("1", "2"):
            probability = 0.5 + (index % 40) / 100 + (0.01 if replicate == "2" else 0)
            lines.append("\t".join([
                contig, "reparation", "CDS", str(start), str(stop),
                f"{probability:.2f}", strand, "0",
                f"ID={identifier};Name={identifier};orf_type={orf_type};length={stop - start + 1};"
                f"ribo_count={100 + index};ribo_rpkm={index}.5;ribo_coverage=0.9;"
                f"sd_score=5.1;sd_pos=-8;prob={probability:.2f};reference=aTIS;"
                f"distance_from_atis=0;condition=A;replicate={replicate};method=reparation;",
            ]))
    Path(path).write_text("\n".join(lines) + "\n")


def write_deepribo_tracks(path):
    """DeepRibo predictions, as create_deepribo_gff.py writes them.

    Distance is a lowercase custom attribute; CDS phase remains strict GFF3.
    """
    lines = ["##gff-version 3"]
    for contig, feature, start, stop, strand, index in _features():
        if feature != "CDS":
            continue
        identifier = f"{contig}:{start}-{stop}:{strand}"
        for replicate in ("1", "2"):
            score = 1.0 + (index % 30) / 10 + (0.05 if replicate == "2" else 0)
            distance = 0 if index % 4 else -1
            lines.append("\t".join([
                contig, "deepribo", "CDS", str(start), str(stop),
                f"{score:.2f}", strand, "0",
                f"ID={identifier};pred_value={score:.2f};deepribo_distance={distance};"
                f"method=deepribo;condition=A;replicate={replicate};",
            ]))
    Path(path).write_text("\n".join(lines) + "\n")


def build(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    write_gff3_annotation(directory / "annotation_gff3.gff")
    write_gtf2_annotation(directory / "annotation_gtf2.gtf")
    # Two variants: one with an empty ORF_type, which is the case that used to
    # crash, and one with it filled in, used for the golden snapshots.
    write_reparation_tracks(directory / "reparation_tracks.gff")
    write_reparation_tracks(directory / "reparation_tracks_typed.gff", orf_type="sORF")
    write_deepribo_tracks(directory / "deepribo_tracks.gff")

    return directory


if __name__ == "__main__":
    import sys

    print(build(sys.argv[1]))
