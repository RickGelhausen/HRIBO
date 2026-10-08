"""
Contains scripts related to annotation processing.
Author: Rick Gelhausen
"""

import interlap
import pandas as pd
import lib.misc as misc


SELECTION_COLUMNS = [
    "contig", "start", "end", "strand", "feature_id", "length_nt", "status", "reason"
]


def _cds_feature_id(row):
    """Prefer GFF ID, then locus_tag, with a stable coordinate fallback."""
    attributes = {}
    if len(row) > 8 and not pd.isna(row[8]):
        for field in str(row[8]).split(";"):
            key, separator, value = field.strip().partition("=")
            if separator:
                attributes.setdefault(key, value.strip())
    for key in ("ID", "locus_tag"):
        if attributes.get(key) not in (None, "", "."):
            return attributes[key]
    return f"{row[0]}:{int(row[3])}-{int(row[4])}:{row[6]}"


def create_annotation_intervals_dict(annotation_df):
    """
    Create interlap instances for each chrom / strand in the annotation
    """

    annotation_intervals_dict = {}

    tmp_dict = {}
    for row in annotation_df.itertuples(index=False):
        chromosome = row[0]
        beginning = int(row[3]) - 1
        end = int(row[4]) - 1
        strand = row[6]
        feature = row[2]

        if feature.lower() != "cds":
            continue

        if (chromosome, strand) not in tmp_dict:
            tmp_dict[(chromosome, strand)] = [(beginning, end)]
        else:
            tmp_dict[(chromosome, strand)].append((beginning, end))

    for key, val in tmp_dict.items():
        inter = interlap.InterLap()
        inter.update(val)
        annotation_intervals_dict[key] = inter

    return annotation_intervals_dict


def metagene_window_bounds(
    beginning,
    end,
    strand,
    positions_out_ORF,
    positions_in_ORF,
):
    """Return inclusive start/stop-profile windows in genomic coordinates.

    ``beginning`` and ``end`` are already zero-based and inclusive.  The
    outside flank points away from the ORF, while the inside flank points into
    it, so the two anchors exchange their genomic geometry on the minus strand.
    Keeping both windows explicit makes asymmetric inside/outside settings and
    the last valid contig coordinate (``genome_length - 1``) unambiguous.
    """
    if strand == "+":
        return {
            "start": (
                beginning - positions_out_ORF,
                beginning + positions_in_ORF - 1,
            ),
            "stop": (
                end - positions_in_ORF + 1,
                end + positions_out_ORF,
            ),
        }

    return {
        "start": (
            end - positions_in_ORF + 1,
            end + positions_out_ORF,
        ),
        "stop": (
            beginning - positions_out_ORF,
            beginning + positions_in_ORF - 1,
        ),
    }


def metagene_windows_fit_contig(windows, genome_length):
    """Whether every inclusive profile window lies on a zero-based contig."""
    last_position = genome_length - 1
    return all(
        window_start >= 0 and window_stop <= last_position
        for window_start, window_stop in windows.values()
    )


def retrieve_annotation_positions(
    annotation_file_path,
    read_intervals_dict,
    total_counts_dict,
    genome_length_dict,
    filtering_methods,
    mapping_method,
    rpkm_threshold,
    overlap_distance,
    positions_out_ORF,
    positions_in_ORF,
    length_cutoff=None,
    *,
    cds_max_length=None,
    return_selection=False,
    required_anchors=("start", "stop"),
):
    """
    Retrieve start/stop positions of annotated genes.
    Filter annotation based on:
        - gene distance
        - gene length (the larger of the in-ORF window and length cutoff)
        - gene type
        - rpkm threshold

    ``cds_max_length`` selects CDSs strictly shorter than the specified number
    of nucleotides, before applying the configured filters.  Neighboring CDSs
    are always indexed from the full annotation.  The optional selection table
    uses zero-based, inclusive ``start`` and ``end`` coordinates and contains
    the first exclusion reason for every CDS.  In selection mode, exact CDS
    duplicates (contig, coordinates, strand) are profiled and reported once,
    with their distinct identifiers joined by commas.  The default two-value
    return and legacy duplicate handling are preserved for existing callers.
    ``required_anchors`` controls which profile windows must fit the contig;
    both are required by default, while start-only analyses can require only
    the start window.
    """

    required_anchors = tuple(required_anchors)
    if not required_anchors or any(anchor not in ("start", "stop") for anchor in required_anchors):
        raise ValueError("required_anchors must contain 'start' and/or 'stop'")

    try:
        annotation_df = pd.read_csv(annotation_file_path, sep="\t", comment="#", header=None)
    except pd.errors.EmptyDataError:
        annotation_df = pd.DataFrame(columns=range(9))
    annotation_intervals_dict = create_annotation_intervals_dict(annotation_df)

    feature_ids = {}
    if return_selection:
        for row in annotation_df.itertuples(index=False):
            if str(row[2]).lower() != "cds":
                continue
            key = (row[0], int(row[3]) - 1, int(row[4]) - 1, row[6])
            identifiers = feature_ids.setdefault(key, [])
            identifier = _cds_feature_id(row)
            if identifier not in identifiers:
                identifiers.append(identifier)

    library_total = None
    if "rpkm" in filtering_methods:
        # RPKM describes abundance within the complete library.  A contig-local
        # denominator makes otherwise identical genes pass or fail depending on
        # which replicon they happen to occupy.
        library_total = misc.library_read_total(total_counts_dict)

    start_codon_dict = {"-" : {}, "+" : {}}
    stop_codon_dict = {"-" : {}, "+" : {}}

    excluded_genes = dict.fromkeys(
        ("cohort", "overlap", "length", "rpkm", "boundary", "type", "error"), 0
    )
    included_genes = 0
    selection_records = []
    seen_coordinates = set()

    def exclude(reason, record):
        excluded_genes[reason] += 1
        if return_selection:
            selection_records.append({**record, "status": "excluded", "reason": reason})

    for row in annotation_df.itertuples(index=False):
        chromosome = row[0]
        feature_type = row[2]
        beginning = int(row[3]) - 1
        end = int(row[4]) - 1
        strand = row[6]

        # check type condition
        if str(feature_type).lower() != "cds":
            excluded_genes["type"] += 1
            continue

        coordinate_key = (chromosome, beginning, end, strand)
        if return_selection:
            if coordinate_key in seen_coordinates:
                continue
            seen_coordinates.add(coordinate_key)

        gene_length = end - beginning + 1
        record = {
            "contig": chromosome,
            "start": beginning,
            "end": end,
            "strand": strand,
            "feature_id": ",".join(feature_ids[coordinate_key]) if return_selection else "",
            "length_nt": gene_length,
        }

        if cds_max_length is not None and gene_length >= cds_max_length:
            exclude("cohort", record)
            continue

        if "overlap" in filtering_methods:
            # check overlap condition
            if (chromosome, strand) in annotation_intervals_dict:
                neighbors = list(annotation_intervals_dict[(chromosome, strand)].find(
                    (beginning-overlap_distance, end+overlap_distance)
                ))
                if return_selection:
                    # Duplicate annotations of this CDS are not another gene.
                    neighbors = set(neighbors)
                if len(neighbors) > 1:
                    exclude("overlap", record)
                    continue
            else:
                exclude("error", record)
                continue

        # check length condition
        if "length" in filtering_methods:
            minimum_gene_length = max(positions_in_ORF, length_cutoff or 0)
            if gene_length < minimum_gene_length:
                exclude("length", record)
                continue

        # check rpkm condition
        if "rpkm" in filtering_methods:
            if (chromosome, strand) in read_intervals_dict:
                gene_read_counts = misc.count_reads(read_intervals_dict, chromosome, strand, beginning, end, mapping_method)
            else:
                exclude("rpkm", record)
                continue
            rpkm = misc.calculate_rpkm(
                gene_length, gene_read_counts, library_total
            )
            if rpkm < rpkm_threshold:
                exclude("rpkm", record)
                continue

        # Every required profile window must lie on the contig; start-only
        # analyses need not support an undisplayed stop window.
        if chromosome not in genome_length_dict or strand not in ("+", "-"):
            exclude("error", record)
            continue
        profile_windows = metagene_window_bounds(
            beginning,
            end,
            strand,
            positions_out_ORF,
            positions_in_ORF,
        )
        if not metagene_windows_fit_contig(
            {anchor: profile_windows[anchor] for anchor in required_anchors},
            genome_length_dict[chromosome],
        ):
            exclude("boundary", record)
            continue

        if strand == "+":
            if chromosome not in start_codon_dict[strand]:
                start_codon_dict[strand][chromosome] = [(beginning, beginning+2)]
            else:
                start_codon_dict[strand][chromosome].append((beginning, beginning+2))

            if chromosome not in stop_codon_dict[strand]:
                stop_codon_dict[strand][chromosome] = [(end-2, end)]
            else:
                stop_codon_dict[strand][chromosome].append((end-2, end))

        else:
            if chromosome not in start_codon_dict[strand]:
                start_codon_dict[strand][chromosome] = [(end-2, end)]
            else:
                start_codon_dict[strand][chromosome].append((end-2, end))

            if chromosome not in stop_codon_dict[strand]:
                stop_codon_dict[strand][chromosome]  = [(beginning, beginning+2)]
            else:
                stop_codon_dict[strand][chromosome].append((beginning, beginning+2))

        included_genes += 1
        if return_selection:
            selection_records.append({**record, "status": "retained", "reason": ""})

    print("Excluded genes:")
    for key, val in excluded_genes.items():
        print(f">>Entry removal based on {key}: {val}")
    print(f"Included genes: {included_genes}")

    if return_selection:
        return start_codon_dict, stop_codon_dict, pd.DataFrame(
            selection_records, columns=SELECTION_COLUMNS
        )
    return start_codon_dict, stop_codon_dict
