"""
Contains scripts related to annotation processing.
Author: Rick Gelhausen
"""

import interlap
import pandas as pd
import lib.misc as misc

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

def retrieve_annotation_positions(annotation_file_path, read_intervals_dict, total_counts_dict, genome_length_dict, filtering_methods, mapping_method, rpkm_threshold,\
                                     overlap_distance, positions_out_ORF, positions_in_ORF):
    """
    Retrieve start/stop positions of annotated genes.
    Filter annotation based on:
        - gene distance
        - gene length
        - gene type
        - rpkm threshold
    """

    annotation_df = pd.read_csv(annotation_file_path, sep="\t", comment="#", header=None)
    annotation_intervals_dict = create_annotation_intervals_dict(annotation_df)

    start_codon_dict = {"-" : {}, "+" : {}}
    stop_codon_dict = {"-" : {}, "+" : {}}

    excluded_genes = { "overlap" : (0, []), "length" : (0, []), "rpkm" : (0, []), "type" : (0, []), "error" : (0, []) }
    included_genes = (0, [])
    for row in annotation_df.itertuples(index=False):
        chromosome = row[0]
        type = row[2]
        beginning = int(row[3]) - 1
        end = int(row[4]) - 1
        strand = row[6]

        # check type condition
        if type.lower() != "cds":
            excluded_genes["type"] = (excluded_genes["type"][0] + 1, excluded_genes["type"][1] + [row])
            continue

        if "overlap" in filtering_methods:
            # check overlap condition
            if (chromosome, strand) in annotation_intervals_dict:
                if len(list(annotation_intervals_dict[(chromosome, strand)].find((beginning-overlap_distance, end+overlap_distance)))) > 1:
                    excluded_genes["overlap"] = (excluded_genes["overlap"][0] + 1, excluded_genes["overlap"][1] + [row])
                    continue
            else:
                excluded_genes["error"] = (excluded_genes["error"][0] + 1, excluded_genes["error"][1] + [row])
                continue

        # check length condition
        gene_length = end - beginning + 1

        if "length" in filtering_methods:
            if gene_length < positions_in_ORF:
                excluded_genes["length"] = (excluded_genes["length"][0] + 1, excluded_genes["length"][1] + [row])
                continue

        # check rpkm condition
        if "rpkm" in filtering_methods:
            if (chromosome, strand) in read_intervals_dict:
                gene_read_counts = misc.count_reads(read_intervals_dict, chromosome, strand, beginning, end, mapping_method)
            else:
                excluded_genes["rpkm"] = (excluded_genes["rpkm"][0] + 1, excluded_genes["rpkm"][1] + [row])
                gene_read_counts = 0
                continue
            rpkm = misc.calculate_rpkm(gene_length, gene_read_counts, total_counts_dict[chromosome])
            if rpkm < rpkm_threshold:
                excluded_genes["rpkm"] = (excluded_genes["rpkm"][0] + 1, excluded_genes["rpkm"][1] + [row])
                continue

        # remove boundary cases
        if strand == "+":
            if beginning - positions_out_ORF < 0:
                continue
            if end + positions_out_ORF > genome_length_dict[chromosome]:
                continue

            if chromosome not in start_codon_dict[strand]:
                start_codon_dict[strand][chromosome] = [(beginning, beginning+2)]
            else:
                start_codon_dict[strand][chromosome].append((beginning, beginning+2))

            if chromosome not in stop_codon_dict[strand]:
                stop_codon_dict[strand][chromosome] = [(end-2, end)]
            else:
                stop_codon_dict[strand][chromosome].append((end-2, end))

        else:
            if beginning + positions_out_ORF < 0:
                continue
            if end - positions_out_ORF > genome_length_dict[chromosome]:
                continue

            if chromosome not in start_codon_dict[strand]:
                start_codon_dict[strand][chromosome] = [(end-2, end)]
            else:
                start_codon_dict[strand][chromosome].append((end-2, end))

            if chromosome not in stop_codon_dict[strand]:
                stop_codon_dict[strand][chromosome]  = [(beginning, beginning+2)]
            else:
                stop_codon_dict[strand][chromosome].append((beginning, beginning+2))

        included_genes = (included_genes[0] + 1, included_genes[1] + [row])

    print("Excluded genes:")
    for key, val in excluded_genes.items():
        print(f">>Entry removal based on {key}: {val[0]}")
    print(f"Included genes: {included_genes[0]}")

    return start_codon_dict, stop_codon_dict