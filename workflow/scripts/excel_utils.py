#!/usr/bin/env python

import sys
import pandas as pd
from collections import Counter, OrderedDict


import gff_utils

from Bio.Seq import Seq
from Bio import SeqIO

class OrderedCounter(Counter, OrderedDict):
    pass


# The attribute parsing lives in gff_utils, which is dependency free so that the
# scripts running in the mergetools environment can share it.
parse_attributes = gff_utils.parse_attributes
first_attribute = gff_utils.first_attribute


def get_te_header(wildcards):
    """
    generate the correct TE_header based on the available data
    """
    te_header = []
    te_header_dict = OrderedDict()
    for card in wildcards:
        method, condition, replicate = card.split("-")
        if method == "RIBO":
            if "%s-%s-%s" %("RNA", condition, replicate) in wildcards:
                if ("RIBO", condition) in  te_header_dict:
                    te_header_dict[("RIBO", condition)].append(replicate)
                else:
                    te_header_dict[("RIBO", condition)] = [replicate]
        elif method == "TIS":
            if "%s-%s-%s" %("RNATIS", condition, replicate) in wildcards:
                if ("TIS", condition) in  te_header_dict:
                    te_header_dict[("TIS", condition)].append(replicate)
                else:
                    te_header_dict[("TIS", condition)] = [replicate]
        elif method == "TTS":
            if "%s-%s-%s" %("RNATTS", condition, replicate) in wildcards:
                if ("TTS", condition) in  te_header_dict:
                    te_header_dict[("TTS", condition)].append(replicate)
                else:
                    te_header_dict[("TTS", condition)] = [replicate]

    for key, val in te_header_dict.items():
        method, condition = key
        if len(val) > 1:
            t_header = ["%s-%s-%s" % (method, condition, x) for x in val] + ["%s-%s-avg" % (method, condition)]
        else:
            t_header = ["%s-%s-%s" % (method, condition, x) for x in val]
        te_header.extend(t_header)

    return te_header

def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    if read_length == 0:
        print("read_length: 0 detected! Setting RPKM to 0!")
        return 0
    elif total_mapped == 0:
        print("total_mapped: 0 detected! Setting RPKM to 0!")
        return 0

    return float("%.2f" % ((read_count * 1000000000) / (total_mapped * read_length)))

def get_unique(in_list):
    seen = set()
    seen_add = seen.add
    return sorted(
        [x for x in in_list if not (x in seen or seen_add(x))],
        key=lambda value: value.lower() if isinstance(value, str) else value,
    )

def retrieve_column_information(attributes):
    """
    check for gff2/gff3 format and generate a list of information for the final tables
    [pred_value, name, product, note, evidence, locus_tag, old_locus_tag]
    """

    if "ORF_type=;" in attributes:
        attributes = attributes.replace("ORF_type=;", "")

    parsed = parse_attributes(attributes)

    return [
        first_attribute(parsed, "pred_value", "prob"),
        first_attribute(parsed, "name"),
        first_attribute(parsed, "product"),
        first_attribute(parsed, "note"),
        first_attribute(parsed, "evidence"),
        first_attribute(parsed, "locus_tag"),
        first_attribute(parsed, "old_locus_tag"),
    ]

def get_genome_information(genome, start, stop, strand):
    """
    retrieve the nucleotide sequence and amino acid sequence
    and the start and stop codons
    """
    if strand == "+":
        nucleotide_seq = genome[0][start:stop+1]
        nt_window = genome[0][start-15:start]
    else:
        nucleotide_seq = genome[1][start:stop+1][::-1]
        nt_window = str(Seq(genome[0][stop+1:stop+16]).reverse_complement())

    start_codon = nucleotide_seq[0:3]
    stop_codon = nucleotide_seq[-3:]

    coding_dna = Seq(nucleotide_seq)
    if len(coding_dna) % 3 != 0:
        aa_seq = ""
    else:
        aa_seq = str(coding_dna.translate(table=11,to_stop=True))

    return start_codon, stop_codon, nucleotide_seq, aa_seq, nt_window

def excel_writer(output_path, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Note", "Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"] + [card + "_rpkm" for card in wildcards]
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        worksheet.freeze_panes(1, 0)
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only or series.empty:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            #print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.close()

def te(ribo_count, rna_count):
    """
    calculate the translational efficiency for one entry
    """

    if ribo_count == 0 and rna_count == 0:
        return "NaN"
    elif rna_count == 0:
        return "NaN"
    else:
        return ribo_count / rna_count

def get_avg(t_eff):
    """
    get the final TE list
    """

    valid_count = 0
    sum = 0
    for t in t_eff:
        if t != "NaN":
            valid_count += 1
            sum += t

    if valid_count == 0:
        t_eff.extend(["NaN"])

    else:
        t_eff.extend([sum / valid_count])

    return t_eff

def calculate_te(read_list, wildcards, conditions):
    """
    calculate the translational efficiency
    """
    read_dict = OrderedDict()
    te_dict = OrderedDict()
    for idx in range(len(wildcards)):
        method, condition, replicate = wildcards[idx].split("-")
        key = (method, condition, replicate)
        if key not in read_dict:
            read_dict[key] = read_list[idx]
        else:
            print("warning: multiple equal keys")

    te_list = []
    for key, val in read_dict.items():
        method, condition, replicate = key
        if method == "RIBO":
            if ("RNA", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNA", condition, replicate)]
                cur_te = te(rpkm_ribo, rpkm_rna)
                if ("RIBO", condition) in te_dict:
                    te_dict[("RIBO", condition)].append(cur_te)
                else:
                    te_dict[("RIBO", condition)] = [cur_te]

        elif method == "TIS":
            if ("RNATIS", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNATIS", condition, replicate)]
                cur_te = te(rpkm_ribo, rpkm_rna)
                if ("TIS", condition) in te_dict:
                    te_dict[("TIS", condition)].append(cur_te)
                else:
                    te_dict[("TIS", condition)] = [cur_te]

        elif method == "TTS":
            if ("RNATTS", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNATTS", condition, replicate)]
                cur_te = te(rpkm_ribo, rpkm_rna)
                if ("TTS", condition) in te_dict:
                    te_dict[("TTS", condition)].append(cur_te)
                else:
                    te_dict[("TTS", condition)] = [cur_te]

    te_list = []
    for key, val in te_dict.items():
        if len(val) > 1:
            t_eff = get_avg(val)
        else:
            t_eff = val
        te_list.extend(t_eff)

    return te_list

# Which columns each differential expression tool contributes to the overview
# table. The tools differ only in these names, not in how the file is read.
DIFFEX_COLUMNS = {
    "riborex": ("log2FoldChange", "pvalue", "padj"),
    "xtail": ("log2FC_TE_final", "pvalue_final", "pvalue_adjust"),
    "deltate": (
        "RIBO_log2FoldChange", "RIBO_pvalue", "RIBO_padj",
        "RNA_log2FoldChange", "RNA_pvalue", "RNA_padj",
        "TE_log2FoldChange", "TE_pvalue", "TE_padj",
    ),
}


def generate_diffex_dict(path, tool):
    """{(gene_id, contrast): (values...)} for one differential expression tool.

    The contrast column is written as "contrast_<name>" upstream, so only the
    part after the underscore is kept.
    """
    frame = pd.read_csv(path, sep=",", comment="#")
    columns = DIFFEX_COLUMNS[tool]

    result = {}
    for row in frame.itertuples(index=False, name="Pandas"):
        gene_id = getattr(row, "gene_id")
        contrast = getattr(row, "contrast").split("_")[1]
        result[(gene_id, contrast)] = tuple(getattr(row, column) for column in columns)

    return result


def generate_riborex_dict(riborex_path):
    return generate_diffex_dict(riborex_path, "riborex")


def generate_xtail_dict(xtail_path):
    return generate_diffex_dict(xtail_path, "xtail")


def generate_deltate_dict(deltate_path):
    return generate_diffex_dict(deltate_path, "deltate")


def _prediction_rows(path):
    """Yield (identifier, row, parsed attributes, read counts) for a prediction GFF."""
    try:
        frame = pd.read_csv(path, header=None, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        return
    prefix_columns = 9

    for row in frame.itertuples(index=False, name="Pandas"):
        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        parsed = parse_attributes(getattr(row, "_8"))
        read_list = [getattr(row, "_%s" % x) for x in range(prefix_columns, len(row))]

        identifier = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
        yield identifier, row, parsed, read_list


def generate_reparation_dict(reparation_path):
    """
    create a dictionary containing all important reparation input
    """
    return {
        identifier: (
            first_attribute(parsed, "prob"),
            first_attribute(parsed, "evidence"),
            read_list,
        )
        for identifier, row, parsed, read_list in _prediction_rows(reparation_path)
    }


def generate_deepribo_dict(deepribo_path):
    """
    create a dictionary containing all important deepribo input
    """
    deepribo_dict = {}
    for identifier, row, parsed, read_list in _prediction_rows(deepribo_path):
        pred_value = first_attribute(parsed, "pred_value")
        # Entries without a prediction value carry no DeepRibo result to report.
        if pred_value == "":
            continue

        deepribo_dict[identifier] = (
            getattr(row, "_5"),
            pred_value,
            first_attribute(parsed, "evidence"),
            read_list,
        )

    return deepribo_dict


def _read_annotation_entries(annotation_path, keep_feature, default):
    """Split an annotation into feature entries and their gene-level parents.

    `keep_feature` selects which features become entries; everything tagged gene
    or pseudogene is collected separately so that a missing locus tag can fall
    back to the gene's.

    `default` is what a missing attribute becomes. It differs between callers --
    "" for the CDS view and None for the non-CDS view -- and that difference is
    load bearing: the fallback below tests against "", so with None defaults a
    missing locus tag is never filled in from the gene. Preserved as it was.
    """
    frame = pd.read_csv(annotation_path, sep="\t", comment="#", header=None)

    entries = {}
    genes = {}
    for row in frame.itertuples(index=False, name="Pandas"):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        parsed = parse_attributes(getattr(row, "_8"))
        read_list = [getattr(row, "_%s" % x) for x in range(9, len(row))]

        key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)

        if keep_feature(feature):
            entries[key] = {
                "feature": feature,
                "gene_id": first_attribute(parsed, "gene_id", "id", default=default),
                "locus_tag": first_attribute(parsed, "locus_tag", default=default),
                "old_locus_tag": first_attribute(parsed, "old_locus_tag", default=default),
                "name": first_attribute(parsed, "name", "gene_name", default=default),
                "product": first_attribute(parsed, "product", default=default),
                "note": first_attribute(parsed, "note", default=default),
                "read_list": read_list,
            }
        elif feature.lower() in ["gene", "pseudogene"]:
            genes[key] = (
                first_attribute(parsed, "name", "gene_name"),
                first_attribute(parsed, "locus_tag", "gene_id"),
                first_attribute(parsed, "old_locus_tag"),
            )

    return entries, genes


def _apply_gene_fallback(entry, genes, key):
    """Fill a missing locus tag from the gene feature at the same coordinates."""
    gene_name = ""
    locus_tag = entry["locus_tag"]
    old_locus_tag = entry["old_locus_tag"]

    if key in genes:
        gene_name, gene_locus_tag, gene_old_locus_tag = genes[key]
        if locus_tag == "":
            locus_tag = gene_locus_tag
        if old_locus_tag == "":
            old_locus_tag = gene_old_locus_tag

    return gene_name, locus_tag, old_locus_tag


def generate_annotation_dict(annotation_path):
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)
    """
    coding = ["cds", "srna", "rrna", "trna", "ncrna"]
    entries, genes = _read_annotation_entries(
        annotation_path, lambda feature: feature.lower() in coding, default=""
    )

    annotation_meta_dict = {}
    for key, entry in entries.items():
        gene_name, locus_tag, old_locus_tag = _apply_gene_fallback(entry, genes, key)
        annotation_meta_dict[key] = (
            entry["gene_id"], locus_tag, entry["name"], gene_name,
            old_locus_tag, entry["product"], entry["note"], entry["read_list"],
        )

    return annotation_meta_dict


def generate_non_cds_dict(annotation_path):
    """
    create dictionary from annotation ignoring cds.
    key : (gene_id, locus_tag, name, gene_name)
    """
    excluded = ["cds", "gene", "pseudogene", "exon"]
    entries, genes = _read_annotation_entries(
        annotation_path, lambda feature: feature.lower() not in excluded, default=None
    )

    annotation_meta_dict = {}
    for key, entry in entries.items():
        gene_name, locus_tag, old_locus_tag = _apply_gene_fallback(entry, genes, key)
        annotation_meta_dict[key] = (
            entry["feature"], entry["gene_id"], locus_tag, entry["name"], gene_name,
            old_locus_tag, entry["product"], entry["note"], entry["read_list"],
        )

    return annotation_meta_dict


def annotation_to_dict(annotation_file):
    """
    read the annotation file into a dictionary
    { ID : (chromosome,start,stop,strand,attributes) }
    """

    annotation_df = pd.read_csv(annotation_file, sep="\t", comment="#", header=None)

    # Gene and pseudogene features carry the locus tags; the CDS-like features
    # that follow reference them through their Parent attribute.
    parent_dict = {}
    for row in annotation_df.itertuples(index=False, name="Pandas"):
        if getattr(row, "_2").lower() not in ["gene", "pseudogene"]:
            continue

        parsed = parse_attributes(getattr(row, "_8"))
        if "id" not in parsed:
            continue

        parent_dict[parsed["id"]] = (
            first_attribute(parsed, "locus_tag"),
            first_attribute(parsed, "old_locus_tag"),
            first_attribute(parsed, "gene"),
        )

    annotation_dict = {}
    for row in annotation_df.itertuples(index=False, name="Pandas"):
        if getattr(row, "_2").lower() not in ["cds", "srna", "rrna", "ncrna", "trna"]:
            continue

        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        parsed = parse_attributes(getattr(row, "_8"))

        parent = first_attribute(parsed, "parent")
        identifier = f"{chromosome}:{start}-{stop}:{strand}"

        if parent == "" or parent not in parent_dict:
            # Without a resolvable parent the feature keeps only its own name.
            name = first_attribute(parsed, "name")
            annotation_dict[identifier] = (chromosome, start, stop, strand, name, "", "")
        else:
            locus_tag, old_locus_tag, name = parent_dict[parent]
            if name == "":
                name = first_attribute(parsed, "name")
            annotation_dict[identifier] = (
                chromosome, start, stop, strand, name, locus_tag, old_locus_tag
            )

    return annotation_dict



# --------------------------------------------------------------------------
# Shared table building
#
# generate_excel.py, generate_excel_reparation.py and generate_excel_deepribo.py
# differed only in which columns they emit; everything below the column list --
# loading the genome, the mapped read totals, the RPKM and translational
# efficiency arithmetic -- was copied between them.
# --------------------------------------------------------------------------


class TableContext:
    """Everything the per-row work needs that does not change between rows."""

    def __init__(self, genome_path, total_mapped_path):
        self.genome = {}
        for entry in SeqIO.parse(genome_path, "fasta"):
            self.genome[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

        self.total_mapped = {}
        wildcards = []
        with open(total_mapped_path) as handle:
            for line in handle:
                wildcard, chromosome, value = line.strip().split("\t")
                self.total_mapped[(wildcard, chromosome)] = int(value)
                wildcards.append(wildcard)

        self.wildcards = get_unique(wildcards)
        self.te_header = get_te_header(self.wildcards)
        self.conditions = get_unique([card.split("-")[1] for card in self.wildcards])

    def rpkm_columns(self):
        return ["%s_rpkm" % card for card in self.wildcards]

    def te_columns(self):
        return ["%s_TE" % condition for condition in self.te_header]


class Row:
    """One annotation row, with everything the column definitions may need."""

    def __init__(self, raw, context, prefix_columns):
        self.raw = raw
        self.chromosome = getattr(raw, "_0")
        self.source = getattr(raw, "_1")
        self.feature = getattr(raw, "_2")
        self.start = getattr(raw, "_3")
        self.stop = getattr(raw, "_4")
        self.score = getattr(raw, "_5")
        self.strand = getattr(raw, "_6")
        self.phase = getattr(raw, "_7")
        self.attributes = getattr(raw, "_8")
        parsed_attributes = parse_attributes(self.attributes)
        self.novel_rank = first_attribute(parsed_attributes, "novel_rank")
        if self.novel_rank == "":
            # Compatibility with mapped DeepRibo files produced before the
            # novel rank moved out of the GFF3 phase column.
            self.novel_rank = self.phase
        else:
            try:
                self.novel_rank = int(self.novel_rank)
            except ValueError:
                pass

        self.identifier = "%s:%s-%s:%s" % (
            self.chromosome, self.start, self.stop, self.strand
        )
        self.length = self.stop - self.start + 1
        self.codon_count = int(self.length / 3)

        (self.pred_value, self.name, self.product, self.note,
         self.evidence, self.locus_tag, self.old_locus_tag) = retrieve_column_information(
            self.attributes
        )

        self.start_codon, self.stop_codon, self.nucleotide_seq, self.aa_seq, self.nt_window = (
            "", "", "", "", ""
        )
        if self.chromosome in context.genome:
            (self.start_codon, self.stop_codon, self.nucleotide_seq,
             self.aa_seq, self.nt_window) = get_genome_information(
                context.genome[self.chromosome], self.start - 1, self.stop - 1, self.strand
            )

        self.read_list = [
            getattr(raw, "_%s" % x) for x in range(prefix_columns, len(raw))
        ]
        self.rpkm_list = []
        for index, value in enumerate(self.read_list):
            key = (context.wildcards[index], self.chromosome)
            if key not in context.total_mapped:
                self.rpkm_list.append(0)
            else:
                self.rpkm_list.append(
                    calculate_rpkm(context.total_mapped[key], value, self.length)
                )

        self.te_list = calculate_te(self.rpkm_list, context.wildcards, context.conditions)


def build_annotation_table(reads_path, context, columns, source=None):
    """Build one DataFrame from a read-count annotation.

    `columns` is a list of (header, accessor) pairs. An accessor is either an
    attribute name on Row, or a callable taking the Row. Two names are expanded
    in place into one column per library: "te_list" and "rpkm_list".

    `source` overrides the Source column when a script writes a fixed value
    rather than passing the file's own second column through.
    """
    header = []
    for name, _ in columns:
        if name == "te_list":
            header.extend(context.te_columns())
        elif name == "rpkm_list":
            header.extend(context.rpkm_columns())
        else:
            header.append(name)

    try:
        read_df = pd.read_csv(reads_path, comment="#", header=None, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=header), []
    prefix_columns = len(read_df.columns) - len(context.wildcards)

    records = []
    for raw in read_df.itertuples(index=False, name="Pandas"):
        row = Row(raw, context, prefix_columns)
        if source is not None:
            row.source = source

        values = []
        for name, accessor in columns:
            value = accessor(row) if callable(accessor) else getattr(row, accessor)
            if name in ("te_list", "rpkm_list"):
                values.extend(value)
            else:
                values.append(value)
        records.append((row, values))

    frame = pd.DataFrame.from_records([values for _, values in records], columns=header)
    return frame, [row for row, _ in records]


# Columns every one of the three tables shares, in the order they appear.
def identity_columns(identifier_header="Identifier"):
    return [
        (identifier_header, "identifier"),
        ("Genome", "chromosome"),
        ("Source", "source"),
        ("Feature", "feature"),
        ("Start", "start"),
        ("Stop", "stop"),
        ("Strand", "strand"),
    ]


LOCUS_COLUMNS = [
    ("Locus_tag", "locus_tag"),
    ("Old_locus_tag", "old_locus_tag"),
    ("Name", "name"),
    ("Length", "length"),
    ("Codon_count", "codon_count"),
]

MEASURE_COLUMNS = [("te_list", "te_list"), ("rpkm_list", "rpkm_list")]

SEQUENCE_COLUMNS = [
    ("Start_codon", "start_codon"),
    ("Stop_codon", "stop_codon"),
    ("Upstream_15nt", "nt_window"),
    ("Nucleotide_seq", "nucleotide_seq"),
    ("Aminoacid_seq", "aa_seq"),
]


# --------------------------------------------------------------------------
# Differential expression tables
#
# generate_excel_riborex.py, generate_excel_xtail.py and generate_excel_deltate.py
# shared everything except which statistics columns they carry and which of
# those drives the sorting and the up/down split.
# --------------------------------------------------------------------------


def resolve_location(unique_id, annotation_dict):
    """Coordinates and names for one feature identifier.

    Identifiers absent from the annotation are predicted ORFs, whose coordinates
    are encoded in the identifier itself as <genome>:<start>-<stop>:<strand>.
    """
    if unique_id in annotation_dict:
        return annotation_dict[unique_id]

    if ":" in str(unique_id) and "-" in str(unique_id):
        chromosome, section, strand = unique_id.split(":")
        start, stop = section.split("-")
        return chromosome, start, stop, strand, "", "", ""

    sys.exit("Error... ID is not novel and not in the annotation!")


# The same leading block as the feature tables, so a reader moving between
# workbooks finds the identifying columns in the same place.
DIFFEX_IDENTITY_HEADER = [
    "Identifier", "Genome", "Start", "Stop", "Strand",
    "Locus_tag", "Old_locus_tag", "Name",
]


def build_diffex_table(frame, annotation_dict, genome_dict, statistics, identifier_field):
    """One differential expression sheet.

    `statistics` is a list of (header, field) pairs naming the tool's own
    columns; `identifier_field` is the attribute holding the feature identifier,
    which differs because the tools' CSVs are written by different code.
    """
    header = (
        DIFFEX_IDENTITY_HEADER
        + [name for name, _ in statistics]
        + ["Length", "Codon_count", "Start_codon", "Stop_codon", "Nucleotide_seq", "Aminoacid_seq"]
    )

    records = []
    for row in frame.itertuples(index=False, name="Pandas"):
        unique_id = getattr(row, identifier_field)
        chromosome, start, stop, strand, gene_name, locus_tag, old_locus_tag = resolve_location(
            unique_id, annotation_dict
        )

        start = int(start)
        stop = int(stop)
        length = stop - start + 1
        codon_count = int(length / 3)

        start_codon, stop_codon, nucleotide_seq, aa_seq = "", "", "", ""
        if chromosome in genome_dict:
            start_codon, stop_codon, nucleotide_seq, aa_seq, _ = get_genome_information(
                genome_dict[chromosome], start - 1, stop - 1, strand
            )

        records.append(
            [unique_id, chromosome, start, stop, strand, locus_tag, old_locus_tag, gene_name]
            + [getattr(row, field) for _, field in statistics]
            + [length, codon_count, start_codon, stop_codon, nucleotide_seq, aa_seq]
        )

    return pd.DataFrame.from_records(records, columns=header)


def split_up_down(all_df, log2fc_column, padj_column, log2fc_cutoff, padj_cutoff):
    """Sort, then split into the significantly up- and down-regulated sheets."""
    all_df = all_df.sort_values(
        by=[padj_column, "Genome", "Start", "Stop", "Strand"], kind="stable"
    )
    significant = all_df[padj_column] <= padj_cutoff
    return {
        "all": all_df,
        "TE_up": all_df[(all_df[log2fc_column] >= log2fc_cutoff) & significant],
        "TE_down": all_df[(all_df[log2fc_column] <= log2fc_cutoff * -1) & significant],
    }


def read_genome_dict(genome_path):
    """{sequence id: (forward, complement)} for the sequence lookups."""
    return {
        str(entry.id): (str(entry.seq), str(entry.seq.complement()))
        for entry in SeqIO.parse(genome_path, "fasta")
    }
