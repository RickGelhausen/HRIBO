import numpy as np
import pysam

from lib import misc

class Feature:
    def __init__(self, gfftype, seqid, start, end, strand):
        self.gfftype = gfftype
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand

    def __str__(self):
        return self.gfftype + "_" + self.seqid + "_" + self.start + "_" + self.end + "_" + self.strand



def init_write_wig(file_handle,library_name, strand):
    tid = "%s_%s" % (library_name, strand)
    file_handle.write(("track type=wiggle_0 name=\"%s\"\n" % (tid)))


def write_wig(file_handle, seq_id, mappings, discard_zeros=True, factor=1.0):
    file_handle.write("variableStep chrom=%s span=1\n" % (seq_id))
    file_handle.write("\n".join(["%s %s" % (pos + 1, mapping * factor) for pos, mapping in filter(lambda pos_and_cov: pos_and_cov[1] != 0.0, enumerate(mappings))]) + "\n")


def add_aln_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    # ``reference_start:reference_end`` paints across deletions and spliced
    # reference skips.  Coverage belongs only to CIGAR-aligned reference blocks.
    for block_start, block_end in entry.get_blocks():
        if ((entry.is_reverse is False and entry.is_read2 is False) or
                (entry.is_reverse is True and entry.is_read2 is True)):
            if strand_swap:
               mappings["reverse"][block_start:block_end] -= increment
            else:
               mappings["forward"][block_start:block_end] += increment
        else:
            if strand_swap:
               mappings["forward"][block_start:block_end] += increment
            else:
               mappings["reverse"][block_start:block_end] -= increment


def add_first_base_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    if ((entry.is_reverse is False and entry.is_read2 is False) or
            (entry.is_reverse is True and entry.is_read2 is True)):
        if strand_swap:
            mappings["reverse"][end-1] -= increment
        else:
            mappings["forward"][start] += increment
    else:
        if strand_swap:
            mappings["forward"][start] += increment
        else:

            mappings["reverse"][end-1] -= increment

def add_last_base_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    if ((entry.is_reverse is False and entry.is_read2 is False) or
            (entry.is_reverse is True and entry.is_read2 is True)):
        if strand_swap:
            mappings["reverse"][start] -= increment
        else:
            mappings["forward"][end-1] += increment
    else:
        if strand_swap:
            mappings["forward"][end-1] += increment
        else:
            mappings["reverse"][start] -= increment


def add_centered_mapping(mappings, strand_swap, entry, increment, start, end, clip_length):
    center_start = start + clip_length
    center_end = end - clip_length
    center_blocks = [
        (max(block_start, center_start), min(block_end, center_end))
        for block_start, block_end in entry.get_blocks()
        if block_start < center_end and block_end > center_start
    ]
    center_length = float(
        sum(block_end - block_start for block_start, block_end in center_blocks)
    )
    if center_length < 1.0:
        return
    for block_start, block_end in center_blocks:
        if ((entry.is_reverse is False and entry.is_read2 is False) or
                (entry.is_reverse is True and entry.is_read2 is True)):
            if strand_swap:
                mappings["reverse"][block_start:block_end] -= (
                    increment / center_length
                )
            else:
                mappings["forward"][block_start:block_end] += (
                    increment / center_length
                )
        else:
            if strand_swap:
                mappings["forward"][block_start:block_end] += (
                    increment / center_length
                )
            else:
                mappings["reverse"][block_start:block_end] -= (
                    increment / center_length
                )


def wig_file_path_setter(path, library_name, normtype, strand):
    resultpath = (path + normtype + "/" + library_name + "." + normtype + "." + strand + ".wig")
    return(resultpath)


def select_mapping_add_function(mapping_style):
    if mapping_style == "first_base_only":
        return add_first_base_mapping
    elif mapping_style == "last_base_only":
        return add_last_base_mapping
    elif mapping_style == "centered":
        return add_centered_mapping
    else:
        return add_aln_mapping


def compute_seqid_mapping(bam_path, read_count_splitting, mapping_add_function, strand_swap, clip_length):
    bam = pysam.Samfile(bam_path)
    results = []
    for ref_seq, length in zip(bam.references, bam.lengths):
        mappings = {}
        for strand in ["forward", "reverse"]:
            mappings[strand] = np.array([0.0] * length)
        for entry in bam.fetch(ref_seq):
            number_of_hits = dict(entry.tags)["NH"]
            start = entry.pos
            end = entry.aend
            if read_count_splitting is True:
                increment = 1.0 / float(number_of_hits)
            else:
                increment = 1.0
            mapping_add_function(mappings, strand_swap, entry, increment, start, end, clip_length)
        results.append((ref_seq, mappings))
    return results


def compute_wig(bam_path, wig_file_path, library_name, library_read_total,
                minimum_library_total, read_count_splitting=True,
                mapping_style="centered", clip_length=11,
                non_strand_specific=False, strand_swap=False):
    if library_read_total <= 0:
        raise ValueError(
            f"library {library_name!r} has no mapped reads; track normalisation "
            "is undefined"
        )
    if minimum_library_total <= 0:
        raise ValueError(
            "the minimum library-wide mapped-read total must be positive"
        )

    mapping_add_function = select_mapping_add_function(mapping_style)
    seqid_mapping = compute_seqid_mapping(bam_path, read_count_splitting, mapping_add_function, strand_swap, clip_length)
    raw_fw_handle=open(wig_file_path_setter(wig_file_path, library_name, "raw", "forward"), "w")
    init_write_wig(raw_fw_handle,library_name,"forward")
    raw_rw_handle=open(wig_file_path_setter(wig_file_path, library_name, "raw", "reverse"), "w")
    init_write_wig(raw_rw_handle,library_name,"reverse")
    min_fw_handle=open(wig_file_path_setter(wig_file_path, library_name, "min", "forward"), "w")
    init_write_wig(min_fw_handle, library_name, "forward")
    min_rw_handle=open(wig_file_path_setter(wig_file_path, library_name, "min" , "reverse"), "w")
    init_write_wig(min_rw_handle,library_name,"reverse")
    mil_fw_handle=open(wig_file_path_setter(wig_file_path, library_name, "mil", "forward"), "w")
    init_write_wig(mil_fw_handle, library_name, "forward")
    mil_rw_handle=open(wig_file_path_setter(wig_file_path, library_name, "mil", "reverse"), "w")
    init_write_wig(mil_rw_handle, library_name, "reverse")
    minimum_factor = minimum_library_total / library_read_total
    per_million_factor = 1000000 / library_read_total
    for (seqid, mappings) in seqid_mapping:
        write_wig(raw_fw_handle, seqid, mappings["forward"])
        write_wig(raw_rw_handle, seqid, mappings["reverse"])
        write_wig(min_fw_handle, seqid, mappings["forward"], factor=minimum_factor)
        write_wig(min_rw_handle, seqid, mappings["reverse"], factor=minimum_factor)
        write_wig(mil_fw_handle, seqid, mappings["forward"], factor=per_million_factor)
        write_wig(mil_rw_handle, seqid, mappings["reverse"], factor=per_million_factor)
    raw_fw_handle.close()
    raw_rw_handle.close()
    min_fw_handle.close()
    min_rw_handle.close()
    mil_fw_handle.close()
    mil_rw_handle.close()


def get_read_count_dict(input_read_filepath, library_name):
    """Return this library's total and the smallest total in the experiment."""
    _, library_totals = misc.read_mapped_read_summary(input_read_filepath)
    if library_name not in library_totals:
        available = ", ".join(sorted(library_totals))
        raise ValueError(
            f"{input_read_filepath}: no mapped-read rows for library "
            f"{library_name!r}; available libraries: {available}"
        )
    return library_totals[library_name], min(library_totals.values())
