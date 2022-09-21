import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.5.1")

ADAPTERS=config["adapter"]
CODONS=config["alternativestartcodons"]
DIFFEXPRESS=config["differentialexpression"]
DEEPRIBO=config["deepribo"]
CONTRASTS=config["contrasts"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["samples"], dtype=str, sep="\t")
samples.sort_values(by=["method", "condition", "replicate"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)
samples.set_index(["method", "condition", "replicate"], drop=False, inplace=True)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")

samples_meta_start = samples.loc[(samples["method"] == "RIBO") | (samples["method"] == "TIS")]
samples_meta_stop = samples.loc[(samples["method"] == "TTS")]

if DIFFEXPRESS.lower() == "on" and len(samples["condition"].unique()) <= 1:
    sys.exit("Differential Expression requested, but only one condition given.\n\
            Please ensure, that you either provide multiple condtions or turn off differential expression in the config.yaml.")

hasRIBO=True
if "RIBO" not in samples["method"].unique():
    hasRIBO=False
    print("No Ribo-seq libraries were detected. No prediction tools for this setup are currently implemented. If you have pure Ribo-seq libraries, please use the method tag RIBO. Continuing...")

report: "report/workflow.rst"

conditions=sorted(samples["condition"].unique(), key=lambda s: s.lower())

if DIFFEXPRESS.lower() == "on":
    if CONTRASTS != "":
        CONTRASTS = CONTRASTS.split(",")
    else:
        CONTRASTS=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(sorted(samples["condition"].unique(), key=lambda s: s.lower()),2)))]  for item in sublist]


    def getContrast(CONTRASTS):
        return [("contrasts/"+((element.replace("[", '')).replace("]", '')).replace("'", '')) for element in CONTRASTS]

    def getContrastXtail(CONTRASTS):
        return [("xtail/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.xlsx") for element in CONTRASTS]

    def getContrastRiborex(CONTRASTS):
        return [("riborex/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.xlsx") for element in CONTRASTS]

    def getContrastDeltaTE(CONTRASTS):
        return [("deltate/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.xlsx") for element in CONTRASTS]

def get_wigfiles(wildcards):
    method=samples["method"]
    condition=samples["condition"]
    replicate=samples["replicate"]
    wilds = zip(method, condition, replicate)

    bigwigs = [["totalmapped", "uniquemapped", "global", "centered", "fiveprime", "threeprime"], ["raw", "mil", "min"], ["forward", "reverse"], list(wilds)]
    bigwigs = list(iter.product(*bigwigs))

    wigfiles = []
    for bw in bigwigs:
        wigfiles.append("%stracks/%s/%s-%s-%s.%s.%s.%s.bw" %(bw[0], bw[1], bw[3][0], bw[3][1], bw[3][2], bw[1], bw[2], bw[0]))

    return wigfiles


# Preprocessing
include: "rules/preprocessing.smk"
# Adaper removal and quality control
include: "rules/trimming.smk"
# removal of reads mapping to ribosomal rna genes
include: "rules/rrnafiltering.smk"
# mapping
include: "rules/mapping.smk"
# Visualization
include: "rules/visualization.smk"
include: "rules/merge.smk"
# reparation
include: "rules/reparation.smk"
# metagene
include: "rules/metageneprofiling.smk"
include: "rules/auxiliary.smk"
# multiqc
include: "rules/qcauxiliary.smk"
include: "rules/qcsingleend.smk"
#readcounts
include: "rules/readcounting.smk"
if DIFFEXPRESS.lower() == "on":
    include: "rules/diffex_contrast.smk"
    include: "rules/diffex_xtail.smk"
    include: "rules/diffex_riborex.smk"
    include: "rules/diffex_deltate.smk"

if DEEPRIBO.lower() == "on":
    #deepribo
    include: "rules/deepribo.smk"
else:
    include: "rules/conditionals.smk"

if hasRIBO:
    if DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "on":
       rule all:
          input:
              expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              get_wigfiles,
              "qc/multi/multiqc_report.html",
              "tracks/potentialStopCodons.gff",
              "tracks/potentialStartCodons.gff",
              "tracks/potentialAlternativeStartCodons.gff",
              "tracks/potentialRibosomeBindingSite.gff",
              "auxiliary/annotation_total.xlsx",
              "auxiliary/annotation_unique.xlsx",
              "auxiliary/total_read_counts.xlsx",
              "auxiliary/unique_read_counts.xlsx",
              "auxiliary/samples.xlsx",
              "auxiliary/predictions_reparation.xlsx",
              "figures/heatmap_SpearmanCorr_readCounts.pdf",
              "auxiliary/predictions_deepribo.xlsx",
              rules.createOverviewTableAll.output,
              unpack(getContrast),
              unpack(getContrastXtail),
              unpack(getContrastRiborex),
              unpack(getContrastDeltaTE),
              "metageneprofiling/merged_offsets.json"


    elif DIFFEXPRESS.lower() == "off" and DEEPRIBO.lower() == "on":
       rule all:
          input:
              expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              get_wigfiles,
              "qc/multi/multiqc_report.html",
              "tracks/potentialStopCodons.gff",
              "tracks/potentialStartCodons.gff",
              "tracks/potentialAlternativeStartCodons.gff",
              "tracks/potentialRibosomeBindingSite.gff",
              "auxiliary/annotation_total.xlsx",
              "auxiliary/annotation_unique.xlsx",
              "auxiliary/total_read_counts.xlsx",
              "auxiliary/unique_read_counts.xlsx",
              "auxiliary/samples.xlsx",
              "auxiliary/predictions_reparation.xlsx",
              "figures/heatmap_SpearmanCorr_readCounts.pdf",
              "auxiliary/predictions_deepribo.xlsx",
              rules.createOverviewTablePredictions.output,
              "metageneprofiling/merged_offsets.json"

    elif DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "off":
       rule all:
          input:
              expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              get_wigfiles,
              "qc/multi/multiqc_report.html",
              "tracks/potentialStopCodons.gff",
              "tracks/potentialStartCodons.gff",
              "tracks/potentialAlternativeStartCodons.gff",
              "tracks/potentialRibosomeBindingSite.gff",
              "auxiliary/annotation_total.xlsx",
              "auxiliary/annotation_unique.xlsx",
              "auxiliary/total_read_counts.xlsx",
              "auxiliary/unique_read_counts.xlsx",
              "auxiliary/samples.xlsx",
              "auxiliary/predictions_reparation.xlsx",
              "figures/heatmap_SpearmanCorr_readCounts.pdf",
              rules.createOverviewTableDiffExpr.output,
              unpack(getContrast),
              unpack(getContrastXtail),
              unpack(getContrastRiborex),
              unpack(getContrastDeltaTE),
              "metageneprofiling/merged_offsets.json"

    else:
       rule all:
          input:
              expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              get_wigfiles,
              "qc/multi/multiqc_report.html",
              "tracks/potentialStopCodons.gff",
              "tracks/potentialStartCodons.gff",
              "tracks/potentialAlternativeStartCodons.gff",
              "tracks/potentialRibosomeBindingSite.gff",
              "auxiliary/annotation_total.xlsx",
              "auxiliary/annotation_unique.xlsx",
              "auxiliary/total_read_counts.xlsx",
              "auxiliary/unique_read_counts.xlsx",
              "auxiliary/samples.xlsx",
              "auxiliary/predictions_reparation.xlsx",
              "figures/heatmap_SpearmanCorr_readCounts.pdf",
              rules.createOverviewTableReparation.output,
              "metageneprofiling/merged_offsets.json"
else:
     if DIFFEXPRESS.lower() == "on":
       rule all:
          input:
              expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              get_wigfiles,
              "tracks/potentialStopCodons.gff",
              "tracks/potentialStartCodons.gff",
              "tracks/potentialAlternativeStartCodons.gff",
              "tracks/potentialRibosomeBindingSite.gff",
              "auxiliary/annotation_total.xlsx",
              "auxiliary/annotation_unique.xlsx",
              "auxiliary/total_read_counts.xlsx",
              "auxiliary/unique_read_counts.xlsx",
              "auxiliary/samples.xlsx",
              "figures/heatmap_SpearmanCorr_readCounts.pdf",
              unpack(getContrast),
              unpack(getContrastXtail),
              unpack(getContrastRiborex),
              unpack(getContrastDeltaTE),
              "metageneprofiling/merged_offsets.json"
     else:
       rule all:
          input:
              expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]),
              expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]),
              get_wigfiles,
              "tracks/potentialStopCodons.gff",
              "tracks/potentialStartCodons.gff",
              "tracks/potentialAlternativeStartCodons.gff",
              "tracks/potentialRibosomeBindingSite.gff",
              "auxiliary/annotation_total.xlsx",
              "auxiliary/annotation_unique.xlsx",
              "auxiliary/total_read_counts.xlsx",
              "auxiliary/unique_read_counts.xlsx",
              "auxiliary/samples.xlsx",
              "figures/heatmap_SpearmanCorr_readCounts.pdf",
              "metageneprofiling/merged_offsets.json"

onsuccess:
    print("Done, no error")
