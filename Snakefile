import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
from validate import validate_config

min_version("7.18.2")

validate_config(config)

ADAPTERS=config["biologySettings"]["adapter"]
CODONS=config["biologySettings"]["alternativeStartCodons"]
DIFFEXPRESS=config["differentialExpressionSettings"]["differentialExpression"]
CONTRASTS=config["differentialExpressionSettings"]["contrasts"]
DEEPRIBO=config["predictionSettings"]["deepribo"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["biologySettings"]["samples"], dtype=str, sep="\t")
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


def get_wigfiles():
    method=samples["method"]
    condition=samples["condition"]
    replicate=samples["replicate"]
    wilds = zip(method, condition, replicate)

    bigwigs = [["global", "centered", "fiveprime", "threeprime"], ["raw", "mil", "min"], ["forward", "reverse"], list(wilds)]
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

hribo_output = []
hribo_output.extend(expand("metageneprofiling/TIS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]))
hribo_output.extend(expand("metageneprofiling/TIS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_start["method"], condition=samples_meta_start["condition"], replicate=samples_meta_start["replicate"]))
hribo_output.extend(expand("metageneprofiling/TTS/raw/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]))
hribo_output.extend(expand("metageneprofiling/TTS/norm/{method}-{condition}-{replicate}", zip, method=samples_meta_stop["method"], condition=samples_meta_stop["condition"], replicate=samples_meta_stop["replicate"]))
hribo_output.append("qc/multi/multiqc_report.html")
hribo_output.append("tracks/potentialStopCodons.gff")
hribo_output.append("tracks/potentialStartCodons.gff")
hribo_output.append("tracks/potentialRibosomeBindingSite.gff")
hribo_output.append("auxiliary/annotation_total.xlsx")
hribo_output.append("auxiliary/annotation_unique.xlsx")
hribo_output.append("auxiliary/total_read_counts.xlsx")
hribo_output.append("auxiliary/unique_read_counts.xlsx")
hribo_output.append("auxiliary/samples.xlsx")
hribo_output.append("figures/heatmap_SpearmanCorr_readCounts.pdf")

hribo_output.extend(get_wigfiles())

if hasRIBO:
    hribo_output.append("auxiliary/predictions_reparation.xlsx")

    conditions=sorted(samples["condition"].unique(), key=lambda s: s.lower())

    if DIFFEXPRESS.lower() == "on":
        if CONTRASTS != "":
            CONTRASTS = CONTRASTS.split(",")
        else:
            CONTRASTS=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(sorted(samples["condition"].unique(), key=lambda s: s.lower()),2)))]  for item in sublist]

        hribo_output.extend([("contrasts/"+((element.replace("[", "")).replace("]", "")).replace("'", "")) for element in CONTRASTS])
        hribo_output.extend([("xtail/" + ((element.replace("[", "")).replace("]", "")).replace("'", "") + "_significant.xlsx") for element in CONTRASTS])
        hribo_output.extend([("riborex/" + ((element.replace("[", "")).replace("]", "")).replace("'", "") + "_significant.xlsx") for element in CONTRASTS])
        hribo_output.extend([("deltate/" + ((element.replace("[", "")).replace("]", "")).replace("'", "") + "_significant.xlsx") for element in CONTRASTS])

    if DEEPRIBO.lower() == "on":
        hribo_output.append("auxiliary/predictions_deepribo.xlsx")


    if DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "on":
        hribo_output.extend(rules.createOverviewTableAll.output)
    elif DIFFEXPRESS.lower() == "off" and DEEPRIBO.lower() == "on":
        hribo_output.extend(rules.createOverviewTablePredictions.output)
    elif DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "off":
        hribo_output.extend(rules.createOverviewTableDiffExpr.output)
    elif DIFFEXPRESS.lower() == "off" and DEEPRIBO.lower() == "off":
        hribo_output.extend(rules.createOverviewTableReparation.output)

rule all:
    input:
        hribo_output

onsuccess:
    print("Done, no error")
