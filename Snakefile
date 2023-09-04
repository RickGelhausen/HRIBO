import os
import sys
import pandas as pd
import itertools as iter
from snakemake.utils import min_version
from validate import validate_config, validate_sample_sheet

min_version("7.18.2")


onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["biologySettings"]["samples"], dtype=str, sep="\t")
validate_sample_sheet(samples)

samples.sort_values(by=["method", "condition", "replicate"], key=lambda col: col.str.lower() if col.dtype==str else col, inplace=True)

samples_metagene = samples.loc[(samples["method"] == "RIBO") | (samples["method"] == "TIS") | (samples["method"] == "TTS")]


conditions=sorted(samples["condition"].unique(), key=lambda s: s.lower())
validate_config(config, conditions)

ADAPTERS_S3=config["biologySettings"]["adapterS3"]
ADAPTERS_S5=config["biologySettings"]["adapterS5"]

ADAPTERS_P3R1=config["biologySettings"]["adapterP3R1"]
ADAPTERS_P5R1=config["biologySettings"]["adapterP5R1"]
ADAPTERS_P3R2=config["biologySettings"]["adapterP3R2"]
ADAPTERS_P5R2=config["biologySettings"]["adapterP5R2"]

CODONS=config["biologySettings"]["alternativeStartCodons"]
DIFFEXPRESS=config["differentialExpressionSettings"]["differentialExpression"]
CONTRASTS=config["differentialExpressionSettings"]["contrasts"]
DEEPRIBO=config["predictionSettings"]["deepribo"]

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
# # Visualization
# include: "rules/visualization.smk"
# include: "rules/merge.smk"
# # reparation
# include: "rules/reparation.smk"
# # metagene
# include: "rules/metageneprofiling.smk"
include: "rules/auxiliary.smk"
# multiqc
include: "rules/qcauxiliary.smk"
include: "rules/qc.smk"
# #readcounts
# include: "rules/readcounting.smk"
# if DIFFEXPRESS.lower() == "on":
#     include: "rules/diffex_contrast.smk"
#     include: "rules/diffex_xtail.smk"
#     include: "rules/diffex_riborex.smk"
#     include: "rules/diffex_deltate.smk"

# if DEEPRIBO.lower() == "on":
#     #deepribo
#     include: "rules/deepribo.smk"
# else:
#     include: "rules/conditionals.smk"

# include: "rules/pca.smk"

hribo_output = []
#hribo_output.extend(["maplink/RIBO-WT-1.bam", "maplink/RIBO-WT-2.bam", "maplink/TIS-WT-1.bam", "maplink/TIS-WT-2.bam"])
# hribo_output.extend(expand("metageneprofiling/{method}-{condition}-{replicate}", zip, method=samples_metagene["method"], condition=samples_metagene["condition"], replicate=samples_metagene["replicate"]))
# hribo_output.append("metageneprofiling/read_length_fractions.html")
hribo_output.append("qc/multi/multiqc_report.html")
# hribo_output.append("tracks/potentialStopCodons.gff")
# hribo_output.append("tracks/potentialStartCodons.gff")
# hribo_output.append("tracks/potentialRibosomeBindingSite.gff")
# hribo_output.append("tracks/potentialAlternativeStartCodons.gff")
# hribo_output.append("auxiliary/annotation_total.xlsx")
# hribo_output.append("auxiliary/annotation_unique.xlsx")
# hribo_output.append("auxiliary/total_read_counts.xlsx")
# hribo_output.append("auxiliary/unique_read_counts.xlsx")
# hribo_output.append("auxiliary/samples.xlsx")
# hribo_output.append("figures/heatmap_SpearmanCorr_readCounts.pdf")
# hribo_output.append("pca/PCA_3D.html")
# hribo_output.extend(get_wigfiles())

# if hasRIBO:
#     hribo_output.append("auxiliary/predictions_reparation.xlsx")

#     if DIFFEXPRESS.lower() == "on":
#         if CONTRASTS == []:
#             CONTRASTS=[item for sublist in [[('-'.join(str(i) for i in x))] for x in list((iter.combinations(sorted(samples["condition"].unique(), key=lambda s: s.lower()),2)))]  for item in sublist]

#         hribo_output.extend([("contrasts/" +((element.replace("[", "")).replace("]", "")).replace("'", "")) for element in CONTRASTS])
#         hribo_output.extend([("xtail/" + ((element.replace("[", "")).replace("]", "")).replace("'", "") + "_sorted.xlsx") for element in CONTRASTS])
#         hribo_output.extend([("riborex/" + ((element.replace("[", "")).replace("]", "")).replace("'", "") + "_sorted.xlsx") for element in CONTRASTS])
#         hribo_output.extend([("deltate/" + ((element.replace("[", "")).replace("]", "")).replace("'", "") + "_sorted.xlsx") for element in CONTRASTS])

#     if DEEPRIBO.lower() == "on":
#         hribo_output.append("auxiliary/predictions_deepribo.xlsx")


#     if DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "on":
#         hribo_output.extend(rules.createOverviewTableAll.output)
#     elif DIFFEXPRESS.lower() == "off" and DEEPRIBO.lower() == "on":
#         hribo_output.extend(rules.createOverviewTablePredictions.output)
#     elif DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "off":
#         hribo_output.extend(rules.createOverviewTableDiffExpr.output)
#     elif DIFFEXPRESS.lower() == "off" and DEEPRIBO.lower() == "off":
#         hribo_output.extend(rules.createOverviewTableReparation.output)

rule all:
    input:
        hribo_output

onsuccess:
    print("Done, no error")
