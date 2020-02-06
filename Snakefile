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

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")
report: "report/workflow.rst"
def getContrast(wildcards):
  conditions=samples["condition"].unique()
  contrastsTupleList=list((iter.combinations(conditions,2)))
  contrasts=[[('-'.join(str(i) for i in x))] for x in contrastsTupleList]
  flat_contrasts= [item for sublist in contrasts for item in sublist]
  elements = [("contrasts/"+((element.replace("[", '')).replace("]", '')).replace("'", '')) for element in flat_contrasts]
  return elements

def getContrastXtail(wildcards):
  conditions=samples["condition"].unique()
  contrastsTupleList=list((iter.combinations(conditions,2)))
  contrasts=[[('-'.join(str(i) for i in x))] for x in contrastsTupleList]
  flat_contrasts= [item for sublist in contrasts for item in sublist]
  elements = [("xtail/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.csv") for element in flat_contrasts]
  return elements

def getContrastRiborex(wildcards):
  conditions=samples["condition"].unique()
  contrastsTupleList=list((iter.combinations(conditions,2)))
  contrasts=[[('-'.join(str(i) for i in x))] for x in contrastsTupleList]
  flat_contrasts= [item for sublist in contrasts for item in sublist]
  elements = [("riborex/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.csv") for element in flat_contrasts]
  return elements

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

if DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "on":
   rule all:
      input:
          expand("metageneprofiling/raw/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/norm/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
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
          "auxiliary/summary.xlsx",
          "figures/heatmap_SpearmanCorr_readCounts.pdf",
          "auxiliary/predictions_deepribo.xlsx"
          unpack(getContrast),
          unpack(getContrastXtail),
          unpack(getContrastRiborex)

elif DIFFEXPRESS.lower() == "off" and DEEPRIBO.lower() == "on":
   rule all:
      input:
          expand("metageneprofiling/raw/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/norm/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
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
          "auxiliary/summary.xlsx",
          "figures/heatmap_SpearmanCorr_readCounts.pdf",
          "auxiliary/predictions_deepribo.xlsx"

elif DIFFEXPRESS.lower() == "on" and DEEPRIBO.lower() == "off":
   rule all:
      input:
          expand("metageneprofiling/raw/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/norm/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
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
          "auxiliary/summary.xlsx",
          "figures/heatmap_SpearmanCorr_readCounts.pdf",
          unpack(getContrast),
          unpack(getContrastXtail),
          unpack(getContrastRiborex)

else:
   rule all:
      input:
          expand("metageneprofiling/raw/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/norm/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
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
          "auxiliary/summary.xlsx",
          "figures/heatmap_SpearmanCorr_readCounts.pdf",

onsuccess:
    print("Done, no error")

# Preprocessing
include: "rules/preprocessing.smk"
# Adaper removal and quality control
include: "rules/trimming.smk"
# removal of reads mapping to ribosomal rna genes
include: "rules/rrnafiltering.smk"
# mapping
include: "rules/mappingsingleend.smk"
include: "rules/mappingauxiliary.smk"
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

if DIFFEXPRESS.lower() == "on":
    # xtail
    include: "rules/xtail.smk"

if DEEPRIBO.lower() == "on":
    #deepribo
    include: "rules/deepribo.smk"
