import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.4.2")

ADAPTERS=config["adapter"]
INDEXPATH=config["genomeindexpath"]
CODONS=config["alternativestartcodons"]
TISHMODE=config["tishmode"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")
report: "report/workflow.rst"
def getContrast(wildcards):
  conditions=samples["condition"].unique()
  print(conditions)
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
  elements = [("figures/fc_" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + ".jpg") for element in flat_contrasts]
  return elements

def getContrastRiborex(wildcards):
  conditions=samples["condition"].unique()
  contrastsTupleList=list((iter.combinations(conditions,2)))
  contrasts=[[('-'.join(str(i) for i in x))] for x in contrastsTupleList]
  flat_contrasts= [item for sublist in contrasts for item in sublist]
  elements = [("riborex/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.csv") for element in flat_contrasts]
  return elements

if TISHMODE == "TISONLY":
   rule all:
      input:
         expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "TIS", "condition"]),
         expand("tracks/{condition}.ribotish.gff", zip, condition=samples["condition"]),
         expand("coverage/{method}-{condition}-{replicate}.bed", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         expand("tracks/{method}-{condition}-{replicate}.fwd.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         expand("tracks/{method}-{condition}-{replicate}.rev.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         unpack(getContrast),
         unpack(getContrastXtail),
         unpack(getContrastRiborex),
         "qc/multi/multiqc_report.html",
         expand("tracks/{condition}.merged.gff", zip, condition=samples["condition"]),
         "xtail/newAnnotation.gff",
         "figures/heatmap_SpearmanCorr_readCounts.pdf"

else:
   rule all:
      input:
         expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
         expand("tracks/{condition}.ribotish.gff", zip, condition=samples["condition"]),
         expand("tracks/{method}-{condition}-{replicate}.fwd.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         expand("tracks/{method}-{condition}-{replicate}.rev.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         expand("tracks/{condition}.reparation.gff", zip, condition=samples["condition"]),
         expand("figures/{condition}-{replicate}-qual.jpg", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
         expand("coverage/{method}-{condition}-{replicate}.bed", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
         unpack(getContrast),
         unpack(getContrastXtail),
         unpack(getContrastRiborex),
         "qc/multi/multiqc_report.html",
         expand("figures/{condition}-{replicate}_metagene.jpg", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples.loc[samples["method"] == "RIBO", "replicate"]),
         expand("tracks/{condition}.merged.gff", zip, condition=samples["condition"]),
         "xtail/newAnnotation.gff",
         "figures/heatmap_SpearmanCorr_readCounts.pdf"

onsuccess:
    print("Done, no error")

#Preprocessing
include: "rules/preprocessing.smk"
#Adaper removal and quality control
include: "rules/trimming.smk"
#removal of reads mapping to ribosomal rna genes
include: "rules/rrnafiltering.smk"
#mapping
include: "rules/mapping.smk"
#Visualization
include: "rules/visualization.smk"
#ribotish
include: "rules/ribotishauxiliary.smk"
if TISHMODE == "TISONLY":
   # merge
   include: "rules/mergetis.smk"
   include: "rules/ribotishtis.smk"
elif TISHMODE == "RIBOONLY":
   #merging
   include: "rules/merge.smk"
   #ribotish
   include: "rules/ribotish.smk"
   #include: "rules/ribotishauxiliary.smk"
else:
   #merging
   include: "rules/merge.smk"
   #ribotish
   include: "rules/ribotishall.smk"
   #include: "rules/ribotishauxiliary.smk"
#reparation
include: "rules/reparation.smk"
#xtail
include: "rules/xtail.smk"
#xtail
#include: "rules/xtailclassic.smk"
#multiqc
include: "rules/qc.smk"
