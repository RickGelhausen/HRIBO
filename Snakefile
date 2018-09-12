import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.2.4")

ADAPTERS=config["adapter"]
INDEXPATH=config["genomeindexpath"]

onstart:
   if not os.path.exists("logs"):
     os.makedirs("logs")

samples = pd.read_table(config["samples"], dtype=str).set_index(["method", "condition", "replicate"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])
validate(samples, schema="schemas/samples.schema.yaml")
report: "report/workflow.rst"
def getContrast(wildcards):
  conditionsRaw=samples.loc[(wildcards.method, wildcards.condition, wildcards.replicate), ["condition"]].dropna()
  conditions=conditionsRaw.values.tolist()
  contrasts=iter.combinations(conditions, 2)
  print(contrasts)
  return contrasts

rule all:
   input:
       expand("fastqc/raw/{sample.method}-{sample.condition}-{sample.replicate}-raw.html", sample=samples.itertuples()),
       expand("fastqc/trimmed/{sample.method}-{sample.condition}-{sample.replicate}-trimmed.html", sample=samples.itertuples()),
       expand("fastqc/norRNA/{sample.method}-{sample.condition}-{sample.replicate}-norRNA.html", sample=samples.itertuples()),
       expand("ribotish/{sample.condition}-{sample.replicate}-newORFs.tsv_all.txt", sample=samples.itertuples()),
       expand("tracks/{sample.method}-{sample.condition}-{sample.replicate}.bw", sample=samples.itertuples()),
       expand("xtailclassic/{sample.method}-{sample.condition}-{sample.replicate}.csv", sample=samples.itertuples()),
       expand("xtail/{contrast}.csv", contrast=getContrast, sample=samples.itertuples()),
       "xtail/sfactors.csv",
       "xtail/xtail.csv",

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
include: "rules/ribotish.smk"
#xtail
include: "rules/xtail.smk"
#xtail
include: "rules/xtailclassic.smk"

