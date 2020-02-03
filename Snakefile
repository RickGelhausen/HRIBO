import os
import re
import pandas as pd
import itertools as iter
from snakemake.utils import validate, min_version
min_version("5.5.1")

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
  elements = [("xtail/" + ((element.replace("[", '')).replace("]", '')).replace("'", '') + "_significant.csv") for element in flat_contrasts]
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
          #expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
          #expand("tracks/{condition}.ribotish.gff", zip, condition=samples["condition"]),
          #expand("tracks/{method}-{condition}-{replicate}.fwd.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          #expand("tracks/{method}-{condition}-{replicate}.rev.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/raw/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/norm/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/min/{method}-{condition}-{replicate}.min.forward.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/min/{method}-{condition}-{replicate}.min.reverse.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/min/{method}-{condition}-{replicate}.min.forward.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/min/{method}-{condition}-{replicate}.min.reverse.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/min/{method}-{condition}-{replicate}.min.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/min/{method}-{condition}-{replicate}.min.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/min/{method}-{condition}-{replicate}.min.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("coverage/{method}-{condition}-{replicate}.bed", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          "qc/multi/multiqc_report.html",
          "xtail/newAnnotation.gff",
          "figures/heatmap_SpearmanCorr_readCounts.pdf",
          "tracks/potentialStopCodons.gff",
          "tracks/potentialStartCodons.gff",
          "tracks/potentialAlternativeStartCodons.gff",
          "tracks/potentialRibosomeBindingSite.gff",
          "offsets/metagene_start_rois.txt",
          "offsets/metagene_stop_rois.txt",
          #expand("offsets/start/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          #expand("offsets/stop/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          #expand("tracks/color/{method}-{condition}-{replicate}.fwd.bedgraph.gz", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          "tracks/color/potentialStartCodons.gff",
          "auxiliary/annotation_total.xlsx",
          "auxiliary/annotation_unique.xlsx",
          "auxiliary/total_read_counts.xlsx",
          "auxiliary/unique_read_counts.xlsx",
          "auxiliary/samples.xlsx",
          "auxiliary/summary.xlsx",
          unpack(getContrast),
          unpack(getContrastXtail),
          unpack(getContrastRiborex)
else:
   rule all:
      input:
          #expand("ribotish/{condition}-newORFs.tsv_all.txt", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"]),
          #expand("tracks/{condition}.ribotish.gff", zip, condition=samples["condition"]),
          #expand("tracks/{method}-{condition}-{replicate}.fwd.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          #expand("tracks/{method}-{condition}-{replicate}.rev.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/raw/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("metageneprofiling/norm/{method}-{condition}-{replicate}", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/min/{method}-{condition}-{replicate}.min.forward.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("totalmappedtracks/min/{method}-{condition}-{replicate}.min.reverse.totalmapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.forward.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/raw/{method}-{condition}-{replicate}.raw.reverse.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.forward.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/mil/{method}-{condition}-{replicate}.mil.reverse.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/min/{method}-{condition}-{replicate}.min.forward.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("uniquemappedtracks/min/{method}-{condition}-{replicate}.min.reverse.uniquemapped.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/raw/{method}-{condition}-{replicate}.raw.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/mil/{method}-{condition}-{replicate}.mil.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/min/{method}-{condition}-{replicate}.min.forward.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("globaltracks/min/{method}-{condition}-{replicate}.min.reverse.global.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/raw/{method}-{condition}-{replicate}.raw.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/mil/{method}-{condition}-{replicate}.mil.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/min/{method}-{condition}-{replicate}.min.forward.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("centeredtracks/min/{method}-{condition}-{replicate}.min.reverse.centered.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.forward.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("fiveprimetracks/min/{method}-{condition}-{replicate}.min.reverse.fiveprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/raw/{method}-{condition}-{replicate}.raw.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/mil/{method}-{condition}-{replicate}.mil.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.forward.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("threeprimetracks/min/{method}-{condition}-{replicate}.min.reverse.threeprime.bw", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("coverage/{method}-{condition}-{replicate}.bed", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          expand("tracks/{condition}.reparation.gff", zip, condition=samples["condition"]),
          expand("maplink/{method}-{condition}-{replicate}.bam", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          "qc/multi/multiqc_report.html",
          "xtail/newAnnotation.gff",
          #expand("figures/{condition}-{replicate}_metagene.jpg", zip, condition=samples.loc[samples["method"] == "RIBO", "condition"], replicate=samples["replicate"]),
          "figures/heatmap_SpearmanCorr_readCounts.pdf",
          "tracks/potentialStopCodons.gff",
          "tracks/potentialStartCodons.gff",
          "tracks/potentialAlternativeStartCodons.gff",
          "tracks/potentialRibosomeBindingSite.gff",
          "offsets/metagene_start_rois.txt",
          "offsets/metagene_stop_rois.txt",
          #expand("offsets/start/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          #expand("offsets/stop/{method}-{condition}-{replicate}_p_offsets.txt", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          #expand("tracks/color/{method}-{condition}-{replicate}.fwd.bedgraph.gz", zip, method=samples["method"], condition=samples["condition"], replicate=samples["replicate"]),
          "tracks/color/potentialStartCodons.gff",
          "auxiliary/annotation_total.xlsx",
          "auxiliary/annotation_unique.xlsx",
          "auxiliary/total_read_counts.xlsx",
          "auxiliary/unique_read_counts.xlsx",
          "auxiliary/samples.xlsx",
          "auxiliary/summary.xlsx",
          unpack(getContrast),
          unpack(getContrastXtail),
          unpack(getContrastRiborex)

onsuccess:
    print("Done, no error")

#Preprocessing
include: "rules/preprocessing.smk"
#Adaper removal and quality control
include: "rules/trimming.smk"
#removal of reads mapping to ribosomal rna genes
include: "rules/rrnafiltering.smk"
#mapping
include: "rules/mappingsingleend.smk"
include: "rules/mappingauxiliary.smk"
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
else:
   #merging
   include: "rules/merge.smk"
   #ribotish
   include: "rules/ribotishall.smk"
#reparation
include: "rules/reparation.smk"
#xtail
include: "rules/xtail.smk"
#metagene
include: "rules/metageneprofiling.smk"
include: "rules/auxiliary.smk"
#multiqc
include: "rules/qcauxiliary.smk"
include: "rules/qcsingleend.smk"
#report
#include: "rules/report.smk"
