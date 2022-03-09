#!/bin/bash
idx=$1

# timestamp
timestamp=$(date +"%y-%m-%d")
echo "${timestamp}"

report="${idx}report_HRIBO1.5.1_${timestamp}"

mkdir -p $report
mkdir -p "${report}/genome-browser/features/"
mkdir -p "${report}/genome-browser/coverage/global/min"
mkdir -p "${report}/genome-browser/coverage/global/mil"
mkdir -p "${report}/genome-browser/coverage/centered/min"
mkdir -p "${report}/genome-browser/coverage/centered/mil"
mkdir -p "${report}/genome-browser/coverage/fiveprime/min"
mkdir -p "${report}/genome-browser/coverage/fiveprime/mil"
mkdir -p "${report}/genome-browser/coverage/threeprime/min"
mkdir -p "${report}/genome-browser/coverage/threeprime/mil"
mkdir -p "${report}/quality-control"
mkdir -p "${report}/ORF-predictions"
mkdir -p "${report}/differential-expression/xtail/"
mkdir -p "${report}/differential-expression/riborex/"

cp -r HRIBO/manual/Manual_HRIBO.pdf "${report}/manual.pdf"
cp auxiliary/samples.xlsx "${report}/samples.xlsx"

cp -r "tracks/reparation_annotated.gff" "${report}/ORF-predictions/predictions_reparation.gff"
cp -r "tracks/deepribo_merged.gff" "${report}/ORF-predictions/predictions_deepribo.gff"
cp -r auxiliary/predictions_deepribo.xlsx "${report}/ORF-predictions/predictions_deepribo.xlsx"
cp -r auxiliary/predictions_reparation.xlsx "${report}/ORF-predictions/predictions_reparation.xlsx"

cp -r tracks/potential*.gff "${report}/genome-browser/features/"
cp -r annotation/annotation.gff "${report}/genome-browser/annotation.gff"
cp -r genomes/genome.fa "${report}/genome-browser/genome.fa"

cp -r figures/heatmap_SpearmanCorr_readCounts.pdf "${report}/quality-control"
cp qc/multi/multiqc_report.html "${report}/quality-control/multiqc_report.html"
cp auxiliary/total_read_counts.xlsx "${report}/quality-control/total_read_counts.xlsx"
cp auxiliary/unique_read_counts.xlsx "${report}/quality-control/unique_read_counts.xlsx"
cp auxiliary/annotation_total.xlsx "${report}/quality-control/annotation_total.xlsx"
cp auxiliary/annotation_unique.xlsx "${report}/quality-control/annotation_unique.xlsx"

cp -r globaltracks/min/*.bw "${report}/genome-browser/coverage/global/min"
cp -r globaltracks/mil/*.bw "${report}/genome-browser/coverage/global/mil"
cp -r centeredtracks/min/*.bw "${report}/genome-browser/coverage/centered/min"
cp -r centeredtracks/mil/*.bw "${report}/genome-browser/coverage/centered/mil"
cp -r fiveprimetracks/min/*.bw "${report}/genome-browser/coverage/fiveprime/min"
cp -r fiveprimetracks/mil/*.bw "${report}/genome-browser/coverage/fiveprime/mil"
cp -r threeprimetracks/min/*.bw "${report}/genome-browser/coverage/threeprime/min"
cp -r threeprimetracks/mil/*.bw "${report}/genome-browser/coverage/threeprime/mil"

cp -r xtail/*_sorted.xlsx "${report}/differential-expression/xtail/"
cp -r xtail/*_significant.xlsx "${report}/differential-expression/xtail/"
cp -r xtail/*.pdf "${report}/differential-expression/xtail/"

cp -r auxiliary/overview.xlsx "${report}"
cp -r auxiliary/overview.gff "${report}/genome-browser/"
cp -r metageneprofiling "${report}"

zip -r "${report}.zip" $report
