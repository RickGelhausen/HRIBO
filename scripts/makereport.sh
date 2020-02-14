#!/bin/bash
idx=$1

# timestamp
timestamp=$(date +"%d-%m-%y")
echo "${timestamp}"

report="${idx}report_HRIBO1.3.2_${timestamp}"

mkdir -p $report
mkdir -p "${report}/genome-browser/features/"
mkdir -p "${report}/supplementary/metagene/centered/min"
mkdir -p "${report}/supplementary/metagene/centered/mil"
mkdir -p "${report}/supplementary/metagene/fiveprime/min"
mkdir -p "${report}/supplementary/metagene/fiveprime/mil"
mkdir -p "${report}/supplementary/metagene/threeprime/min"
mkdir -p "${report}/supplementary/metagene/threeprime/mil"
mkdir -p "${report}/quality-control"
mkdir -p "${report}/ORF-predictions"
mkdir -p "${report}/genome-browser/coverage/global/min"
mkdir -p "${report}/genome-browser/coverage/global/mil"
mkdir -p "${report}/differential-expression/xtail/"
mkdir -p "${report}/differential-expression/riborex/"

cp -r HRIBO/manual/Manual_HRIBO.pdf "${report}/manual.pdf"
cp auxiliary/samples.xlsx "${report}/samples.xlsx"

cp -r "tracks/combined_annotated.gff" "${report}/ORF-predictions/predictions_reparation.gff"
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
cp -r centeredtracks/min/*.bw "${report}/supplementary/metagene/centered/min"
cp -r centeredtracks/mil/*.bw "${report}/supplementary/metagene/centered/mil"
cp -r fiveprimetracks/min/*.bw "${report}/supplementary/metagene/fiveprime/min"
cp -r fiveprimetracks/mil/*.bw "${report}/supplementary/metagene/fiveprime/mil"
cp -r threeprimetracks/min/*.bw "${report}/supplementary/metagene/threeprime/min"
cp -r threeprimetracks/mil/*.bw "${report}/supplementary/metagene/threeprime/mil"

cp -r xtail/*_sorted.csv "${report}/differential-expression/xtail/"
cp -r xtail/*_significant.csv "${report}/differential-expression/xtail/"
cp -r xtail/*.pdf "${report}/differential-expression/xtail/"
cp -r riborex/*_sorted.csv "${report}/differential-expression/riborex/"
cp -r riborex/*_significant.csv "${report}/differential-expression/riborex/"

cp -r metageneprofiling "${report}"

zip -r "${report}.zip" $report
