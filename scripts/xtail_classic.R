#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-r", "--read_counts_csv_path"), type = "character", default = NULL,
              help = "Path to read counts table", metavar = "character"),
  make_option(c("-t", "--sample_file_path"), type = "character", default = NULL,
              help = "Path to sample.tsv", metavar = "character"),
  make_option(c("-x", "--xtail_result_path"), type = "character", default = "NULL",
              help = "Path for writing xtail result file", metavar = "character"),
  make_option(c("-f", "--xtail_fcplot_path"), type = "character", default = "NULL",
              help = "Path for writing xtail fc plot file", metavar = "character"),
  make_option(c("-p", "--xtail_rplot_path"), type = "character", default = "NULL",
              help = "Path for writing xtail rplot file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$normalized_read_counts_csv_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -t, -x), see --help \n", call.=FALSE)
}

library(xtail)

# read table with mormalized read counts
counts <- read.csv(options$read_counts_csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path ,header=TRUE, sep="\t", stringsAsFactors=FALSE)

# split data frame into RIBO and RNA
RIBO <- counts[, sampleSheet$method == "RIBO"]
RNA <- counts[, sampleSheet$method == "RNA"]

#create condition vextor
condition <- sampleSheet$condition[which(sampleSheet$method == "RIBO")]

# run xtail analysis
test_results <- xtail(RNA, RIBO, condition)

# turn results into table
test_tab <- resultsTable(test_results, log2FCs = TRUE)

# write results into file
write.csv(test_tab, options$xtail_result_path, quote = F)

#plot results
pdf(file=options$xtail_fcplot_path, paper = "a4r", height = 10, width = 13)
plotFCs(test_results)
dev.off()
pdf(file=options$xtail_rplot_path, paper = "a4r", height = 10, width = 13)
plotRs(test_results)
dev.off()

