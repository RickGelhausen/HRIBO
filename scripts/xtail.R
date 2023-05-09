#!/usr/bin/env Rscript

library(optparse)
library(plyr)

option_list = list(
  make_option(c("-r", "--rpf_in"), type = "character", default = NULL,
              help = "Path to RPF read table", metavar = "character"),
  make_option(c("-m", "--mrna_in"), type = "character", default = NULL,
              help = "Path to mRNA read table", metavar = "character"),
  make_option(c("-c", "--condition_vector_in"), type = "character", default = NULL,
              help = "Contrast, pair of conditions ", metavar = "character"),
  make_option(c("-x", "--xtail_result_path"), type = "character", default = "NULL",
              help = "Path for writing xtail result file", metavar = "character"),
  make_option(c("-f", "--xtail_fcplot_path"), type = "character", default = "NULL",
              help = "Path for writing xtail fc plot file", metavar = "character"),
  make_option(c("-p", "--xtail_rplot_path"), type = "character", default = "NULL",
              help = "Path for writing xtail rplot file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$rpf_in) || is.null(options$mrna_in) || is.null(options$condition_vector_in) || is.null(options$xtail_result_path) || is.null(options$xtail_fcplot_path) || is.null(options$xtail_rplot_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -m, -c, -x -f -p), see --help \n", call.=FALSE)
}

library(xtail)

# read the tsv file, and convert to a data frame with the first column as rownames
RNA <- read.table(options$mrna_in, sep="\t", header=TRUE, row.names=1)
RIBO <- read.table(options$rpf_in, sep="\t", header=TRUE, row.names=1)

# read the conditions vector text file and split it into a list
first_line <- readLines(options$condition_vector_in, n=1)
contrastconditionsvector <- unlist(strsplit(first_line, split = ","))

# run xtail analysis
test.results <- xtail(RNA, RIBO, contrastconditionsvector)

# turn results into table
test.tab <- resultsTable(test.results, log2FCs = TRUE)
# write results into file
write.csv(test.tab, options$xtail_result_path, quote = F)

#plot results
pdf(file=options$xtail_fcplot_path, paper = "a4r", height = 10, width = 13)
plotFCs(test.results)
dev.off()
pdf(file=options$xtail_rplot_path, paper = "a4r", height = 10, width = 13)
plotRs(test.results)
dev.off()
