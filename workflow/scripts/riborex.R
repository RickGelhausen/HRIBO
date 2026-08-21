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
  make_option(c("-x", "--riborexdeseq_result_path"), type = "character", default = "NULL",
              help = "Path for writing deseq2 result file", metavar = "character"),
  make_option(c("-y", "--riborexedgeR_result_path"), type = "character", default = "NULL",
              help = "Path for writing edgeR result file", metavar = "character"),
  make_option(c("-z", "--riborexvoom_result_path"), type = "character", default = "NULL",
              help = "Path for writing voom result file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$rpf_in) || is.null(options$mrna_in) || is.null(options$condition_vector_in) || is.null(options$riborexdeseq_result_path) || is.null(options$riborexedgeR_result_path) || is.null(options$riborexvoom_result_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -m, -c, -x, -y, -z), see --help \n", call.=FALSE)
}

library(riborex)
# read the tsv file, and convert to a data frame with the first column as rownames
RNA <- read.table(options$mrna_in, sep="\t", header=TRUE, row.names=1)
RIBO <- read.table(options$rpf_in, sep="\t", header=TRUE, row.names=1)

# read the conditions vector text file and split it into a list
first_line <- readLines(options$condition_vector_in, n=1)
contrastconditionsvector <- unlist(strsplit(first_line, split = ","))

# run riborex analysis
results.deseq2 <- riborex(RNA, RIBO, contrastconditionsvector, contrastconditionsvector)

summary(results.deseq2)
# write results into file
write.csv(results.deseq2, options$riborexdeseq_result_path, quote = F)

#plot results
