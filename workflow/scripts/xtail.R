#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-r", "--rpf_in"), type = "character", default = NULL,
              help = "Path to RPF read table", metavar = "character"),
  make_option(c("-m", "--mrna_in"), type = "character", default = NULL,
              help = "Path to mRNA read table", metavar = "character"),
  make_option(c("-c", "--condition_vector_in"), type = "character", default = NULL,
              help = "Contrast, pair of conditions ", metavar = "character"),
  make_option(c("-x", "--xtail_result_path"), type = "character", default = NULL,
              help = "Path for writing xtail result file", metavar = "character"),
  make_option(c("-f", "--xtail_fcplot_path"), type = "character", default = NULL,
              help = "Path for writing xtail fc plot file", metavar = "character"),
  make_option(c("-p", "--xtail_rplot_path"), type = "character", default = NULL,
              help = "Path for writing xtail rplot file", metavar = "character"),
  make_option(c("--threads"), type = "integer", default = 1L,
              help = "Number of worker processes used by xTail [default %default]",
              metavar = "integer"),
  make_option(c("--bins"), type = "integer", default = 10000L,
              help = "Probability-density bins used by xTail [default %default]",
              metavar = "integer"),
  make_option(c("--min_mean_count"), type = "double", default = 1,
              help = "Minimum mean RNA/RPF count retained by xTail [default %default]",
              metavar = "number")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$rpf_in) || is.null(options$mrna_in) || is.null(options$condition_vector_in) || is.null(options$xtail_result_path) || is.null(options$xtail_fcplot_path) || is.null(options$xtail_rplot_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -m, -c, -x -f -p), see --help \n", call.=FALSE)
}

if (is.na(options$threads) || options$threads < 1L) {
  stop("--threads must be a positive integer", call.=FALSE)
}

if (is.na(options$bins) || options$bins < 1L) {
  stop("--bins must be a positive integer", call.=FALSE)
}

if (is.na(options$min_mean_count) || !is.finite(options$min_mean_count) || options$min_mean_count < 1) {
  stop("--min_mean_count must be a number greater than or equal to 1", call.=FALSE)
}

library(xtail)

# read the tsv file, and convert to a data frame with the first column as rownames
RNA <- read.table(options$mrna_in, sep="\t", header=TRUE, row.names=1)
RIBO <- read.table(options$rpf_in, sep="\t", header=TRUE, row.names=1)

# read the conditions vector text file and split it into a list
first_line <- readLines(options$condition_vector_in, n=1)
contrastconditionsvector <- unlist(strsplit(first_line, split = ","))

# run xtail analysis
test.results <- xtail(
  RNA,
  RIBO,
  contrastconditionsvector,
  threads = options$threads,
  bins = options$bins,
  minMeanCount = options$min_mean_count
)

# turn results into table
test.tab <- resultsTable(test.results, log2FCs = TRUE)
# xTail 1.2.0 accidentally retains one condition-specific `*_log2TE` column
# when log2Rs is false. Keep HRIBO's stable public result schema and fail clearly
# if a future xTail release removes a required statistic.
result.columns <- c(
  "mRNA_log2FC",
  "RPF_log2FC",
  "log2FC_TE_v1",
  "pvalue_v1",
  "log2FC_TE_v2",
  "pvalue_v2",
  "log2FC_TE_final",
  "pvalue_final",
  "pvalue.adjust"
)
missing.columns <- setdiff(result.columns, colnames(test.tab))
if (length(missing.columns) > 0L) {
  stop(
    paste("xTail result is missing required columns:", paste(missing.columns, collapse = ", ")),
    call. = FALSE
  )
}
test.tab <- test.tab[, result.columns, drop = FALSE]
# write results into file
write.csv(test.tab, options$xtail_result_path, quote = F)

#plot results
pdf(file=options$xtail_fcplot_path, paper = "a4r", height = 10, width = 13)
plotFCs(test.results)
dev.off()
pdf(file=options$xtail_rplot_path, paper = "a4r", height = 10, width = 13)
plotRs(test.results)
dev.off()
