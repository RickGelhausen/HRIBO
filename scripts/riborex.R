#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-r", "--raw_read_counts_csv_path"), type = "character", default = NULL,
              help = "Path to normalized read counts table", metavar = "character"),
  make_option(c("-c", "--contrast"), type = "character", default = NULL,
              help = "Contrast, pair of conditions ", metavar = "character"),
  make_option(c("-t", "--sample_file_path"), type = "character", default = NULL,
              help = "Path to sample.tsv", metavar = "character"),
  make_option(c("-x", "--riborexdeseq_result_path"), type = "character", default = "NULL",
              help = "Path for writing deseq2 result file", metavar = "character"),
  make_option(c("-y", "--riborexedgeR_result_path"), type = "character", default = "NULL",
              help = "Path for writing edgeR result file", metavar = "character"),
  make_option(c("-z", "--riborexvoom_result_path"), type = "character", default = "NULL",
              help = "Path for writing voom result file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$raw_read_counts_csv_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -t, -x, -y, -z), see --help \n", call.=FALSE)
}

library(riborex)

# read table with raw read counts
counts <- read.csv(options$raw_read_counts_csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path ,header=TRUE, sep="\t", stringsAsFactors=FALSE)

#create condition vector
constraststring <- gsub("contrasts/", "", options$contrast)
#print(constraststring)
contrastconditions <- unlist(strsplit(constraststring,"-"))
cond1 <- contrastconditions[1]
cond2 <- contrastconditions[2]

# split data frame into RIBO and RNA
RIBO <- counts[, (sampleSheet$method == "RIBO") & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2)]
RNA <- counts[, (sampleSheet$method == "RNA")  & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2)]

#contrastconditions <- c( "TruDWT", "TruDWT", "TruDDel", "TruDDel")
numberofreplicates <- max(sampleSheet$replicate)
contrastconditionsvector <- rep(contrastconditions,each=numberofreplicates)
# run xtail analysis
results.deseq2 <- riborex(RNA, RIBO, contrastconditionsvector, contrastconditionsvector)
#results.edgeRD <- riborex(RNA, RIBO, contrastconditionsvector, contrastconditionsvector, "edgeRD")
#results.voom <- riborex(RNA, RIBO, contrastconditionsvector, contrastconditionsvector, "Voom")
summary(results.deseq2)
#summary(results.edgeRD)
#summary(results.voom)
# write results into file
write.csv(results.deseq2, options$riborexdeseq_result_path, quote = F)
#write.csv(test_tab, options$riborexedgeR_result_path, quote = F)
#write.csv(test_tab, options$riborexvoom_result_path, quote = F)

#plot results

