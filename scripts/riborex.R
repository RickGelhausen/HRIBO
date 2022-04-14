#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-r", "--raw_read_counts_csv_path"), type = "character", default = NULL,
              help = "Path to raw read counts table", metavar = "character"),
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
print(sampleSheet)
sampleSheet <- sampleSheet[
  order( sampleSheet[,1], sampleSheet[,2], sampleSheet[,3] ),
]
print(sampleSheet)

#create condition vector
constraststring <- gsub("contrasts/", "", options$contrast)

contrastconditions <- unlist(strsplit(constraststring,"-"))
cond1 <- contrastconditions[1]
cond2 <- contrastconditions[2]
print(cond1)
print(cond2)
# split data frame into RIBO and RNA
RIBO <- counts[, (sampleSheet$method == "RIBO") & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2)]
RNA <- counts[, (sampleSheet$method == "RNA")  & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2)]

head(count,5)
print((sampleSheet$method == "RIBO") & ( sampleSheet$condition == cond1 | sampleSheet$condition == cond2))

head(RIBO,5)
head(RNA,5)
countsheader <- colnames(counts)
countsheader <- countsheader[grepl("RIBO", countsheader)]
replicatescondition1 <- length(grep(paste("-",cond1,"-",sep=""), countsheader))
replicatescondition2 <- length(grep(paste("-",cond2,"-",sep=""), countsheader))

conditionsvector1 <- rep(cond1,each=replicatescondition1)
conditionsvector2 <- rep(cond2,each=replicatescondition2)
contrastconditionsvector <- c(conditionsvector1, conditionsvector2)

# run riborex analysis
results.deseq2 <- riborex(RNA, RIBO, contrastconditionsvector, contrastconditionsvector)

summary(results.deseq2)
# write results into file
write.csv(results.deseq2, options$riborexdeseq_result_path, quote = F)

#plot results
