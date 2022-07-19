#!/usr/bin/env Rscript

library(optparse)
library(plyr)

option_list = list(
  make_option(c("-r", "--raw_read_counts_csv_path"), type = "character", default = NULL,
              help = "Path to read counts table", metavar = "character"),
  make_option(c("-c", "--contrast"), type = "character", default = NULL,
              help = "Contrast, pair of conditions ", metavar = "character"),
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

if (is.null(options$raw_read_counts_csv_path)){
  print_help(option_parser)
  stop("Please supply arguments (-r, -t, -x), see --help \n", call.=FALSE)
}

library(xtail)

# read table with raw read counts
counts <- read.csv(options$raw_read_counts_csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
sampleSheet <- sampleSheet[
  order( sampleSheet[,1], sampleSheet[,2], sampleSheet[,3] ),
]

#create condition vector
constraststring <- gsub("contrasts/", "", options$contrast)
contrastconditions <- unlist(strsplit(constraststring, "-"))
cond1 <- contrastconditions[1]
cond2 <- contrastconditions[2]
print(paste("contrastconditions: ", cond1, " vs ", cond2))

  # split data frame into RIBO and RNA
RIBO_c1 <- counts[, (sampleSheet$method == "RIBO") & ( sampleSheet$condition == cond1)]
RIBO_c1 <- cbind(identifier = rownames(RIBO_c1), RIBO_c1)
rownames(RIBO_c1) <- 1:nrow(RIBO_c1)
RIBO_c2 <- counts[, (sampleSheet$method == "RIBO") & ( sampleSheet$condition == cond2)]
RIBO_c2 <- cbind(identifier = rownames(RIBO_c2), RIBO_c2)
rownames(RIBO_c2) <- 1:nrow(RIBO_c2)

RIBO <- join(RIBO_c1, RIBO_c2, by = "identifier")
rownames(RIBO) <- RIBO$identifier
RIBO <- subset(RIBO, select = -c(identifier))


RNA_c1 <- counts[, (sampleSheet$method == "RNA")  & ( sampleSheet$condition == cond1)]
RNA_c1 <- cbind(identifier = rownames(RNA_c1), RNA_c1)
rownames(RNA_c1) <- 1:nrow(RNA_c1)
RNA_c2 <- counts[, (sampleSheet$method == "RNA")  & ( sampleSheet$condition == cond2)]
RNA_c2 <- cbind(identifier = rownames(RNA_c2), RNA_c2)
rownames(RNA_c2) <- 1:nrow(RNA_c2)

RNA <- join(RNA_c1, RNA_c2, by = "identifier")
rownames(RNA) <- RNA$identifier
RNA <- subset(RNA, select = -c(identifier))


countsheader <- colnames(counts)
countsheader <- countsheader[grepl("RIBO", countsheader)]
replicatescondition1 <- length(grep(paste("-",cond1,"-",sep=""), countsheader))
replicatescondition2 <- length(grep(paste("-",cond2,"-",sep=""), countsheader))

#numberofreplicates <- max(sampleSheet$replicate)
#contrastconditionsvector <- rep(contrastconditions,each=numberofreplicates)
conditionsvector1 <- rep(cond1,each=replicatescondition1)
conditionsvector2 <- rep(cond2,each=replicatescondition2)
contrastconditionsvector <- c(conditionsvector1, conditionsvector2)
print(contrastconditionsvector)

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
