#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-b", "--bam_directory_path"), type = "character", default = NULL,
              help = "Path to directory containing .bam files", metavar="character"),
  make_option(c("-a", "--annotation_file_path"), type = "character", default = NULL,
              help = "Path to .gtf file with annotation", metavar = "character"),
  make_option(c("-t", "--sample_file_path"), type = "character", default = NULL,
              help = "Path to sample.tsv", metavar = "character"),
  make_option(c("-r", "--raw_count_out_path"), type = "character", default = NULL,
              help = "Path for writing output normalized count file", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$bam_directory_path)){
  print_help(option_parser)
  stop("Please supply arguments (-b, -a, -t, -r), see --help \n", call.=FALSE)
}

library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(plyr)
library(rtracklayer)

# import longest protein coding transcripts
gencode <- import.gff(options$annotation_file_path)

# keep only cds
sel <- gencode$type == "CDS"
cds <- gencode[which(sel),]

# create empty data frame
gene.counts <- data.frame(gene.id = cds$transcript_id)

# get sample files
sample.files <- paste(options$bam_directory_path, list.files(options$bam_directory_path, pattern = "\\.bam$"), sep = "")

for (i in sample.files) {

  # get sample name
  name.i <- as.character(i)

  # import reads
  reads <- readGAlignments(i)

  # convert to granges
  reads <- granges(reads)

  # keep only first nt
  reads <- flank(reads, -1)

  # count reads into genes
  gene.counts[, name.i] <- countOverlaps(cds, reads)

}

# sum all reads across gene.ids
gene.counts <- ddply(gene.counts,"gene.id",numcolwise(sum))

# change row names and drop column gene.id
rownames(gene.counts) <- gene.counts$gene.id
gene.counts$gene.id <- NULL

# get sample sheet
sampleSheet <- read.csv(file=options$sample_file_path ,header=TRUE, sep="\t", stringsAsFactors=FALSE)

#col names and row names 
sampleName <- function(x) {
  method <- x[1]
  condition <- x[2]
  replicate <- x[3]
  sname <- paste(method, condition, replicate, sep="-")
  return(sname)
}

# generate sampleTable
sampleNames <- apply(sampleSheet,1,sampleName)
sample.files <- paste(sampleNames, ".bam", sep="")
conditions <- sampleSheet[,2]
sampleTable <- data.frame(row.names = sampleNames, fileName = sample.files, condition = conditions)

colnames(gene.counts) <- rownames(sampleTable)


# create DESeq data set
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = sampleTable,
                              design = ~ condition)
# get normalized counts
raw.counts <-  counts(dds, normalized = FALSE)
# save normalized counts
write.csv(raw.counts, options$raw_count_cds_out_path, quote = F)

