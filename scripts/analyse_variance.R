#!/usr/bin/env Rscript

library(ggplot2)
library(DESeq2)

library(optparse)
library(plyr)

option_list = list(
  make_option(c("-r", "--read_count"), type = "character", default = NULL,
              help = "Path to read counts", metavar = "character"),
  make_option(c("-m", "--meta_table"), type = "character", default = NULL,
              help = "Path to meta data table", metavar = "character"),
  make_option(c("-o", "--output_path"), type = "character", default = NULL,
              help = "Output path", metavar = "character")
);

option_parser = OptionParser(option_list = option_list);
options = parse_args(option_parser);

if (is.null(options$read_count) || is.null(options$meta_table) || is.null(options$output_path)) {
  print_help(option_parser)
  stop("Please supply arguments (-r, -m, -o), see --help \n", call.=FALSE)
}

## try annotated genes as deepribo predictions might be wrong for some organisms

# Load data
data <- read.csv(options$read_count, header = TRUE, row.names = 1)
meta <- read.csv(options$meta_table, header = TRUE, sep="\t", row.names = 1)

class(data)
class(meta)


pdf(paste(options$output_path, "raw_count_distributions.pdf",sep=""))
for (i in 1:ncol(data)) {
    print(ggplot(data, aes(data[,i])) + geom_histogram(stat = "bin", bins = 200) + xlab("Raw expression counts") + ylab("Number of genes") + ggtitle(names(data)[i]))
}
dev.off()

pdf(paste(options$output_path, "mean_vs_variance.pdf",sep=""))
col_dfs <- split.default(data, sub("_\\d+", "", names(data)))
for (x in col_dfs) {
    mean_counts <- apply(x, 1, mean)
    variance_counts <- apply(x, 1, var)

    df <- data.frame(mean_counts, variance_counts)

    name_parts <- strsplit(names(x)[1], "_", fixed = TRUE)
    title <- paste(name_parts[[1]][1], name_parts[[1]][2], sep = "-")

    print(ggplot(df) + geom_point(aes(x=mean_counts, y=variance_counts)) + geom_line(aes(x=mean_counts, y=mean_counts, color="red")) + scale_y_log10() + scale_x_log10() + ggtitle(title))
}
dev.off()

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~sampletype)
dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

write.table(data.frame("Identifier"=rownames(normalized_counts), normalized_counts), file = paste(options$output_path, "normalized_counts.tsv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)


rld <- rlog(dds, blind = TRUE)

write_pca_info <- function (object, intgroup = "condition", ntop = 500){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))

  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])

  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))

  write.table(d, file = paste(options$output_path, "rld.tsv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(percentVar, file = paste(options$output_path, "variance_percentages.tsv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_pca_info(rld, intgroup = "sampletype")
