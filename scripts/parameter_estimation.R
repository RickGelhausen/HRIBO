#!/usr/bin/env Rscript

library("optparse")

# load the functions from the script
source('/usr/local/bin/s_curve_cutoff_estimation.R')

# commandline parser
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input / output file).n", call.=FALSE)
}

# list the dataset and the path to which the png figure is stored
parameters <- get_cutoff_values(path=opt$file, dest="figure")

output <- paste(toString(parameters[["min_RPKM"]]),toString(parameters[["min_coverage"]]), sep=",", collapse=NULL)

write(output, file=opt$out)
