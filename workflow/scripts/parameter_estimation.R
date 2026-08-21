#!/usr/bin/env Rscript

# Estimates the RPKM and coverage cutoffs DeepRibo needs, using the S-curve
# method that ships inside the DeepRibo container. The sourced path only exists
# there, which is why the calling rule uses container: rather than conda:.

library("optparse")

source('/usr/local/bin/s_curve_cutoff_estimation.R')

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="data_list.csv produced by DataParser.py", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="file to write the two cutoffs to", metavar="character"),
  make_option(c("-d", "--dest"), type="character", default=NULL,
              help="prefix for the S-curve diagnostic figure. Must differ between
                    libraries: a shared prefix makes concurrent jobs overwrite
                    each other's figure.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file) || is.null(opt$out)){
  print_help(opt_parser)
  stop("Both an input (-f) and an output (-o) file must be supplied.", call.=FALSE)
}

# Default the figure next to the output, so that two libraries running at the
# same time cannot write to the same path.
dest <- if (is.null(opt$dest)) file.path(dirname(opt$out), "s_curve") else opt$dest

parameters <- get_cutoff_values(path=opt$file, dest=dest)

output <- paste(toString(parameters[["min_RPKM"]]),
                toString(parameters[["min_coverage"]]),
                sep=",", collapse=NULL)

write(output, file=opt$out)
