#!/usr/bin/env Rscript
library("optparse")
library("dplyr")


option_list = list(
  make_option("--col1", type="character", default=NULL,
              help="exon start site", metavar="character"),
    make_option("--col2", type="character", default=NULL,
              help="exon end site", metavar="character"),
        make_option("--out", type="character", default=NULL,
              help="metagene coordinates of start and end sites", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$col1)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

c1 <- read.delim(opt$col1, header = FALSE)
c2 <- read.delim(opt$col2, header = FALSE)

convert <- function(c1, c2){
    diff <- as.vector(unlist(c2-c1))
    cum_diff <- cumsum(diff+1)
    c4 <- cum_diff - 1
    print(cum_diff)
    c3 <- c(1,cum_diff[1:length(cum_diff)-1])
    return(data.frame(c3, c4))
}

d <- convert(c1, c2)

write.table(d, file = opt$out, sep = "\t", row.names = F, col.names = F, quote = F)
