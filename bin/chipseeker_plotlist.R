#!/usr/bin/env Rscript
# chipseeker plots that accept list of annotation objects
load_package <- function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
  invisible()
}
lapply(c("ChIPseeker", "dplyr", "ggplot2", "glue"), load_package)

parser <- argparse::ArgumentParser()
parser$add_argument("-a", "--annotations", required = TRUE, type = "character", dest = "annots", nargs = "+", help = "space-delimited list of Rds files of peak annotation objects from chipseeker")
parser$add_argument("-o", "--outfile", required = TRUE, type = "character", help = "output filenames")

args <- parser$parse_args()

peak_anno_list <- args$annotations %>%
  stringr::str_split(" ") %>%
  lapply(readRDS)
anno_bar_plot <- plotAnnoBar(peak_anno_list)
ggsave(filename = args$outfile, plot = anno_bar_plot)
