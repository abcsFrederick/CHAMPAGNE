#!/usr/bin/env Rscript
load_package <- function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
  invisible(x)
}
messages <- lapply(c("ChIPseeker", "dplyr", "glue", "ggplot2"), load_package)

parser <- argparse::ArgumentParser()
parser$add_argument("-p", "--peak", required = TRUE, help = "peak file")
parser$add_argument("-o", "--outfile-prefix", required = TRUE, type = "character", dest = "outfile_prefix", help = "prefix for output filenames")
parser$add_argument("--genome-txdb", dest = "txdb", required = TRUE, help = "BioConductor TxDb package, e.g. TxDb.Hsapiens.UCSC.hg38.knownGene")
parser$add_argument("--genome-annot", dest = "adb", required = TRUE, help = "BioConductor annotation package, e.g. org.Hs.eg.db")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
outfile_prefix <- args$outfile_prefix
adb <- args$adb
load_package(adb)
load_package(args$txdb)
txdb <- args$txdb %>%
  rlang::sym() %>%
  eval()

np <- read.table(args$peak, sep = "\t")
peak_colnames <- c(
  "chrom",
  "chromStart",
  "chromEnd",
  "name",
  "score",
  "strand",
  "signalValue",
  "pValue",
  "qValue"
)
num_columns <- ncol(np)
if (num_columns == 9) {
  colnames(np) <- peak_colnames
  np <- np %>% mutate(peak = NA)
} else if (num_columns == 10) {
  colnames(np) <- c(peak_colnames, "peak")
} else {
  stop(paste("Expected 9 or 10 columns in peak file, but", num_columns, "given"))
}
# plots for individual peak file
peaks <- GenomicRanges::GRanges(
  seqnames = np$chrom,
  ranges = IRanges(np$chromStart, np$chromEnd),
  qValue = np$qValue
)
plots <- list(
  covplot = covplot(peaks, weightCol = "qValue"),
  plotPeakProf2 = plotPeakProf2(
    peak = peaks, upstream = rel(0.2), downstream = rel(0.2),
    conf = 0.95, by = "gene", type = "body", nbin = 800,
    TxDb = txdb, weightCol = "qValue", ignore_strand = F
  )
)

names(plots) %>%
  lapply(function(plot_name) {
    ggsave(
      filename = glue("{outfile_prefix}_{plot_name}.png"),
      plot = plots[[plot_name]],
      device = "png"
    )
  })
