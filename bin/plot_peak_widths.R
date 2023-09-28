#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(readr)

set.seed(20230926)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: not enough positional arguments.\n  Example usage:\n   Rscript plot_peak_widths.R path/to/peak_widths.tsv")
}
tsv_filename <- args[1]
peak_dat <- read_tsv(tsv_filename)

peak_widths_hist <- peak_dat %>%
  mutate(peak_width = chromEnd - chromStart) %>%
  ggplot(aes(peak_width, fill = tool)) +
  geom_histogram(alpha = 0.7, position = "identity") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(
    label.position = "bottom",
    title = "Peak caller",
    title.position = "top"
  )) +
  facet_wrap("sample_id") +
  labs(x = expression("" * log[10] * " Peak Widths")) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(filename = "peak_widths_histogram.png", plot = peak_widths_hist, device = "png")
