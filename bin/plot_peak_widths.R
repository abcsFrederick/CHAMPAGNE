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
peak_dat <- read_tsv(tsv_filename) %>%
  mutate(peak_width = chromEnd - chromStart)
xmin <- peak_dat %>%
  pull(peak_width) %>%
  min()
xmax <- peak_dat %>%
  pull(peak_width) %>%
  max()

peak_widths_hist <- peak_dat %>%
  ggplot(aes(peak_width, fill = tool)) +
  geom_histogram(alpha = 0.7, position = "identity") +
  scale_x_log10(
    limits = c(xmin - 10^2, xmax + 10^2),
    labels = scales::label_log(digits = 2)
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

ggsave(filename = "peak_widths_histogram.png", plot = peak_widths_hist, device = "png", dpi = 300, height = 4, width = 5)
