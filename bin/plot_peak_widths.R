#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(readr)
library(scales)

set.seed(20230926)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: not enough positional arguments.\n  Example usage:\n   Rscript plot_peak_widths.R path/to/peak_widths.tsv")
}
tsv_filename <- args[1]
peak_dat <- read_tsv(tsv_filename) %>%
  mutate(
    peak_width = chromEnd - chromStart,
    peak_type = case_when(
      tool == "sicer" | tool == "macs_broad" ~ "broad",
      tool == "gem" | tool == "macs_narrow" ~ "narrow",
      TRUE ~ NA_character_
    )
  )

tool_colors <- viridisLite::viridis(n = 4, alpha = 0.7, begin = 0, end = 1, direction = 1, option = "D")
names(tool_colors) <- peak_dat %>%
  pull(tool) %>%
  unique()

plot_hist <- function(peak_dat) {
  xmin <- peak_dat %>%
    pull(peak_width) %>%
    min()
  xmax <- peak_dat %>%
    pull(peak_width) %>%
    max()
  peak_dat %>%
    ggplot(aes(peak_width, fill = tool)) +
    geom_histogram(alpha = 0.7, position = "identity") +
    scale_x_log10(
      limits = c(xmin - 10^2, xmax + 10^2),
      labels = label_log(digits = 2)
    ) +
    scale_fill_manual(
      values = tool_colors,
      breaks = names(tool_colors)
    ) +
    guides(fill = guide_legend(
      label.position = "bottom",
      title = "Peak caller",
      title.position = "top"
    )) +
    facet_wrap(~sample_id) +
    labs(x = expression("" * log[10] * " Peak Widths")) +
    theme_bw() +
    theme(legend.position = "top")
}

hist_broad <- peak_dat %>%
  filter(peak_type == "broad") %>%
  plot_hist() +
  scale_y_log10(labels = label_log(digits = 2))
hist_narrow <- peak_dat %>%
  filter(peak_type == "narrow") %>%
  plot_hist() +
  scale_y_log10(labels = label_log(digits = 2))

ggsave(
  filename = "peak_widths_broad_histogram.png", plot = hist_broad,
  device = "png", dpi = 300, height = 4, width = 5
)
ggsave(
  filename = "peak_widths_narrow_histogram.png", plot = hist_narrow,
  device = "png", dpi = 300, height = 4, width = 5
)
