#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(readr)

set.seed(20230922)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: not enough positional arguments.\n  Example usage:\n   Rscript plot_frip.R path/to/frips.txt")
}

frip_filename <- args[1]
frip_dat <- read_tsv(frip_filename)

sample_frip_plot <- frip_dat %>%
  ggplot(aes(FRiP, bedsample, color = bedtool)) +
  geom_jitter(alpha = 0.8, height = 0.1) +
  geom_hline(
    yintercept = seq(1.5, length(unique(frip_dat %>% pull(bedsample))) - 0.5, 1),
    lwd = 0.5, color = "grey92"
  ) +
  guides(color = guide_legend(title = "Peak caller")) +
  labs(x = "Fraction of Reads in Peaks", y = "") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())
ggsave(filename = "samples_FRiP_mqc.png", plot = sample_frip_plot, device = "png")

nbasesM_frip_plot <- frip_dat %>%
  ggplot(aes(n_basesM, FRiP, color = bedsample, shape = bedtool)) +
  geom_point(alpha = 0.8) +
  guides(
    color = guide_legend(title = "Sample"),
    shape = guide_legend(title = "Peak caller")
  ) +
  labs(x = "Number of Bases (millions) in Peaks", y = "Fraction of Reads in Peaks") +
  theme_bw()
ggsave(filename = "nbasesM_FRiP_mqc.png", plot = nbasesM_frip_plot, device = "png")
