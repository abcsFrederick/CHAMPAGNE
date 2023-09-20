#!/usr/bin/env Rscript --vanilla

library(broom)
library(dplyr)
library(ggplot2)
library(glue)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: not enough positional arguments.\n\tExample usage:\n\t\tRscript plot_jaccard.R path/to/jaccard.txt")
}
jaccard_tsv <- args[1]
jaccard_dat <- read_tsv(jaccard_tsv)
peak_callers <- c(
  jaccard_dat %>% pull(toolA),
  jaccard_dat %>% pull(toolB)
) %>%
  unique()

# plot PCA per peak caller, jaccard on different samples
per_tool <- peak_callers %>% sapply(function(tool) {
  dat <- jaccard_dat %>%
    filter(toolA == tool & toolB == tool & fileA != fileB)
  pca <- dat %>%
    select(jaccard) %>%
    prcomp(scale = TRUE)
  pct_var <- pca %>%
    tidy(matrix = "eigenvalues")
  plot_pca <- pca %>%
    augment(dat) %>%
    ggplot(aes(.fittedPC1, .fittedPC2)) +
    geom_point() +
    labs(
      x = glue("PC1 ({pct_var %>% filter(PC==1) %>% pull(percent)})"),
      y = glue("PC2 ({pct_var %>% filter(PC==2) %>% pull(percent)})")
    )
  ggsave(plot = plot_pca, filename = glue("jaccard_pca_{tool}_mqc.png"), device = "png")
})

# one plot with all tools on same samples
dat <- jaccard_dat %>%
  filter(toolA != toolB & fileA == fileB)
pca <- dat %>%
  select(jaccard) %>%
  prcomp(scale = TRUE)
pct_var <- pca %>%
  tidy(matrix = "eigenvalues")
plot_pca <- pca %>%
  augment(dat) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color = labelA)) +
  geom_point() +
  labs(
    x = glue("PC1 ({pct_var %>% filter(PC==1) %>% pull(percent)})"),
    y = glue("PC2 ({pct_var %>% filter(PC==2) %>% pull(percent)})")
  )
ggsave(plot = plot_pca, filename = glue("jaccard_pca_all_mqc.png"), device = "png")
