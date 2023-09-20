#!/usr/bin/env Rscript --vanilla

library(broom)
library(dplyr)
library(ggplot)
library(glue)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: not enough positional arguments.\n\tExample usage:\n\t\tRscript plot_jaccard.R path/to/jaccard.txt")
}
jaccard_tsv <- args[1]
jaccard_dat <- read_tsv(jaccard_tsv)
pca_all <- jaccard_dat %>%
  select(jaccard) %>%
  prcomp(scale = TRUE)
pct_var <- pca_all %>%
  tidy(matrix = "eigenvalues")
plot_pca_all <- pca_all %>%
  augment(jaccard_dat) %>%
  ggplot(aes(.fittedPC1, .fittedPC2)) +
  geom_point() +
  labs(
    x = glue("PC1 ({pct_var %>% filter(PC==1) %>% pull(percent)})"),
    y = glue("PC2 ({pct_var %>% filter(PC==2) %>% pull(percent)})")
  )
ggsave(plot = plot_pca_all, filename = "jaccard_all_pca.png", device = "png")
