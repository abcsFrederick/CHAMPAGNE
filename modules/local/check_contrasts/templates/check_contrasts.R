#!/usr/bin/env Rscript
library(assertthat)
library(dplyr)
library(glue)

main <- function(contrasts_filename = "${contrasts}",
                 samplesheet_filename = "${samplesheet}",
                 output_basename = "${output_basename}",
                 versions_filename = "versions.yml") {
  write_version(versions_filename)
  contrasts_lst <- yaml::read_yaml(contrasts_filename)
  samples_df <- readr::read_csv(samplesheet_filename)

  assert_that(
    all(unlist(contrasts_lst) %in% (samples_df %>% dplyr::pull(sample))),
    glue("All sample names in contrasts must also be in the sample sheet")
  )
  names(contrasts_lst) %>% lapply(check_contrast, contrasts_lst)
  names(contrasts_lst) %>% lapply(
    write_contrast_samplesheet,
    contrasts_lst, samples_df, output_basename
  )
}

check_contrast <- function(contrast_name, contrasts_lst) {
  contrast_len <- length(contrasts_lst[[contrast_name]])
  assert_that(
    contrast_len == 2,
    glue("Contrasts must have only two groups per comparison, but {contrast_name} has {contrast_len}")
  )
}

write_contrast_samplesheet <- function(contrast_name, contrasts_lst, samples_df, output_basename) {
  contrast_df <- contrasts_lst[[contrast_name]] %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "group", values_to = "sample"
    )
  samples_df %>%
    dplyr::inner_join(contrast_df) %>%
    readr::write_csv(glue("{output_basename}.{contrast_name}.csv"))
}

get_version <- function() {
  return(paste0(R.version[["major"]], ".", R.version[["minor"]]))
}

write_version <- function(version_file) {
  write_lines(get_version(), version_file)
}

main()
