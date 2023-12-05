#!/usr/bin/env Rscript
library(assertthat)
library(dplyr)
library(glue)

main <- function(contrasts_filename = "${contrasts}",
                 samplesheet_filename = "${samplesheet}",
                 output_basename = "${output_basename}",
                 versions_filename = "versions.yml",
                 process_name = "${task.process}") {
  write_version(versions_filename, process_name = process_name)
  contrasts_lst <- yaml::read_yaml(contrasts_filename)
  samples_df <- readr::read_csv(samplesheet_filename)

  assert_that(
    all(unlist(contrasts_lst) %in% (samples_df %>% dplyr::pull(sample_basename))),
    msg = glue("All sample names in contrasts must also be in the sample sheet")
  )
  names(contrasts_lst) %>% lapply(check_contrast, contrasts_lst)
  names(contrasts_lst) %>% lapply(
    write_contrast_samplesheet,
    contrasts_lst, samples_df, output_basename
  )
}

#' Check an individual contrast comparison
check_contrast <- function(contrast_name, contrasts_lst) {
  # Ensure contrast has exactly two groups to compare
  contrast_len <- length(contrasts_lst[[contrast_name]])
  assert_that(
    contrast_len == 2,
    msg = glue("Contrasts must have only two groups per comparison, but {contrast_name} has {contrast_len}")
  )
  # TODO: check that each sample only appears in one group
}

#' Combine sample sheet with contrast group
write_contrast_samplesheet <- function(contrast_name, contrasts_lst, samples_df, output_basename) {
  contrast_df <- contrasts_lst[[contrast_name]] %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "group", values_to = "sample_basename"
    ) %>%
    mutate(contrast = contrast_name)
  print(head(contrast_df))
  samples_df %>%
    dplyr::inner_join(contrast_df, by = "sample_basename") %>%
    readr::write_csv(glue("{output_basename}.{contrast_name}.csv"))
}

#' Get R version as a semantic versioning string
get_r_version <- function() {
  return(paste0(R.version[["major"]], ".", R.version[["minor"]]))
}

#' Write R version to a yaml file
write_version <- function(version_file, process_name = "${task.process}") {
  readr::write_lines(
    glue('"{process_name}":',
      "R: {get_r_version()}",
      .sep = "\\n\\t"
    ),
    version_file
  )
}

main()
