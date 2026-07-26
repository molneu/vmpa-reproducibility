#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

data_dir <- here::here("reproducibility", "figure6", "data", "wb_correlation", "vmpa_wb_correlations")
out_dir <- here::here("reproducibility", "figure6", "results", "target_selection", "vmpa_wb_correlations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "vmpa_glioma_n250_sample_pop_sd_unique_FALSE_raw_result.csv")
output_file <- file.path(out_dir, "vmpa_glioma_n250_sample_pop_sd_unique_FALSE_activity_matrix.csv")

if (!file.exists(input_file)) {
  stop("Missing VMPA raw result file: ", input_file)
}

raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("target", "conf", "signature")
if (!all(required %in% names(raw))) {
  stop("Raw result file must contain columns: ", paste(required, collapse = ", "))
}

sample_columns <- setdiff(names(raw), required)
activity <- data.frame(
  activity_target = make.unique(as.character(raw$target)),
  raw[, sample_columns, drop = FALSE],
  check.names = FALSE
)

write.csv(activity, output_file, row.names = FALSE)
message("Wrote activity matrix: ", output_file)
