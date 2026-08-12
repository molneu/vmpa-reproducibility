#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript run_supplementary_figure9.R <raw-counts.csv> <ensembl-symbol-mapping.rds> <western-blot.csv> <repository-root>",
    call. = FALSE
  )
}

raw_file <- normalizePath(args[1], mustWork = TRUE)
mapping_file <- normalizePath(args[2], mustWork = TRUE)
wb_file <- normalizePath(args[3], mustWork = TRUE)
repo_root <- normalizePath(args[4], mustWork = TRUE)

figure6_dir <- file.path(repo_root, "reproducibility", "figure6")
score_dir <- file.path(figure6_dir, "data", "wb_correlation", "vmpa_wb_correlations")
wb_dir <- file.path(figure6_dir, "data", "wb_correlation")

dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wb_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(wb_file, file.path(wb_dir, "260725 WB_corrected_data_with Ponceau.csv"), overwrite = TRUE)

score_script <- file.path(figure6_dir, "scripts", "generate_temsirolimus_vmpa_scores.R")
status <- system2(
  "Rscript",
  shQuote(c(score_script, raw_file, mapping_file, score_dir))
)
if (status != 0L) stop("VMPA score generation failed.", call. = FALSE)

old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(repo_root)

source(file.path(figure6_dir, "scripts", "wb_correlation", "analyze_original_wb_mtor_s6_across_cell_lines.R"))
source(file.path(figure6_dir, "scripts", "wb_correlation", "make_four_mtor_validation_scatterplots.R"))
