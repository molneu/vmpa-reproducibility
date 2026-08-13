#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript run_supplementary_figure9.R <annotated-raw-counts.csv> <western-blot.csv> <repository-root>",
    call. = FALSE
  )
}

raw_file <- normalizePath(args[1], mustWork = TRUE)
wb_file <- normalizePath(args[2], mustWork = TRUE)
repo_root <- normalizePath(args[3], mustWork = TRUE)

figure6_dir <- file.path(repo_root, "reproducibility", "figure6")
score_dir <- file.path(figure6_dir, "data", "wb_correlation", "vmpa_wb_correlations")
wb_dir <- file.path(figure6_dir, "data", "wb_correlation")

dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wb_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(wb_file, file.path(wb_dir, "260725 WB_corrected_data_with Ponceau.csv"), overwrite = TRUE)

score_script <- file.path(figure6_dir, "scripts", "generate_temsirolimus_vmpa_scores.R")
status <- system2(
  "Rscript",
  shQuote(c(score_script, raw_file, score_dir))
)
if (status != 0L) stop("VMPA score generation failed.", call. = FALSE)

old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(repo_root)

source(file.path(figure6_dir, "scripts", "wb_correlation", "analyze_original_wb_mtor_s6_across_cell_lines.R"))
source(file.path(figure6_dir, "scripts", "wb_correlation", "make_four_mtor_validation_scatterplots.R"))

selected_dir <- file.path(
  figure6_dir, "results", "wb_correlation",
  "original_wb_labels_mtor_s6_n6_n12", "SELECTED_FINAL_PLOTS"
)
plot_dir <- file.path(repo_root, "reproducibility", "supplementary_figure9", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
files <- c(
  "BN118_MTOR_vs_pMTOR_Ponceau_50mm.pdf",
  "BN118_MTOR_vs_pS6_Ponceau_50mm.pdf",
  "BN118_RPS6.1_vs_pS6_Ponceau_50mm.pdf",
  "BN91_MTOR_vs_pMTOR_Ponceau_50mm.pdf",
  "BN91_MTOR_vs_pS6_Ponceau_50mm.pdf",
  "BN91_RPS6.1_vs_pS6_Ponceau_50mm.pdf",
  "supplementary_six_regression_results.csv"
)
copied <- file.copy(
  file.path(selected_dir, files),
  file.path(plot_dir, files),
  overwrite = TRUE
)
if (!all(copied)) stop("Failed to update Supplementary Figure 9 outputs.", call. = FALSE)
