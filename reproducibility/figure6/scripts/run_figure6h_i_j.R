#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript run_figure6h_i_j.R <raw-counts.csv> <ensembl-symbol-mapping.rds> <western-blot.csv> <repository-root>",
    call. = FALSE
  )
}

raw_file <- normalizePath(args[1], mustWork = TRUE)
mapping_file <- normalizePath(args[2], mustWork = TRUE)
wb_file <- normalizePath(args[3], mustWork = TRUE)
repo_root <- normalizePath(args[4], mustWork = TRUE)
figure6_dir <- file.path(repo_root, "reproducibility", "figure6")
score_dir <- file.path(
  figure6_dir, "data", "wb_correlation", "vmpa_wb_correlations"
)
wb_dir <- file.path(figure6_dir, "data", "wb_correlation")

dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wb_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(
  wb_file,
  file.path(wb_dir, "260725 WB_corrected_data_with Ponceau.csv"),
  overwrite = TRUE
)

run_script <- function(path, arguments = character()) {
  status <- system2("Rscript", shQuote(c(path, arguments)))
  if (status != 0L) {
    stop("Script failed: ", path, call. = FALSE)
  }
}

run_script(
  file.path(figure6_dir, "scripts", "generate_temsirolimus_vmpa_scores.R"),
  c(raw_file, mapping_file, score_dir)
)

old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(repo_root)

run_script(file.path(
  figure6_dir, "scripts", "wb_correlation",
  "analyze_original_wb_mtor_s6_across_cell_lines.R"
))
run_script(file.path(
  figure6_dir, "scripts", "wb_correlation",
  "make_four_mtor_validation_scatterplots.R"
))

target_scripts <- file.path(
  figure6_dir,
  "scripts",
  "target_selection",
  c(
    "prepare_activity_matrix_from_raw_result.R",
    "make_vmpa_treatment_response_tables.R",
    "make_vmpa_average_temsi_response.R",
    "rank_vmpa_candidates_n250_sample_pop_sd.R"
  )
)
invisible(lapply(target_scripts, run_script))
run_script(
  file.path(
    figure6_dir, "scripts", "target_selection",
    "vmpa_initial_unbiased_targetable_screen_n250_sample_pop_sd.R"
  ),
  "positive"
)
run_script(file.path(
  figure6_dir, "scripts", "target_selection",
  "plot_figure6i_temsirolimus_candidates.R"
))
run_script(file.path(
  figure6_dir, "scripts", "apoptosis_pathway",
  "make_compass_intrinsic_apoptosis_pathway.R"
))
