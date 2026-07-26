#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
})

data_dir <- here::here("reproducibility", "figure6", "data", "wb_correlation")
results_dir <- here::here("reproducibility", "figure6", "results", "wb_correlation")

wb_file <- file.path(data_dir, "260725 WB_corrected_data_with Ponceau.csv")
vmpa_dir <- file.path(data_dir, "vmpa_wb_correlations")
out_dir <- file.path(
  results_dir,
  "original_wb_labels_mtor_s6_n6_n12"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scalings <- c("none", "sample_z", "sample_pop_sd", "signature_z")
vmpa_targets <- c("MTOR", "RPS6", "RPS6KB2", "RPS6KA1")
cell_order <- c("BN91", "BN91R", "BN118", "BN118R")
treatment_order <- c("DMSO", "50nM", "500nM")

clean_text <- function(x) trimws(as.character(x))

condition_to_treatment <- function(x) {
  result <- toupper(clean_text(x))
  result[result == "T50"] <- "50nM"
  result[result == "T500"] <- "500nM"
  result
}

z_score <- function(x) {
  if (length(x) < 2L || any(!is.finite(x)) || sd(x) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

pearson_values <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(c(n = sum(keep), r = NA_real_, p = NA_real_))
  }
  test <- cor.test(x[keep], y[keep], method = "pearson")
  c(n = sum(keep), r = unname(test$estimate), p = test$p.value)
}

permutations_three <- rbind(
  c(1L, 2L, 3L),
  c(1L, 3L, 2L),
  c(2L, 1L, 3L),
  c(2L, 3L, 1L),
  c(3L, 1L, 2L),
  c(3L, 2L, 1L)
)

blocked_permutation_p <- function(x, y, cell_line) {
  cell_levels <- unique(cell_line)
  if (any(table(cell_line) != 3L) ||
      any(!is.finite(x)) ||
      any(!is.finite(y))) {
    return(NA_real_)
  }
  observed <- abs(cor(x, y))
  choices <- expand.grid(rep(
    list(seq_len(nrow(permutations_three))),
    length(cell_levels)
  ))
  n_permutations <- nrow(choices)
  permuted_x <- matrix(
    rep(x, each = n_permutations),
    nrow = n_permutations,
    ncol = length(x)
  )
  for (j in seq_along(cell_levels)) {
    idx <- which(cell_line == cell_levels[j])
    permutation_indices <- permutations_three[choices[[j]], , drop = FALSE]
    permuted_x[, idx] <- vapply(
      seq_len(ncol(permutation_indices)),
      function(column) x[idx][permutation_indices[, column]],
      numeric(n_permutations)
    )
  }
  centered_x <- sweep(permuted_x, 1L, rowMeans(permuted_x), FUN = "-")
  centered_y <- y - mean(y)
  numerator <- as.vector(centered_x %*% centered_y)
  denominator <- sqrt(rowSums(centered_x^2) * sum(centered_y^2))
  permuted_r <- abs(numerator / denominator)
  mean(permuted_r >= observed - sqrt(.Machine$double.eps))
}

raw_wb <- read.csv(
  wb_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw_wb <- raw_wb[, 1:5, drop = FALSE]
names(raw_wb) <- c("target", "membrane", "condition", "cell_line", "band_value")
raw_wb$target <- clean_text(raw_wb$target)
raw_wb$membrane <- clean_text(raw_wb$membrane)
raw_wb$condition <- clean_text(raw_wb$condition)
raw_wb$cell_line <- clean_text(raw_wb$cell_line)
raw_wb$band_value <- as.numeric(raw_wb$band_value)
raw_wb$treatment <- condition_to_treatment(raw_wb$condition)
raw_wb$sample_id <- paste(raw_wb$cell_line, raw_wb$treatment, sep = "_")
raw_wb$batch <- ifelse(grepl("^BN118", raw_wb$cell_line), "BN118", "BN91")

ponceau <- raw_wb[raw_wb$target == "Ponceau", c(
  "cell_line", "condition", "membrane", "band_value"
), drop = FALSE]
names(ponceau)[names(ponceau) == "band_value"] <- "ponceau_value"
ponceau_key <- paste(
  ponceau$cell_line,
  ponceau$condition,
  ponceau$membrane,
  sep = "\r"
)
if (anyDuplicated(ponceau_key)) {
  stop("Ponceau matching keys are duplicated.")
}

bands <- raw_wb[raw_wb$target != "Ponceau", , drop = FALSE]
bands <- merge(
  bands,
  ponceau,
  by = c("cell_line", "condition", "membrane"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(bands$ponceau_value)) {
  stop("At least one band could not be matched to same-membrane Ponceau.")
}
bands$ponceau_normalized <- bands$band_value / bands$ponceau_value
bands$log2_ponceau_normalized <- log2(bands$ponceau_normalized)
bands <- bands[order(
  match(bands$cell_line, cell_order),
  match(bands$treatment, treatment_order),
  bands$membrane,
  bands$target
), , drop = FALSE]
write.csv(
  bands,
  file.path(out_dir, "original_wb_ponceau_normalized_long.csv"),
  row.names = FALSE
)

get_band <- function(target_name) {
  selected <- bands[bands$target == target_name, c(
    "sample_id", "cell_line", "treatment", "batch", "membrane",
    "band_value", "ponceau_value", "ponceau_normalized",
    "log2_ponceau_normalized"
  ), drop = FALSE]
  suffix <- gsub("[^A-Za-z0-9]+", "_", target_name)
  names(selected)[-(1:4)] <- paste0(names(selected)[-(1:4)], "_", suffix)
  selected
}

pmtor <- get_band("pMTOR")
total_mtor <- get_band("Total mTOR")
ps6 <- get_band("pS6")
total_s6 <- get_band("Total S6")

wb_wide <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = c("sample_id", "cell_line", "treatment", "batch"),
    all = TRUE,
    sort = FALSE
  ),
  list(pmtor, total_mtor, ps6, total_s6)
)
if (anyNA(wb_wide)) {
  stop("A required MTOR/S6 band or Ponceau value is missing.")
}

make_readout <- function(name, unlogged, log2_value) {
  data.frame(
    wb_wide[, c("sample_id", "cell_line", "treatment", "batch")],
    wb_readout = name,
    wb_unlogged = unlogged,
    wb_log2 = log2_value,
    stringsAsFactors = FALSE
  )
}

pmtor_ratio <- (
  wb_wide$ponceau_normalized_pMTOR /
    wb_wide$ponceau_normalized_Total_mTOR
)
ps6_ratio <- (
  wb_wide$ponceau_normalized_pS6 /
    wb_wide$ponceau_normalized_Total_S6
)
wb_readouts <- rbind(
  make_readout(
    "pMTOR/Ponceau",
    wb_wide$ponceau_normalized_pMTOR,
    wb_wide$log2_ponceau_normalized_pMTOR
  ),
  make_readout(
    "pMTOR/Total mTOR",
    pmtor_ratio,
    log2(pmtor_ratio)
  ),
  make_readout(
    "pS6/Ponceau",
    wb_wide$ponceau_normalized_pS6,
    wb_wide$log2_ponceau_normalized_pS6
  ),
  make_readout(
    "pS6/Total S6",
    ps6_ratio,
    log2(ps6_ratio)
  )
)
write.csv(
  wb_readouts,
  file.path(out_dir, "original_wb_mtor_s6_readouts.csv"),
  row.names = FALSE
)

sample_metadata <- unique(raw_wb[, c(
  "sample_id", "cell_line", "treatment", "batch"
)])
sample_metadata <- sample_metadata[order(
  match(sample_metadata$cell_line, cell_order),
  match(sample_metadata$treatment, treatment_order)
), , drop = FALSE]
sample_ids <- sample_metadata$sample_id

read_vmpa_scaling <- function(scaling) {
  input_file <- file.path(
    vmpa_dir,
    paste0(
      "vmpa_glioma_n250_",
      scaling,
      "_unique_FALSE_raw_result.csv"
    )
  )
  raw <- read.csv(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  raw$activity_target <- make.unique(as.character(raw$target))
  raw <- raw[toupper(raw$target) %in% vmpa_targets, , drop = FALSE]
  if (!all(sample_ids %in% names(raw))) {
    stop("VMPA sample columns are missing from: ", input_file)
  }

  rows <- lapply(seq_len(nrow(raw)), function(i) {
    data.frame(
      sample_id = sample_ids,
      target = toupper(raw$target[i]),
      activity_target = raw$activity_target[i],
      confidence = raw$conf[i],
      signature = raw$signature[i],
      vmpa_normalization = scaling,
      vmpa_value = as.numeric(raw[i, sample_ids]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

vmpa <- do.call(rbind, lapply(scalings, read_vmpa_scaling))
vmpa <- merge(
  sample_metadata,
  vmpa,
  by = "sample_id",
  all.y = TRUE,
  suffixes = c("", "_vmpa"),
  sort = FALSE
)
if (any(vmpa$cell_line != vmpa$cell_line_vmpa) ||
    any(vmpa$treatment != vmpa$treatment_vmpa) ||
    any(vmpa$batch != vmpa$batch_vmpa)) {
  stop("VMPA and WB sample metadata do not agree.")
}
vmpa <- vmpa[, !grepl("_vmpa$", names(vmpa)), drop = FALSE]
write.csv(
  vmpa,
  file.path(out_dir, "vmpa_mtor_s6_candidates_all_scalings.csv"),
  row.names = FALSE
)

readout_map <- list(
  MTOR = c(
    "pMTOR/Ponceau",
    "pMTOR/Total mTOR",
    "pS6/Ponceau",
    "pS6/Total S6"
  ),
  RPS6 = c("pS6/Ponceau", "pS6/Total S6"),
  RPS6KB2 = c("pS6/Ponceau", "pS6/Total S6"),
  RPS6KA1 = c("pS6/Ponceau", "pS6/Total S6")
)

matched_rows <- lapply(seq_len(nrow(vmpa)), function(i) {
  allowed <- readout_map[[vmpa$target[i]]]
  readout_rows <- wb_readouts[
    wb_readouts$sample_id == vmpa$sample_id[i] &
      wb_readouts$wb_readout %in% allowed,
    ,
    drop = FALSE
  ]
  cbind(
    vmpa[rep(i, nrow(readout_rows)), , drop = FALSE],
    readout_rows[, c("wb_readout", "wb_unlogged", "wb_log2")],
    stringsAsFactors = FALSE
  )
})
matched <- do.call(rbind, matched_rows)
rownames(matched) <- NULL

transform_key <- interaction(
  matched$target,
  matched$activity_target,
  matched$vmpa_normalization,
  matched$wb_readout,
  drop = TRUE,
  lex.order = TRUE
)
transformed <- do.call(rbind, lapply(split(matched, transform_key), function(x) {
  x <- x[order(
    match(x$cell_line, cell_order),
    match(x$treatment, treatment_order)
  ), , drop = FALSE]
  x$vmpa_centered <- ave(
    x$vmpa_value,
    x$cell_line,
    FUN = function(z) z - mean(z)
  )
  x$vmpa_within_cell_z <- ave(x$vmpa_value, x$cell_line, FUN = z_score)
  x$wb_unlogged_centered <- ave(
    x$wb_unlogged,
    x$cell_line,
    FUN = function(z) z - mean(z)
  )
  x$wb_unlogged_within_cell_z <- ave(
    x$wb_unlogged,
    x$cell_line,
    FUN = z_score
  )
  x$wb_log2_centered <- ave(
    x$wb_log2,
    x$cell_line,
    FUN = function(z) z - mean(z)
  )
  x$wb_log2_within_cell_z <- ave(
    x$wb_log2,
    x$cell_line,
    FUN = z_score
  )
  x
}))
rownames(transformed) <- NULL
transform_key <- interaction(
  transformed$target,
  transformed$activity_target,
  transformed$vmpa_normalization,
  transformed$wb_readout,
  drop = TRUE,
  lex.order = TRUE
)
write.csv(
  transformed,
  file.path(out_dir, "matched_transformed_values_all_candidates.csv"),
  row.names = FALSE
)

cohorts <- list(
  BN91_pair_n6 = c("BN91", "BN91R"),
  BN118_pair_n6 = c("BN118", "BN118R"),
  all_four_n12 = cell_order
)
analysis_definitions <- data.frame(
  method = c("raw", "raw", "within_cell_centered", "within_cell_centered",
             "within_cell_z", "within_cell_z"),
  wb_transformation = rep(c("unlogged", "log2"), 3L),
  vmpa_column = c(
    "vmpa_value", "vmpa_value",
    "vmpa_centered", "vmpa_centered",
    "vmpa_within_cell_z", "vmpa_within_cell_z"
  ),
  wb_column = c(
    "wb_unlogged", "wb_log2",
    "wb_unlogged_centered", "wb_log2_centered",
    "wb_unlogged_within_cell_z", "wb_log2_within_cell_z"
  ),
  stringsAsFactors = FALSE
)

stats_rows <- list()
row_index <- 0L
for (group in split(transformed, transform_key)) {
  for (cohort_name in names(cohorts)) {
    cohort_data <- group[group$cell_line %in% cohorts[[cohort_name]], , drop = FALSE]
    for (definition_index in seq_len(nrow(analysis_definitions))) {
      definition <- analysis_definitions[definition_index, , drop = FALSE]
      x <- cohort_data[[definition$vmpa_column]]
      y <- cohort_data[[definition$wb_column]]
      result <- pearson_values(x, y)
      permutation_p <- if (definition$method == "raw") {
        NA_real_
      } else {
        blocked_permutation_p(x, y, cohort_data$cell_line)
      }
      row_index <- row_index + 1L
      stats_rows[[row_index]] <- data.frame(
        target = cohort_data$target[1],
        activity_target = cohort_data$activity_target[1],
        confidence = cohort_data$confidence[1],
        signature = cohort_data$signature[1],
        vmpa_normalization = cohort_data$vmpa_normalization[1],
        wb_readout = cohort_data$wb_readout[1],
        cohort = cohort_name,
        method = definition$method,
        wb_transformation = definition$wb_transformation,
        n = unname(result["n"]),
        pearson_r = unname(result["r"]),
        naive_pearson_p = unname(result["p"]),
        blocked_exact_permutation_p = permutation_p,
        stringsAsFactors = FALSE
      )
    }
  }
}
stats <- do.call(rbind, stats_rows)
stats <- stats[order(
  match(stats$target, vmpa_targets),
  stats$activity_target,
  match(stats$vmpa_normalization, scalings),
  stats$wb_readout,
  match(stats$cohort, names(cohorts)),
  match(stats$method, c("raw", "within_cell_centered", "within_cell_z")),
  match(stats$wb_transformation, c("unlogged", "log2"))
), , drop = FALSE]
write.csv(
  stats,
  file.path(out_dir, "correlations_mtor_s6_n6_n12_all_scalings.csv"),
  row.names = FALSE
)

primary <- stats[
  stats$vmpa_normalization == "signature_z" &
    stats$method == "within_cell_centered" &
    stats$wb_transformation == "log2",
  ,
  drop = FALSE
]
write.csv(
  primary,
  file.path(out_dir, "primary_signature_z_centered_log2_summary.csv"),
  row.names = FALSE
)

format_p <- function(p) {
  if (!is.finite(p)) {
    return("NA")
  }
  if (p < 0.001) {
    return(formatC(p, format = "e", digits = 2))
  }
  formatC(p, format = "f", digits = 3)
}

cell_colors <- c(
  BN91 = "#0072B2",
  BN91R = "#D55E00",
  BN118 = "#009E73",
  BN118R = "#CC79A7"
)
treatment_shapes <- c(DMSO = 16, `50nM` = 17, `500nM` = 15)

make_mtor_plot <- function(cohort_name) {
  selected_cells <- cohorts[[cohort_name]]
  plot_data <- transformed[
    transformed$target == "MTOR" &
      transformed$vmpa_normalization == "signature_z" &
      transformed$cell_line %in% selected_cells,
    ,
    drop = FALSE
  ]
  annotation <- primary[
    primary$target == "MTOR" &
      primary$cohort == cohort_name,
    ,
    drop = FALSE
  ]
  annotation$label <- paste0(
    "r = ", formatC(annotation$pearson_r, format = "f", digits = 2),
    "\nblocked p = ",
    vapply(annotation$blocked_exact_permutation_p, format_p, character(1))
  )

  plot <- ggplot(
    plot_data,
    aes(
      x = vmpa_centered,
      y = wb_log2_centered,
      color = cell_line,
      shape = treatment
    )
  ) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.45) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.45) +
    geom_smooth(
      data = plot_data,
      aes(x = vmpa_centered, y = wb_log2_centered, group = wb_readout),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_point(size = 4.3, stroke = 0.8) +
    geom_label(
      data = annotation,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.08,
      size = 4.2,
      lineheight = 0.95,
      fill = "white",
      linewidth = 0
    ) +
    facet_wrap(~wb_readout, nrow = 1, scales = "free_y") +
    scale_color_manual(values = cell_colors, name = "Cell line") +
    scale_shape_manual(values = treatment_shapes, name = "Treatment") +
    labs(
      title = paste0("MTOR activity versus WB | ", cohort_name),
      subtitle = paste(
        "Original WB labels; within-cell-line centered;",
        "signature_z; log2 WB readouts"
      ),
      x = "MTOR VMPA activity, centered within cell line",
      y = "WB readout, centered within cell line"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(size = 19, face = "bold"),
      plot.subtitle = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.spacing = grid::unit(0.9, "lines")
    )

  stub <- paste0("mtor_signature_z_centered_", cohort_name)
  ggsave(
    file.path(out_dir, paste0(stub, ".pdf")),
    plot,
    width = 15,
    height = 5.6,
    device = "pdf"
  )
  ggsave(
    file.path(out_dir, paste0(stub, ".png")),
    plot,
    width = 15,
    height = 5.6,
    dpi = 300
  )
}

make_ps6_candidate_plot <- function() {
  plot_data <- transformed[
    transformed$target %in% c("RPS6", "RPS6KB2", "RPS6KA1") &
      transformed$vmpa_normalization == "signature_z" &
      transformed$wb_readout == "pS6/Total S6",
    ,
    drop = FALSE
  ]
  plot_data$panel <- paste0(
    plot_data$activity_target,
    " (",
    plot_data$confidence,
    "): ",
    sub("^.*:", "", plot_data$signature)
  )
  panel_levels <- unique(plot_data$panel)
  plot_data$panel <- factor(plot_data$panel, levels = panel_levels)

  annotation <- primary[
    primary$target %in% c("RPS6", "RPS6KB2", "RPS6KA1") &
      primary$wb_readout == "pS6/Total S6" &
      primary$cohort == "all_four_n12",
    ,
    drop = FALSE
  ]
  annotation$panel <- paste0(
    annotation$activity_target,
    " (",
    annotation$confidence,
    "): ",
    sub("^.*:", "", annotation$signature)
  )
  annotation$panel <- factor(annotation$panel, levels = panel_levels)
  annotation$label <- paste0(
    "r = ", formatC(annotation$pearson_r, format = "f", digits = 2),
    "\nblocked p = ",
    vapply(annotation$blocked_exact_permutation_p, format_p, character(1))
  )

  plot <- ggplot(
    plot_data,
    aes(
      x = vmpa_centered,
      y = wb_log2_centered,
      color = cell_line,
      shape = treatment
    )
  ) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.45) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.45) +
    geom_smooth(
      data = plot_data,
      aes(x = vmpa_centered, y = wb_log2_centered, group = panel),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_point(size = 4.1, stroke = 0.8) +
    geom_label(
      data = annotation,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.08,
      size = 3.8,
      lineheight = 0.95,
      fill = "white",
      linewidth = 0
    ) +
    facet_wrap(~panel, nrow = 2, scales = "free_x") +
    scale_color_manual(values = cell_colors, name = "Cell line") +
    scale_shape_manual(values = treatment_shapes, name = "Treatment") +
    labs(
      title = "Exploratory pS6-associated VMPA signatures | all four, n=12",
      subtitle = paste(
        "Original WB labels; within-cell-line centered;",
        "signature_z; log2 pS6/Total S6"
      ),
      x = "VMPA activity, centered within cell line",
      y = "Log2 pS6/Total S6, centered within cell line"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(size = 19, face = "bold"),
      plot.subtitle = element_text(size = 13),
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.spacing = grid::unit(0.9, "lines")
    )

  stub <- "ps6_candidate_signatures_signature_z_centered_all_four_n12"
  ggsave(
    file.path(out_dir, paste0(stub, ".pdf")),
    plot,
    width = 14,
    height = 9,
    device = "pdf"
  )
  ggsave(
    file.path(out_dir, paste0(stub, ".png")),
    plot,
    width = 14,
    height = 9,
    dpi = 300
  )
}

# Exploratory multi-panel figures were intentionally retired. The separate
# publication plotting script generates only the selected final PDFs.

message("Wrote original-label MTOR/S6 analysis to: ", normalizePath(out_dir))
