# Pathway-style schematic for intrinsic apoptosis COMPASS scores.
#
# Biological models: BN91, BN91R, BN118, BN118R.
# Nodes are colored by absolute sample-wise population SD normalized activity
# scores for DMSO, Temsi50, and Temsi500. MCL1 is from SigCom LINCS
# bottom-300 because it is not available in glioma COMPASS/protivity.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(script_path) == 1 && file.exists(script_path)) {
  dirname(normalizePath(script_path))
} else {
  getwd()
}

find_repo_root <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "reproducibility"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find repository root from: ", script_dir)
    }
    path <- parent
  }
}

repo_root <- find_repo_root(script_dir)
data_dir <- file.path(repo_root, "reproducibility", "figure6", "data", "apoptosis_pathway")
res_dir <- file.path(repo_root, "reproducibility", "figure6", "results", "apoptosis_pathway")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

models_to_plot <- c("BN91", "BN91R", "BN118", "BN118R")

clean_sample_name <- function(x) {
  x <- sub("^[0-9]+_", "", x)
  sub("_S[0-9]+_L001[.]dedup$", "", x)
}

sample_metadata <- function(samples) {
  sample_clean <- clean_sample_name(samples)
  treatment <- rep(NA_character_, length(sample_clean))
  treatment[grepl("DMSO", sample_clean)] <- "DMSO"
  treatment[grepl("Temsi_50", sample_clean)] <- "Temsi50"
  treatment[grepl("Temsi_500", sample_clean)] <- "Temsi500"

  annotated_cell_line <- ifelse(grepl("BN91", sample_clean), "BN91", "BN118")
  recurrence <- ifelse(grepl("_R_", sample_clean), "R", "P")
  annotated_model <- paste(annotated_cell_line, recurrence, sep = "_")
  annotated_label <- c(BN91_P = "BN91", BN91_R = "BN91R", BN118_P = "BN118", BN118_R = "BN118R")[annotated_model]
  biological_label <- c(BN91 = "BN118", BN91R = "BN118R", BN118 = "BN91", BN118R = "BN91R")[annotated_label]

  data.frame(
    sample = samples,
    sample_clean = sample_clean,
    annotated_label = unname(annotated_label),
    biological_label = unname(biological_label),
    treatment = factor(treatment, levels = c("DMSO", "Temsi50", "Temsi500")),
    stringsAsFactors = FALSE
  )
}

node_info <- data.frame(
  node = c("Temsirolimus", "PUMA", "Bcl-xL", "MCL1", "BAX", "MOMP /\ncyto c release", "Caspase-3", "Apoptosis"),
  protein = c(NA, "BBC3", "BCL2L1", "MCL1", "BAX", NA, "CASP3", NA),
  source = c(NA, "COMPASS", "COMPASS", "SigCom", "COMPASS", NA, "COMPASS", NA),
  x = c(0.55, 2.05, 2.05, 3.7, 3.7, 5.25, 6.95, 8.55),
  y = c(2.65, 2.65, 1.25, 1.25, 2.65, 2.65, 2.65, 2.65),
  node_type = c("input", "protein", "protein", "protein", "protein", "process", "protein", "outcome"),
  stringsAsFactors = FALSE
)

activation_edges <- data.frame(
  x = c(0.95, 2.5, 4.15, 5.75, 7.4),
  y = c(2.65, 2.65, 2.65, 2.65, 2.65),
  xend = c(1.6, 3.25, 4.75, 6.45, 8.1),
  yend = c(2.65, 2.65, 2.65, 2.65, 2.65),
  label = c("induces", "activates", "MOMP", "activates", "executes"),
  stringsAsFactors = FALSE
)

inhibition_edges <- data.frame(
  x = c(2.05, 2.45, 3.7),
  y = c(1.62, 1.45, 1.62),
  xend = c(2.05, 3.38, 3.7),
  yend = c(2.28, 2.42, 2.28),
  target = c("PUMA", "BAX", "BAX"),
  stringsAsFactors = FALSE
)

make_inhibition_bars <- function(edges, bar_width = 0.34) {
  out <- edges
  dx <- out$xend - out$x
  dy <- out$yend - out$y
  edge_len <- sqrt(dx^2 + dy^2)
  px <- -dy / edge_len
  py <- dx / edge_len
  out$bar_x <- out$xend - px * bar_width / 2
  out$bar_y <- out$yend - py * bar_width / 2
  out$bar_xend <- out$xend + px * bar_width / 2
  out$bar_yend <- out$yend + py * bar_width / 2
  out
}

inhibition_bars <- make_inhibition_bars(inhibition_edges)

gsva_data <- read.csv(
  file.path(data_dir, "COMPASS_glioma_protivity_GSVA_scores_sample_wise_population_sd.csv"),
  row.names = 1,
  check.names = FALSE
)
gsva_data <- as.matrix(gsva_data)
storage.mode(gsva_data) <- "numeric"
meta <- sample_metadata(colnames(gsva_data))

sigcom_file <- file.path(data_dir, "SigComLINCS_CRISPR_bottom300_GSVA_scores_sample_wise_population_sd.csv")
if (!file.exists(sigcom_file)) {
  stop("Missing SigCom LINCS bottom-300 score matrix: ", sigcom_file)
}
sigcom_data <- read.csv(sigcom_file, row.names = 1, check.names = FALSE)
sigcom_data <- as.matrix(sigcom_data)
storage.mode(sigcom_data) <- "numeric"

score_for <- function(protein, sample) {
  if (protein == "MCL1") {
    if (!protein %in% rownames(sigcom_data)) return(NA_real_)
    return(sigcom_data[protein, sample])
  }
  sigs <- rownames(gsva_data)[sub("_.*", "", rownames(gsva_data)) == protein]
  if (protein %in% c("BBC3", "BCL2L1", "CASP3")) {
    c3 <- sigs[grepl("_c3$", sigs)]
    if (length(c3) > 0) sigs <- c3
  }
  if (length(sigs) == 0) return(NA_real_)
  mean(gsva_data[sigs, sample], na.rm = TRUE)
}

plot_model <- function(model_label) {
  samples_i <- meta[meta$biological_label == model_label, ]
  if (nrow(samples_i) == 0) {
    warning("No samples found for ", model_label)
    return(invisible(NULL))
  }

  plot_df <- do.call(rbind, lapply(seq_len(nrow(samples_i)), function(i) {
    sample_i <- samples_i$sample[i]
    treatment_i <- as.character(samples_i$treatment[i])
    out <- node_info
    out$treatment <- treatment_i
    out$score <- NA_real_
    protein_idx <- which(!is.na(out$protein))
    out$score[protein_idx] <- vapply(out$protein[protein_idx], score_for, numeric(1), sample = sample_i)
    out
  }))
  plot_df$treatment <- factor(plot_df$treatment, levels = c("DMSO", "Temsi50", "Temsi500"))

  write.csv(
    plot_df,
    file.path(res_dir, paste0("COMPASS_intrinsic_apoptosis_pathway_", model_label, "_scores.csv")),
    row.names = FALSE
  )

  p <- ggplot() +
    geom_segment(
      data = activation_edges,
      aes(x, y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.14, "inches"), type = "closed"),
      linewidth = 0.45,
      color = "grey30"
    ) +
    geom_segment(
      data = inhibition_edges,
      aes(x, y, xend = xend, yend = yend),
      linewidth = 0.45,
      linetype = "dashed",
      color = "grey30"
    ) +
    geom_segment(
      data = inhibition_bars,
      aes(x = bar_x, y = bar_y, xend = bar_xend, yend = bar_yend),
      linewidth = 0.95,
      lineend = "round",
      color = "#8B1A1A"
    ) +
    geom_point(
      data = plot_df[plot_df$node_type == "protein", ],
      aes(x, y, fill = score),
      shape = 21,
      size = 15,
      color = "grey25",
      stroke = 0.5
    ) +
    geom_label(
      data = plot_df[plot_df$node_type != "protein", ],
      aes(x, y, label = node),
      size = 3,
      label.size = 0.3,
      fill = "grey94",
      color = "grey20"
    ) +
    geom_text(
      data = plot_df[plot_df$node_type == "protein", ],
      aes(x, y, label = node),
      size = 3,
      lineheight = 0.9
    ) +
    geom_text(
      data = plot_df[plot_df$node_type == "protein", ],
      aes(x, y - 0.34, label = sprintf("%.2f", score)),
      size = 2.8,
      color = "grey15"
    ) +
    facet_wrap(~ treatment, ncol = 1) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      na.value = "grey90",
      name = "Activity score"
    ) +
    annotate("text", x = 2.95, y = 1.7, label = "blocks", size = 2.7, color = "#8B1A1A") +
    coord_cartesian(xlim = c(0, 9.2), ylim = c(0.6, 3.35), expand = FALSE) +
    theme_void(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0(model_label, ": intrinsic apoptosis pathway activity"),
      subtitle = "PUMA, Bcl-xL, BAX and Caspase-3 are COMPASS scores; MCL1 is SigCom LINCS bottom-300; process/outcome nodes are unmeasured"
    )

  out_pdf <- file.path(res_dir, paste0("COMPASS_intrinsic_apoptosis_pathway_", model_label, ".pdf"))
  ggsave(out_pdf, p, width = 10, height = 8)
  message("Output: ", out_pdf)
}

invisible(lapply(models_to_plot, plot_model))
