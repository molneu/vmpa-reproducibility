#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
})

score_file <- here::here(
  "reproducibility", "figure6", "data", "wb_correlation",
  "vmpa_wb_correlations",
  "vmpa_glioma_n250_sample_pop_sd_unique_FALSE_raw_result.csv"
)
out_dir <- here::here(
  "reproducibility", "figure6", "results", "apoptosis_pathway"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(score_file)) {
  stop("Missing VMPA sample-population-SD score file: ", score_file, call. = FALSE)
}

scores <- read.csv(score_file, check.names = FALSE, stringsAsFactors = FALSE)
scores$activity_target <- make.unique(as.character(scores$target))

selected <- data.frame(
  node = c("BCL-xL", "PUMA", "BAX"),
  target = c("BCL2L1", "BBC3", "BAX"),
  stringsAsFactors = FALSE
)

select_rows <- function(target) {
  rows <- which(scores$target == target)
  if (target %in% c("BCL2L1", "BBC3")) {
    confidence_three <- rows[scores$conf[rows] == "c3"]
    if (length(confidence_three) > 0L) rows <- confidence_three
  }
  if (length(rows) == 0L) {
    stop("No VMPA signature found for ", target, call. = FALSE)
  }
  rows
}

make_panel <- function(cell_line) {
  treatments <- c("DMSO", "50nM", "500nM")
  values <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
    rows <- select_rows(selected$target[i])
    data.frame(
      node = selected$node[i],
      target = selected$target[i],
      treatment = treatments,
      score = vapply(
        treatments,
        function(treatment) {
          mean(scores[rows, paste0(cell_line, "_", treatment)])
        },
        numeric(1)
      ),
      n_signatures = length(rows),
      signatures = paste(scores$signature[rows], collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(
    values,
    file.path(out_dir, paste0("Figure6j_", cell_line, "_scores.csv")),
    row.names = FALSE
  )

  values$x <- c(`BCL-xL` = 1, PUMA = 1, BAX = 2)[values$node]
  values$y <- c(`BCL-xL` = 2, PUMA = 1, BAX = 1.5)[values$node]
  values$treatment <- factor(
    values$treatment,
    levels = treatments,
    labels = c("0", "50", "500")
  )

  plot <- ggplot(values, aes(x, y)) +
    geom_segment(
      data = data.frame(x = 1.35, y = 1.5, xend = 1.65, yend = 1.5),
      aes(x, y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(1.5, "mm")),
      linewidth = 0.4
    ) +
    geom_label(
      aes(label = node, fill = score),
      size = 2.2,
      linewidth = 0.25,
      label.padding = grid::unit(1.2, "mm")
    ) +
    facet_grid(treatment ~ .) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      name = "Activity"
    ) +
    coord_cartesian(xlim = c(0.6, 2.4), ylim = c(0.55, 2.4), expand = FALSE) +
    labs(title = cell_line, x = NULL, y = NULL) +
    theme_void(base_size = 7) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 7),
      strip.text.y = element_text(size = 6),
      legend.position = "bottom",
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5)
    )

  ggsave(
    file.path(out_dir, paste0("Figure6j_", cell_line, ".pdf")),
    plot,
    width = 48,
    height = 75,
    units = "mm"
  )
}

make_panel("BN118R")
make_panel("BN91R")
