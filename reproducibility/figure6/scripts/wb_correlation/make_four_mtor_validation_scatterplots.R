#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
})

analysis_dir <- file.path(
  here::here("reproducibility", "figure6", "results", "wb_correlation"),
  "original_wb_labels_mtor_s6_n6_n12"
)
data_file <- file.path(
  analysis_dir,
  "matched_transformed_values_all_candidates.csv"
)
stats_file <- file.path(
  analysis_dir,
  "primary_signature_z_centered_log2_summary.csv"
)
out_dir <- file.path(analysis_dir, "SELECTED_FINAL_PLOTS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(
  data_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stats <- read.csv(
  stats_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

make_definition <- function(
    target,
    activity_target,
    cohort,
    wb_readout,
    title,
    file_stub) {
  data.frame(
    target = target,
    activity_target = activity_target,
    cohort = cohort,
    wb_readout = wb_readout,
    title = title,
    file_stub = file_stub,
    stringsAsFactors = FALSE
  )
}

panel_definitions <- do.call(rbind, list(
  make_definition(
    "MTOR", "MTOR", "BN91_pair_n6", "pMTOR/Ponceau",
    "BN91 | pMTOR", "BN91_MTOR_vs_pMTOR_Ponceau"
  ),
  make_definition(
    "MTOR", "MTOR", "BN118_pair_n6", "pMTOR/Ponceau",
    "BN118 | pMTOR", "BN118_MTOR_vs_pMTOR_Ponceau"
  ),
  make_definition(
    "MTOR", "MTOR", "BN91_pair_n6", "pS6/Ponceau",
    "BN91 | pS6", "BN91_MTOR_vs_pS6_Ponceau"
  ),
  make_definition(
    "MTOR", "MTOR", "BN118_pair_n6", "pS6/Ponceau",
    "BN118 | pS6", "BN118_MTOR_vs_pS6_Ponceau"
  ),
  make_definition(
    "MTOR", "MTOR", "BN91_pair_n6", "pMTOR/Total mTOR",
    "BN91 | pMTOR/mTOR", "BN91_MTOR_vs_pMTOR_Total_mTOR"
  ),
  make_definition(
    "MTOR", "MTOR", "BN118_pair_n6", "pMTOR/Total mTOR",
    "BN118 | pMTOR/mTOR", "BN118_MTOR_vs_pMTOR_Total_mTOR"
  ),
  make_definition(
    "MTOR", "MTOR", "BN91_pair_n6", "pS6/Total S6",
    "BN91 | pS6/S6", "BN91_MTOR_vs_pS6_Total_S6"
  ),
  make_definition(
    "MTOR", "MTOR", "BN118_pair_n6", "pS6/Total S6",
    "BN118 | pS6/S6", "BN118_MTOR_vs_pS6_Total_S6"
  ),
  make_definition(
    "MTOR", "MTOR", "all_four_n12", "pMTOR/Ponceau",
    "All 12 | pMTOR", "all12_MTOR_vs_pMTOR_Ponceau"
  ),
  make_definition(
    "MTOR", "MTOR", "all_four_n12", "pS6/Ponceau",
    "All 12 | pS6", "all12_MTOR_vs_pS6_Ponceau"
  ),
  make_definition(
    "MTOR", "MTOR", "all_four_n12", "pMTOR/Total mTOR",
    "All 12 | pMTOR/mTOR", "all12_MTOR_vs_pMTOR_Total_mTOR"
  ),
  make_definition(
    "MTOR", "MTOR", "all_four_n12", "pS6/Total S6",
    "All 12 | pS6/S6", "all12_MTOR_vs_pS6_Total_S6"
  )
))

rps6_signatures <- data.frame(
  activity_target = c("RPS6", "RPS6.1", "RPS6.2"),
  signature_short = c("K09", "O13", "I22"),
  stringsAsFactors = FALSE
)
for (signature_index in seq_len(nrow(rps6_signatures))) {
  for (cohort in c("BN91_pair_n6", "BN118_pair_n6", "all_four_n12")) {
    family_label <- switch(
      cohort,
      BN91_pair_n6 = "BN91",
      BN118_pair_n6 = "BN118",
      all_four_n12 = "All 12"
    )
    file_family <- switch(
      cohort,
      BN91_pair_n6 = "BN91",
      BN118_pair_n6 = "BN118",
      all_four_n12 = "all12"
    )
    for (readout in c("pS6/Ponceau", "pS6/Total S6")) {
      readout_short <- if (readout == "pS6/Ponceau") "pS6" else "pS6/S6"
      panel_definitions <- rbind(
        panel_definitions,
        make_definition(
          "RPS6",
          rps6_signatures$activity_target[signature_index],
          cohort,
          readout,
          paste0(
            family_label,
            " | ",
            rps6_signatures$signature_short[signature_index],
            "-",
            readout_short
          ),
          paste(
            file_family,
            rps6_signatures$activity_target[signature_index],
            "vs",
            gsub("[^A-Za-z0-9]+", "_", readout),
            sep = "_"
          )
        )
      )
    }
  }
}

selected_file_stubs <- c(
  "all12_MTOR_vs_pMTOR_Ponceau",
  "BN91_MTOR_vs_pMTOR_Ponceau",
  "BN118_MTOR_vs_pMTOR_Ponceau",
  "BN91_MTOR_vs_pS6_Ponceau",
  "BN118_MTOR_vs_pS6_Ponceau",
  "BN91_RPS6.1_vs_pS6_Ponceau",
  "BN118_RPS6.1_vs_pS6_Ponceau"
)
panel_definitions <- panel_definitions[
  match(selected_file_stubs, panel_definitions$file_stub),
  ,
  drop = FALSE
]

cohort_cells <- list(
  BN91_pair_n6 = c("BN91", "BN91R"),
  BN118_pair_n6 = c("BN118", "BN118R"),
  all_four_n12 = c("BN91", "BN91R", "BN118", "BN118R")
)
treatment_shapes <- c(DMSO = 16, `50nM` = 16, `500nM` = 16)

main_figure_mask <- (
  panel_definitions$cohort %in% c("BN91_pair_n6", "BN118_pair_n6")
)
main_definitions <- panel_definitions[main_figure_mask, , drop = FALSE]
main_regression_rows <- lapply(seq_len(nrow(main_definitions)), function(i) {
  definition <- main_definitions[i, , drop = FALSE]
  selected_cells <- cohort_cells[[definition$cohort]]
  plot_data <- data[
    data$target == definition$target &
      data$activity_target == definition$activity_target &
      data$vmpa_normalization == "signature_z" &
      data$wb_readout == definition$wb_readout &
      data$cell_line %in% selected_cells,
    ,
    drop = FALSE
  ]
  pearson_r <- cor(
    plot_data$vmpa_centered,
    plot_data$wb_log2_centered
  )
  model <- lm(
    wb_log2_centered ~ vmpa_centered,
    data = plot_data
  )
  model_summary <- summary(model)
  slope_t <- model_summary$coefficients["vmpa_centered", "t value"]
  data.frame(
    file_stub = definition$file_stub,
    target = definition$target,
    activity_target = definition$activity_target,
    wb_readout = definition$wb_readout,
    cohort = definition$cohort,
    n = nrow(plot_data),
    pearson_r = pearson_r,
    r_squared = model_summary$r.squared,
    one_sided_regression_p = pt(
      slope_t,
      df.residual(model),
      lower.tail = FALSE
    ),
    stringsAsFactors = FALSE
  )
})
main_regressions <- do.call(rbind, main_regression_rows)
main_regressions$fdr_adjusted_p <- p.adjust(
  main_regressions$one_sided_regression_p,
  method = "fdr"
)

format_p <- function(p) {
  if (p < 0.001) {
    formatC(p, format = "e", digits = 1)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

make_panel <- function(definition, compact = FALSE, compact_mm = 14) {
  selected_cells <- cohort_cells[[definition$cohort]]
  plot_data <- data[
    data$target == definition$target &
      data$activity_target == definition$activity_target &
      data$vmpa_normalization == "signature_z" &
      data$wb_readout == definition$wb_readout &
      data$cell_line %in% selected_cells,
    ,
    drop = FALSE
  ]
  panel_stats <- stats[
    stats$target == definition$target &
      stats$activity_target == definition$activity_target &
      stats$cohort == definition$cohort &
      stats$wb_readout == definition$wb_readout,
    ,
    drop = FALSE
  ]
  expected_n <- if (definition$cohort == "all_four_n12") 12L else 6L
  if (nrow(plot_data) != expected_n || nrow(panel_stats) != 1L) {
    stop("Unexpected data size for panel: ", definition$file_stub)
  }

  main_regression <- main_regressions[
    main_regressions$file_stub == definition$file_stub,
    ,
    drop = FALSE
  ]
  panel_model <- lm(
    wb_log2_centered ~ vmpa_centered,
    data = plot_data
  )
  panel_model_summary <- summary(panel_model)
  if (nrow(main_regression) == 1L) {
    annotation <- ""
    display_title <- paste0(
      definition$title,
      "\nR2=",
      formatC(main_regression$r_squared, format = "f", digits = 2),
      "; FDR p=",
      format_p(main_regression$fdr_adjusted_p)
    )
  } else {
    annotation <- paste0(
      "R2=",
      formatC(panel_model_summary$r.squared, format = "f", digits = 2),
      "\np=",
      format_p(panel_stats$blocked_exact_permutation_p)
    )
    display_title <- definition$title
  }

  if (compact) {
    base_size <- 4.5
    title_size <- 4.8
    axis_title_size <- 4.3
    axis_text_size <- 3.8
    annotation_size <- 1.3
    point_size <- if (compact_mm >= 18) 1.3 else 1.15
    line_width <- 0.32
    x_label <- "VMPA"
    y_label <- "log2 WB"
    plot_margin <- grid::unit(c(0.6, 0.6, 0.5, 0.6), "mm")
  } else {
    base_size <- 9
    title_size <- 7
    axis_title_size <- 9
    axis_text_size <- 8
    annotation_size <- 2.8
    point_size <- 3.2
    line_width <- 0.65
    x_label <- paste0(
      definition$activity_target,
      " VMPA activity\n(centered within cell line)"
    )
    y_label <- paste0(
      "Log2 ",
      definition$wb_readout,
      "\n(centered within cell line)"
    )
    plot_margin <- grid::unit(c(2, 2, 1.5, 2), "mm")
  }

  ggplot(
    plot_data,
    aes(
      x = vmpa_centered,
      y = wb_log2_centered
    )
  ) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
    geom_smooth(
      data = plot_data,
      aes(x = vmpa_centered, y = wb_log2_centered),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      level = 0.95,
      color = "black",
      fill = "gray",
      alpha = 0.6,
      linewidth = 0.5
    ) +
    geom_point(color = "darkred", shape = 16, size = 2) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.08,
      vjust = 1.08,
      size = annotation_size,
      lineheight = 0.9,
      color = "black"
    ) +
    scale_x_continuous(
      breaks = function(x) pretty(x, n = if (compact) 2 else 4),
      expand = expansion(mult = c(0.12, 0.10))
    ) +
    scale_y_continuous(
      breaks = function(x) pretty(x, n = if (compact) 2 else 4),
      expand = expansion(mult = c(0.12, 0.12))
    ) +
    labs(
      title = display_title,
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = grid::unit(2.5, "mm"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.margin = plot_margin
    )
}

for (i in seq_len(nrow(panel_definitions))) {
  definition <- panel_definitions[i, , drop = FALSE]

  readable_plot <- make_panel(definition, compact = FALSE)
  ggsave(
    file.path(out_dir, paste0(definition$file_stub, "_50mm.pdf")),
    readable_plot,
    width = 50,
    height = 50,
    units = "mm",
    device = "pdf"
  )
}

write.csv(
  panel_definitions,
  file.path(out_dir, "panel_definitions.csv"),
  row.names = FALSE
)
write.csv(
  main_regressions,
  file.path(out_dir, "supplementary_six_regression_results.csv"),
  row.names = FALSE
)

message("Wrote selected MTOR validation panels to: ", normalizePath(out_dir))
