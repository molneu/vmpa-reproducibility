#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
})

out_dir <- here::here(
  "reproducibility", "figure6", "results", "target_selection",
  "vmpa_wb_correlations"
)
screen_file <- file.path(
  out_dir,
  "vmpa_initial_unbiased_targetable_screen_n250_sample_pop_sd_positive_3of4.csv"
)
response_file <- file.path(
  out_dir,
  "vmpa_temsi_50_500_average_response_by_setting.csv"
)

if (!all(file.exists(c(screen_file, response_file)))) {
  stop("Required VMPA candidate-screen input is missing.", call. = FALSE)
}

screen <- read.csv(screen_file, check.names = FALSE, stringsAsFactors = FALSE)
response <- read.csv(response_file, check.names = FALSE, stringsAsFactors = FALSE)

screen <- screen[
  screen$unique_signatures == FALSE &
    screen$passes_targetability &
    screen$n_lines_up_both_doses == 4,
  ,
  drop = FALSE
]

response <- response[
  response$gene_set_size == 250 &
    response$scaling == "sample_pop_sd" &
    response$unique_signatures == FALSE &
    response$activity_target %in% screen$target_id,
  ,
  drop = FALSE
]

summarise_target <- function(target_id) {
  x <- response[response$activity_target == target_id, , drop = FALSE]
  data.frame(
    target_id = target_id,
    median_response_50 = median(x$delta_50nM_vs_DMSO),
    median_response_500 = median(x$delta_500nM_vs_DMSO),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, lapply(unique(screen$target_id), summarise_target))
summary_table <- merge(
  summary_table,
  screen[, c(
    "target_id", "gene", "median_mean_delta", "targetability_tier",
    "targetability_basis", "clinical_tractability_status", "screen_rank"
  )],
  by = "target_id",
  all.x = TRUE,
  sort = FALSE
)
summary_table$status <- ifelse(
  summary_table$clinical_tractability_status == "Approved drug",
  "Approved",
  ifelse(
    summary_table$clinical_tractability_status %in% c(
      "Phase I clinical", "Advanced clinical"
    ),
    "Clinical phase",
    "Tractable"
  )
)

make_plot <- function(dose) {
  response_column <- paste0("median_response_", dose)
  x <- summary_table[order(
    -summary_table$median_mean_delta,
    -summary_table[[response_column]],
    -summary_table$targetability_tier,
    summary_table$screen_rank
  ), , drop = FALSE]
  x <- head(x, 10)
  x$plot_value <- x[[response_column]]
  x <- x[order(-x$plot_value), , drop = FALSE]
  x$plot_label <- factor(x$target_id, levels = rev(x$target_id))

  write.csv(
    x,
    file.path(out_dir, paste0("Figure6i_", dose, "nM_plot_data.csv")),
    row.names = FALSE
  )

  plot <- ggplot(
    x,
    aes(
      x = plot_value,
      y = plot_label,
      color = status,
      size = targetability_tier
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_manual(
      values = c(
        Tractable = "#6baed6",
        `Clinical phase` = "#d6604d",
        Approved = "#7f0000"
      )
    ) +
    scale_size_continuous(range = c(2.5, 5), breaks = 3:4, limits = c(3, 4)) +
    labs(
      x = "VMPA activity score",
      y = NULL,
      color = "Status",
      size = "Tractability tier"
    ) +
    theme_classic(base_size = 7) +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5)
    )

  ggsave(
    file.path(out_dir, paste0("Figure6i_temsirolimus_", dose, "nM.pdf")),
    plot,
    width = 42,
    height = 65,
    units = "mm"
  )
}

make_plot("50")
make_plot("500")
