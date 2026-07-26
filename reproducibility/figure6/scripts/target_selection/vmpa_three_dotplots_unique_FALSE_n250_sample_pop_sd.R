#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

out_dir <- here::here("reproducibility", "figure6", "results", "target_selection", "vmpa_wb_correlations")
response_file <- file.path(out_dir, "vmpa_temsi_50_500_average_response_by_setting.csv")

args <- commandArgs(trailingOnly = TRUE)
screen_rule <- if (length(args) > 0) tolower(args[1]) else "both"
if (!screen_rule %in% c("both", "either", "positive")) screen_rule <- "both"
file_tag <- if (screen_rule == "either") "_top5_either" else if (screen_rule == "positive") "_positive_3of4" else ""
display_n <- if (screen_rule == "positive" && length(args) > 1) as.integer(args[2]) else if (screen_rule == "positive") 10 else 12
if (!is.finite(display_n) || display_n < 1) display_n <- 10
display_tag <- if (screen_rule == "positive") paste0("_top", display_n) else ""
size_mode <- if (length(args) > 2 && tolower(args[3]) == "dose") "dose" else "both"
size_tag <- if (size_mode == "dose") "_dose_positive" else ""
selection_mode <- if (length(args) > 3 && tolower(args[4]) == "activity") "activity" else "robustness"
selection_tag <- if (selection_mode == "activity") "_activity_ranked" else ""
screen_file <- file.path(
  out_dir,
  paste0("vmpa_initial_unbiased_targetable_screen_n250_sample_pop_sd", file_tag, ".csv")
)
plot_file <- file.path(
  out_dir,
  paste0("vmpa_initial_unbiased_targetable_plot_data_n250_sample_pop_sd", file_tag, display_tag, ".csv")
)

required_files <- c(screen_file, response_file)
if (!(screen_rule == "positive" && display_n != 10)) required_files <- c(required_files, plot_file)
if (!all(file.exists(required_files))) {
  stop("Required input file is missing.")
}

suppressPackageStartupMessages(library(ggplot2))

screen <- read.csv(screen_file, check.names = FALSE, stringsAsFactors = FALSE)
plot_candidates <- if (file.exists(plot_file)) {
  read.csv(plot_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  NULL
}
response <- read.csv(response_file, check.names = FALSE, stringsAsFactors = FALSE)

# Keep a data-derived candidate set from the averaged-response screen.
# For exploratory top-N plots, select directly from the full screen table.
if (screen_rule == "positive" && display_n != 10) {
  candidate_pool <- screen[
    screen$unique_signatures == FALSE & screen$passes_targetability,
    ,
    drop = FALSE
  ]
  if (selection_mode == "activity") {
    candidate_pool <- candidate_pool[order(
      -candidate_pool$median_mean_delta,
      -candidate_pool$n_lines_up_both_doses,
      candidate_pool$rank
    ), , drop = FALSE]
  } else {
    candidate_pool <- candidate_pool[order(
      -candidate_pool$n_lines_up_both_doses,
      -candidate_pool$median_mean_delta,
      candidate_pool$rank
    ), , drop = FALSE]
  }
  candidate_ids <- head(candidate_pool$target_id, display_n)
} else {
  candidate_ids <- unique(plot_candidates$target_id[plot_candidates$unique_signatures == FALSE &
    plot_candidates$passes_targetability])
}
candidates <- screen[
  screen$unique_signatures == FALSE &
    screen$passes_targetability &
    screen$target_id %in% candidate_ids,
  ,
  drop = FALSE
]
candidates <- candidates[match(candidate_ids, candidates$target_id), , drop = FALSE]
candidates <- candidates[!is.na(candidates$target_id), , drop = FALSE]

response <- response[
  response$gene_set_size == 250 &
    response$scaling == "sample_pop_sd" &
    response$unique_signatures == FALSE &
    response$activity_target %in% candidates$target_id,
  ,
  drop = FALSE
]

summarise_response <- function(target_id, response_type) {
  z <- response[response$activity_target == target_id, , drop = FALSE]
  if (nrow(z) == 0) {
    return(data.frame(
      target_id = target_id,
      response_type = response_type,
      median_response = NA_real_,
      minimum_response = NA_real_,
      n_lines_positive = NA_integer_,
      n_lines_positive_both_doses = NA_integer_,
      n_lines_top5 = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  if (response_type == "Temsi 50/500 mean") {
    value <- z$mean_temsi_delta_vs_DMSO
    positive <- z$both_doses_up
    top5 <- if (screen_rule == "positive" || screen_rule == "either") {
      z$up_percentile_50nM >= 0.95 | z$up_percentile_500nM >= 0.95
    } else {
      z$both_doses_top5pct
    }
  } else if (response_type == "Temsi 50 nM") {
    value <- z$delta_50nM_vs_DMSO
    positive <- z$delta_50nM_vs_DMSO > 0
    top5 <- z$up_percentile_50nM >= 0.95
  } else if (response_type == "Temsi 500 nM") {
    value <- z$delta_500nM_vs_DMSO
    positive <- z$delta_500nM_vs_DMSO > 0
    top5 <- z$up_percentile_500nM >= 0.95
  } else {
    stop("Unknown response type: ", response_type)
  }

  keep <- is.finite(value)
  data.frame(
    target_id = target_id,
    response_type = response_type,
    median_response = median(value[keep]),
    minimum_response = min(value[keep]),
    n_lines_positive = sum(positive[keep], na.rm = TRUE),
    n_lines_positive_both_doses = sum(z$both_doses_up[keep], na.rm = TRUE),
    n_lines_top5 = sum(top5[keep], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

response_types <- c("Temsi 50/500 mean", "Temsi 50 nM", "Temsi 500 nM")
summary_table <- do.call(rbind, lapply(response_types, function(response_type) {
  do.call(rbind, lapply(candidates$target_id, summarise_response, response_type = response_type))
}))

summary_table <- merge(
  summary_table,
  candidates[, c(
    "target_id", "gene", "targetability_tier", "targetability_basis",
    "clinical_tractability_status", "screen_rank"
  )],
  by = "target_id",
  all.x = TRUE,
  sort = FALSE
)
summary_table <- summary_table[order(
  match(summary_table$response_type, response_types),
  -summary_table$median_response,
  -summary_table$n_lines_positive_both_doses,
  -summary_table$n_lines_top5
), , drop = FALSE]
summary_table$clinical_plot_group <- ifelse(
  summary_table$clinical_tractability_status %in% c(
    "Phase I clinical", "Advanced clinical", "Approved drug"
  ),
  "Clinical evidence (Phase I+)",
  "Direct/PROTAC/antibody tractability"
)

write.csv(
  summary_table,
  file.path(out_dir, paste0("vmpa_unique_FALSE_three_response_dotplot_data_n250_sample_pop_sd", file_tag, display_tag, size_tag, selection_tag, ".csv")),
  row.names = FALSE
)

make_plot <- function(response_type, file_stub) {
  z <- summary_table[summary_table$response_type == response_type, , drop = FALSE]
  size_column <- if (size_mode == "dose") "n_lines_positive" else "n_lines_positive_both_doses"
  z <- z[order(-z$median_response, -z[[size_column]], -z$n_lines_top5), , drop = FALSE]
  z$plot_label <- factor(z$target_id, levels = rev(z$target_id))
  x_label <- if (response_type == "Temsi 50 nM") {
    "Median VMPA activity-score delta\nTemsirolimus 50 nM vs DMSO"
  } else if (response_type == "Temsi 500 nM") {
    "Median VMPA activity-score delta\nTemsirolimus 500 nM vs DMSO"
  } else {
    "Median VMPA activity-score delta\nmean Temsirolimus 50/500 nM vs DMSO"
  }
  size_label <- if (size_mode == "dose") {
    "Cell lines positive at plotted dose (out of 4)"
  } else if (screen_rule == "positive") {
    "Cell lines positive at both doses (out of 4)"
  } else if (response_type == "Temsi 50/500 mean") {
    "Cell lines in top 5% at either dose"
  } else {
    paste0("Cell lines in top 5% at ", response_type)
  }

  p <- ggplot(z, aes(
    x = median_response,
    y = plot_label,
    color = clinical_plot_group,
    size = .data[[size_column]]
  )) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(shape = 16, alpha = 0.9) +
    scale_color_manual(
      values = c(
        "Direct/PROTAC/antibody tractability" = "#6baed6",
        "Clinical evidence (Phase I+)" = "#b2182b"
      ),
      name = "Targetability"
    ) +
    scale_size_continuous(range = c(3, 7), breaks = 0:4) +
    labs(
      x = x_label,
      y = NULL,
      size = size_label,
      title = paste0("Unbiased VMPA candidates: ", response_type),
      subtitle = if (screen_rule == "positive") {
        paste0(
          "unique=FALSE | 250 genes | sample_pop_sd\n",
          "Two-dose positive in >=3/4 lines | ",
          ifelse(size_mode == "dose", "size = dose-specific positives", "size = positives at both doses"),
          ifelse(selection_mode == "activity", " | activity-ranked", "")
        )
      } else {
        paste0(
          "unique=FALSE | 250 genes | sample_pop_sd\nTwo-dose positive in >=3/4 lines + top 5% at ",
          ifelse(screen_rule == "either", "either dose", "both doses")
        )
      }
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 8.5),
      legend.title = element_text(size = 9),
      axis.text.y = element_text(size = 10.5),
      axis.title.x = element_text(size = 10.5),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), p, width = 6.2, height = 5.0)
  ggsave(file.path(out_dir, paste0(file_stub, ".png")), p, width = 6.2, height = 5.0, dpi = 260)
}

if (screen_rule != "positive") {
  make_plot(
    "Temsi 50/500 mean",
    paste0("vmpa_unique_FALSE_mean_temsi_50_500_dotplot_n250_sample_pop_sd", file_tag, display_tag)
  )
}
make_plot(
  "Temsi 50 nM",
  if (screen_rule == "positive" && display_n == 10) {
    paste0("temsi_50", size_tag)
  } else if (screen_rule == "positive") {
    paste0("temsi_50_top", display_n, size_tag, selection_tag)
  } else {
    paste0("vmpa_unique_FALSE_temsi_50_dotplot_n250_sample_pop_sd", file_tag, display_tag)
  }
)
make_plot(
  "Temsi 500 nM",
  if (screen_rule == "positive" && display_n == 10) {
    paste0("temsi_500", size_tag)
  } else if (screen_rule == "positive") {
    paste0("temsi_500_top", display_n, size_tag, selection_tag)
  } else {
    paste0("vmpa_unique_FALSE_temsi_500_dotplot_n250_sample_pop_sd", file_tag, display_tag)
  }
)

cat("Wrote three unique=FALSE response dotplots and their data table.\n")
