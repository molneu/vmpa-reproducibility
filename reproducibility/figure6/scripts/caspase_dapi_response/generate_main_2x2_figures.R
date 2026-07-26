#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("generate_main_2x2_figures.R")
script_dir <- dirname(script_path)

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
data_dir <- file.path(repo_root, "reproducibility", "figure6", "data", "caspase_dapi_response")
res_dir <- file.path(repo_root, "reproducibility", "figure6", "results", "caspase_dapi_response")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

data_file <- file.path(data_dir, "source_data_well_summed_BN118R_BN91R.csv")

source_data <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(
    cell_line = factor(cell_line, levels = c("BN118R", "BN91R")),
    Temsi_num = as.numeric(temsirolimus_nM),
    DT_num = as.numeric(DT2216_nM),
    Temsi = factor(Temsi_num, levels = c(0, 50, 500), labels = c("Temsi 0", "Temsi 50", "Temsi 500")),
    DT = factor(DT_num, levels = c(0, 500, 1000), labels = c("0", "500", "1000")),
    ratio_pct = as.numeric(red_blue_ratio_strict_pct),
    dapi_raw = as.numeric(blue_cells_sum)
  ) %>%
  group_by(cell_line, Temsi) %>%
  mutate(dapi_relative = dapi_raw / mean(dapi_raw[DT_num == 0], na.rm = TRUE)) %>%
  ungroup()

write_csv(source_data, file.path(res_dir, "main_plot_source_values.csv"))

compute_combined_temsi_interactions <- function(dat) {
  dat <- dat %>%
    mutate(
      temsi_present = factor(ifelse(Temsi_num == 0, "No temsi", "Temsi present"),
                             levels = c("No temsi", "Temsi present")),
      DTcat = factor(DT_num, levels = c(0, 500, 1000))
    )

  endpoints <- c(
    ratio_pct = "Caspase-3 red/blue ratio (%)",
    dapi_raw = "DAPI raw summed nuclei count",
    dapi_relative = "DAPI relative nuclei count"
  )

  rows <- list()
  k <- 1
  for (cl in levels(dat$cell_line)) {
    for (ep in names(endpoints)) {
      dd <- dat %>% filter(cell_line == cl)
      fit_cat <- lm(as.formula(paste(ep, "~ DTcat * temsi_present")), data = dd)
      a_cat <- anova(fit_cat)
      fit_num <- lm(as.formula(paste(ep, "~ DT_num * temsi_present")), data = dd)
      a_num <- anova(fit_num)

      rows[[k]] <- tibble(
        cell_line = cl,
        endpoint = endpoints[[ep]],
        endpoint_id = ep,
        categorical_interaction_F = a_cat["DTcat:temsi_present", "F value"],
        categorical_interaction_df1 = a_cat["DTcat:temsi_present", "Df"],
        categorical_interaction_df2 = a_cat["Residuals", "Df"],
        categorical_interaction_p = a_cat["DTcat:temsi_present", "Pr(>F)"],
        linear_slope_interaction_F = a_num["DT_num:temsi_present", "F value"],
        linear_slope_interaction_df1 = a_num["DT_num:temsi_present", "Df"],
        linear_slope_interaction_df2 = a_num["Residuals", "Df"],
        linear_slope_interaction_p = a_num["DT_num:temsi_present", "Pr(>F)"]
      )
      k <- k + 1
    }
  }

  bind_rows(rows) %>%
    mutate(
      categorical_BH_all = p.adjust(categorical_interaction_p, method = "BH"),
      linear_slope_BH_all = p.adjust(linear_slope_interaction_p, method = "BH")
    )
}

interactions <- compute_combined_temsi_interactions(source_data)
write_csv(interactions, file.path(res_dir, "combined_temsi_present_vs_absent_interaction_tests_recomputed.csv"))

label_data <- interactions %>%
  filter(endpoint_id %in% c("ratio_pct", "dapi_relative")) %>%
  mutate(
    assay = ifelse(endpoint_id == "ratio_pct", "ratio", "dapi"),
    label = paste0("interaction p = ", signif(categorical_interaction_p, 3),
                   "\nBH q = ", signif(categorical_BH_all, 3)),
    cell_line = factor(cell_line, levels = c("BN118R", "BN91R"))
  )
write_csv(label_data, file.path(res_dir, "interaction_pvalue_labels.csv"))

pal <- c("Temsi 0" = "#4B5563", "Temsi 50" = "#0072B2", "Temsi 500" = "#D55E00")

base_theme <- theme_classic(base_size = 5) +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5),
    axis.title = element_text(size = 5, colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    plot.title = element_text(face = "bold", size = 5),
    plot.subtitle = element_text(size = 5),
    plot.margin = margin(2, 2, 2, 2)
  )

error_df <- function(kind = c("SD", "SEM")) {
  kind <- match.arg(kind)
  function(x) {
    e <- sd(x)
    if (kind == "SEM") e <- e / sqrt(length(x))
    data.frame(y = mean(x), ymin = mean(x) - e, ymax = mean(x) + e)
  }
}

make_panel <- function(dat, labels, cell, y, assay_id, ylab, title, error_kind,
                       show_legend = FALSE, hline = FALSE, ylim = NULL,
                       show_title = TRUE) {
  dd <- dat %>% filter(cell_line == cell)
  y_top <- if (is.null(ylim)) max(dd[[y]], na.rm = TRUE) * 1.18 else ylim[2]
  lab <- labels %>% filter(cell_line == cell, assay == assay_id) %>% mutate(x = 1.08, y = y_top * 0.96)

  p <- ggplot(dd, aes(x = DT, y = .data[[y]], colour = Temsi, group = Temsi)) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.25) +
    stat_summary(fun = mean, geom = "point", size = 0.9) +
    stat_summary(fun.data = error_df(error_kind), geom = "errorbar", width = 0.11, linewidth = 0.15) +
    geom_text(data = lab, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 1.4, lineheight = 0.9) +
    scale_colour_manual(values = pal) +
    coord_cartesian(clip = "off", ylim = if (is.null(ylim)) c(0, y_top) else ylim) +
    labs(title = if (show_title) title else NULL, x = NULL, y = ylab) +
    base_theme

  if (hline) p <- p + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.15)
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

save_2x2 <- function(error_kind = c("SD", "SEM")) {
  error_kind <- match.arg(error_kind)
  suffix <- paste0("_", error_kind)

  p11 <- make_panel(source_data, label_data, "BN118R", "ratio_pct", "ratio",
                    "Caspase-3+ / nuclei (%)", "BN118R: Caspase-3 activation",
                    error_kind, show_legend = FALSE, ylim = c(0, 25), show_title = FALSE)
  p12 <- make_panel(source_data, label_data, "BN118R", "dapi_relative", "dapi",
                    "Relative DAPI nuclei count", "BN118R: Cell number",
                    error_kind, hline = TRUE, ylim = c(0, 1.5), show_title = FALSE)
  p21 <- make_panel(source_data, label_data, "BN91R", "ratio_pct", "ratio",
                    "Caspase-3+ / nuclei (%)", "BN91R: Caspase-3 activation",
                    error_kind, ylim = c(0, 25), show_title = FALSE)
  p22 <- make_panel(source_data, label_data, "BN91R", "dapi_relative", "dapi",
                    "Relative DAPI nuclei count", "BN91R: Cell number",
                    error_kind, hline = TRUE, ylim = c(0, 1.5), show_title = FALSE)

  pdf(file.path(res_dir, paste0("main_2x2_DT_x_axis_rows_celllines_cols_assays", suffix, ".pdf")),
      width = 6.4, height = 8.4)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(p11, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p12, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p21, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p22, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  dev.off()

  ggsave(file.path(res_dir, paste0("BN118R_ratio_DT_x_axis", suffix, ".pdf")), p11, width = 35, height = 25, units = "mm", device = "pdf")
  ggsave(file.path(res_dir, paste0("BN118R_dapi_relative_DT_x_axis", suffix, ".pdf")), p12, width = 35, height = 25, units = "mm", device = "pdf")
  ggsave(file.path(res_dir, paste0("BN91R_ratio_DT_x_axis", suffix, ".pdf")), p21, width = 35, height = 25, units = "mm", device = "pdf")
  ggsave(file.path(res_dir, paste0("BN91R_dapi_relative_DT_x_axis", suffix, ".pdf")), p22, width = 35, height = 25, units = "mm", device = "pdf")
}

save_2x2("SEM")

summary_values <- source_data %>%
  group_by(cell_line, Temsi, DT) %>%
  summarise(
    n = n(),
    ratio_mean = mean(ratio_pct),
    ratio_sd = sd(ratio_pct),
    ratio_sem = sd(ratio_pct) / sqrt(n()),
    dapi_relative_mean = mean(dapi_relative),
    dapi_relative_sd = sd(dapi_relative),
    dapi_relative_sem = sd(dapi_relative) / sqrt(n()),
    dapi_raw_mean = mean(dapi_raw),
    dapi_raw_sd = sd(dapi_raw),
    dapi_raw_sem = sd(dapi_raw) / sqrt(n()),
    .groups = "drop"
  )
write_csv(summary_values, file.path(res_dir, "main_plot_summary_mean_SD_SEM.csv"))

message("Wrote figures and recomputed statistics to: ", res_dir)
