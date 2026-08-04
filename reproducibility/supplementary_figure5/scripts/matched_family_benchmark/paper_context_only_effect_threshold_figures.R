required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "forcats", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install missing packages first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
})

set.seed(1)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("paper_context_only_effect_threshold_figures.R")
}
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
data_dir <- file.path(repo_root, "reproducibility", "supplementary_figure5", "data", "matched_family_benchmark")
out_dir <- file.path(repo_root, "reproducibility", "supplementary_figure5", "results", "matched_family_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

excluded_contexts <- c("breast_ductal", "breast_invasive_unspecified")
excluded_prefix <- "^(esophageal|gastric|headneck)"

context_order <- c(
  "glioblastoma_idh_wt", "glioma_other",
  "nsclc_large_cell", "nsclc_adenocarcinoma", "nsclc_squamous", "nsclc_other",
  "melanoma_cutaneous", "melanoma_skin",
  "ovarian_hgsc", "ovarian_other",
  "pdac_adenocarcinoma",
  "crc_colon", "crc_other",
  "prostate_other",
  "breast_er_pr_positive_her2_negative", "breast_her2_positive_non_tnbc",
  "breast_tnbc", "breast_other"
)

context_labels <- c(
  glioblastoma_idh_wt = "Glioblastoma, IDH-wildtype",
  glioma_other = "Glioma, other",
  nsclc_large_cell = "Large-cell lung carcinoma",
  nsclc_adenocarcinoma = "Lung adenocarcinoma",
  nsclc_squamous = "Lung squamous cell carcinoma",
  nsclc_other = "NSCLC, other",
  melanoma_cutaneous = "Cutaneous melanoma",
  melanoma_skin = "Melanoma, skin",
  ovarian_hgsc = "Ovarian HGSC",
  ovarian_other = "Ovarian, other",
  pdac_adenocarcinoma = "PDAC adenocarcinoma",
  crc_colon = "CRC colon",
  crc_other = "CRC, other",
  prostate_other = "Prostate",
  breast_er_pr_positive_her2_negative = "Breast ER/PR+ HER2-",
  breast_her2_positive_non_tnbc = "Breast HER2+ non-TNBC",
  breast_tnbc = "Breast TNBC",
  breast_other = "Breast, other"
)

parent_labels <- c(
  glioma = "Glioma",
  nsclc = "NSCLC",
  melanoma = "Melanoma",
  ovarian = "Ovarian",
  pdac = "PDAC",
  crc = "CRC",
  prostate = "Prostate",
  breast = "Breast"
)

parent_order <- c("glioma", "nsclc", "melanoma", "ovarian", "pdac", "crc", "prostate", "breast")
parent_palette <- c(
  glioma = "#D95F5F",
  nsclc = "#B77A00",
  melanoma = "#2EA43B",
  ovarian = "#0AA99A",
  pdac = "#13A3C7",
  crc = "#4F86E8",
  prostate = "#8B6CE8",
  breast = "#D84FA3"
)

context_palette <- c(
  glioblastoma_idh_wt = "#D95F5F",
  glioma_other = "#C94B4B",
  nsclc_large_cell = "#C98300",
  nsclc_adenocarcinoma = "#A96F00",
  nsclc_squamous = "#8E6000",
  nsclc_other = "#D49A24",
  melanoma_cutaneous = "#2EA43B",
  melanoma_skin = "#66B85F",
  ovarian_hgsc = "#0AA99A",
  ovarian_other = "#43BFB4",
  pdac_adenocarcinoma = "#13A3C7",
  crc_colon = "#6FA2FF",
  crc_other = "#4F86E8",
  prostate_other = "#8B6CE8",
  breast_er_pr_positive_her2_negative = "#D84FA3",
  breast_her2_positive_non_tnbc = "#C94491",
  breast_tnbc = "#E06BB6",
  breast_other = "#B83D86"
)

context_to_parent <- function(x) {
  case_when(
    str_starts(x, "glioblastoma") | str_starts(x, "glioma") ~ "glioma",
    str_starts(x, "nsclc") ~ "nsclc",
    str_starts(x, "melanoma") ~ "melanoma",
    str_starts(x, "ovarian") ~ "ovarian",
    str_starts(x, "pdac") ~ "pdac",
    str_starts(x, "crc") ~ "crc",
    str_starts(x, "prostate") ~ "prostate",
    str_starts(x, "breast") ~ "breast",
    TRUE ~ "other"
  )
}

effect_levels <- c("Weak", "Moderate", "Strong")
effect_labels <- c(
  Weak = "Weak\nabs(r) >= 0.20",
  Moderate = "Moderate\nabs(r) >= 0.30",
  Strong = "Strong\nabs(r) >= 0.50"
)

source_levels <- c("COMPASS", "SigCom LINCS")
source_colors <- c("Random" = "#B8B8B8", "COMPASS" = "#527AA3", "SigCom LINCS" = "#B56B45")

median_delta_ci <- function(x, n_boot = 5000) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(low = NA_real_, high = NA_real_))
  boot <- replicate(n_boot, median(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  stats::quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

format_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), signif(p, 2))
  )
}

wide_all <- read_csv(
  file.path(data_dir, "compass_effect_threshold_sets_matched_values_wide.csv"),
  show_col_types = FALSE
)

wide <- wide_all %>%
  filter(
    !str_detect(compass_context, excluded_prefix),
    !compass_context %in% excluded_contexts
  ) %>%
  mutate(
    set_short = factor(set_short, levels = effect_levels),
    source_context = compass_context
  )

long <- wide %>%
  select(
    set_id, set_label, set_short, threshold, family_id, association_label,
    compass_context, comp_p, sigcom_p, comp_concordant, sigcom_concordant,
    paired_delta_abs_r, comp_abs_r, sigcom_abs_r
  ) %>%
  pivot_longer(
    cols = c(comp_abs_r, sigcom_abs_r),
    names_to = "source",
    values_to = "abs_r"
  ) %>%
  mutate(
    source = recode(source, comp_abs_r = "COMPASS", sigcom_abs_r = "SigCom LINCS"),
    source = factor(source, levels = source_levels),
    set_short = factor(set_short, levels = effect_levels)
  )

stats <- wide %>%
  group_by(set_id, set_label, set_short, threshold) %>%
  summarise(
    n_families = n(),
    median_compass_abs_r = median(comp_abs_r, na.rm = TRUE),
    iqr_compass_abs_r = IQR(comp_abs_r, na.rm = TRUE),
    median_sigcom_abs_r = median(sigcom_abs_r, na.rm = TRUE),
    iqr_sigcom_abs_r = IQR(sigcom_abs_r, na.rm = TRUE),
    median_paired_delta_abs_r = median(paired_delta_abs_r, na.rm = TRUE),
    wilcoxon_p = suppressWarnings(wilcox.test(comp_abs_r, sigcom_abs_r, paired = TRUE)$p.value),
    sigcom_same_direction_n = sum(sigcom_same_direction, na.rm = TRUE),
    sigcom_same_threshold_n = sum(sigcom_same_threshold, na.rm = TRUE),
    sigcom_same_threshold_fraction = sigcom_same_threshold_n / n_families,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ci = list(median_delta_ci(
      wide$paired_delta_abs_r[wide$set_id == set_id],
      n_boot = 5000
    )),
    median_delta_ci95_low = ci[[1]],
    median_delta_ci95_high = ci[[2]]
  ) %>%
  ungroup() %>%
  select(-ci) %>%
  mutate(
    wilcoxon_p_bh = p.adjust(wilcoxon_p, method = "BH"),
    annot = paste0(
      "n=", n_families,
      "\nmedian delta=", sprintf("%.2f", median_paired_delta_abs_r),
      " [", sprintf("%.2f", median_delta_ci95_low),
      ", ", sprintf("%.2f", median_delta_ci95_high), "]",
      "\nWilcoxon P=", format_p(wilcoxon_p),
      "\nSigCom same threshold=", sigcom_same_threshold_n, "/", n_families
    )
  )

write_csv(wide, file.path(out_dir, "paper_contexts_compass_effect_threshold_sets_matched_values_wide.csv"))
write_csv(long, file.path(out_dir, "paper_contexts_compass_effect_threshold_sets_matched_values_long.csv"))
write_csv(stats, file.path(out_dir, "paper_contexts_compass_effect_threshold_sets_paired_statistics.csv"))

p_box <- long %>%
  ggplot(aes(x = source, y = abs_r, fill = source)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.85, linewidth = 0.28) +
  geom_point(
    aes(group = family_id),
    position = position_jitter(width = 0.08, height = 0, seed = 1),
    size = 0.9,
    alpha = 0.42,
    stroke = 0
  ) +
  geom_line(
    aes(group = family_id),
    color = "grey70",
    linewidth = 0.18,
    alpha = 0.28
  ) +
  geom_text(
    data = stats,
    aes(x = 1.5, y = 0.98, label = annot),
    inherit.aes = FALSE,
    size = 2.75,
    lineheight = 0.92
  ) +
  facet_wrap(~ set_short, nrow = 1, labeller = labeller(set_short = effect_labels)) +
  scale_fill_manual(values = source_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.02, 0.06))) +
  labs(
    title = "COMPASS-positive families have higher matched RPPA correlations than SigCom LINCS",
    subtitle = "Paper contexts only; boxplots show matched absolute Pearson r across RPPA marker x activity x context families",
    x = NULL,
    y = "Absolute Pearson r"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    axis.text.x = element_text(color = "black", size = 8.5),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "grey94", color = "grey70", linewidth = 0.25),
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold"),
    panel.spacing.x = unit(0.8, "lines")
  )

ggsave(
  file.path(out_dir, "01_compass_effect_threshold_sets_compass_vs_sigcom_abs_r_boxplot_pairs_paper_contexts_only.pdf"),
  p_box,
  width = 9.2,
  height = 4.6,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "01_compass_effect_threshold_sets_compass_vs_sigcom_abs_r_boxplot_pairs_paper_contexts_only.png"),
  p_box,
  width = 9.2,
  height = 4.6,
  units = "in",
  dpi = 350
)

perm <- read_csv(
  file.path(data_dir, "compass_effect_threshold_sets_compass_sigcom_per_family_permutation_results.csv"),
  show_col_types = FALSE
) %>%
  filter(
    !str_detect(compass_context, excluded_prefix),
    !compass_context %in% excluded_contexts
  ) %>%
  mutate(
    source = factor(source, levels = source_levels),
    set_short = recode(set_id, weak = "Weak", moderate = "Moderate", strong = "Strong"),
    set_short = factor(set_short, levels = effect_levels)
  )

perm_summary <- perm %>%
  group_by(set_id, set_label, set_short, threshold, source) %>%
  summarise(
    n_families = n(),
    observed_count = sum(observed_success, na.rm = TRUE),
    expected_count = sum(null_probability_success, na.rm = TRUE),
    observed_expected_ratio = observed_count / expected_count,
    .groups = "drop"
  )

write_csv(perm, file.path(out_dir, "paper_contexts_compass_effect_threshold_sets_compass_sigcom_per_family_permutation_results.csv"))
write_csv(perm_summary, file.path(out_dir, "paper_contexts_compass_effect_threshold_sets_compass_sigcom_permutation_summary.csv"))

p_perm <- perm_summary %>%
  mutate(source = as.character(source)) %>%
  bind_rows(
    perm_summary %>%
      distinct(set_id, set_label, set_short, threshold) %>%
      mutate(
        source = "Random",
        n_families = NA_integer_,
        observed_count = NA_integer_,
        expected_count = NA_real_,
        observed_expected_ratio = 1
      )
  ) %>%
  mutate(source = factor(source, levels = c("Random", source_levels))) %>%
  ggplot(aes(x = set_short, y = observed_expected_ratio, fill = source)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed", color = "grey35") +
  geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "grey20", linewidth = 0.22) +
  geom_text(
    aes(label = paste0(sprintf("%.1f", observed_expected_ratio), "x")),
    position = position_dodge(width = 0.78),
    vjust = -0.3,
    size = 3.2
  ) +
  scale_fill_manual(values = source_colors) +
  scale_x_discrete(labels = effect_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Observed COMPASS-positive families exceed the permutation null",
    subtitle = "Paper contexts only; expected counts are summed from sample-label permutation probabilities",
    x = NULL,
    y = "Observed / expected by null",
    fill = "Signature source"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "02_compass_effect_threshold_sets_restricted_family_permutation_ratio_paper_contexts_only.pdf"),
  p_perm,
  width = 7.2,
  height = 4.7,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "02_compass_effect_threshold_sets_restricted_family_permutation_ratio_paper_contexts_only.png"),
  p_perm,
  width = 7.2,
  height = 4.7,
  units = "in",
  dpi = 350
)

select_top_direction <- function(df, direction = c("up", "down"), n = 5) {
  direction <- match.arg(direction)
  ranked <- if (direction == "up") {
    df %>%
      filter(comp_r >= 0.30) %>%
      arrange(desc(comp_r), comp_p)
  } else {
    df %>%
      filter(comp_r <= -0.30) %>%
      arrange(comp_r, comp_p)
  }
  if (nrow(ranked) == 0) return(ranked[0, , drop = FALSE])

  diverse <- ranked %>%
    distinct(compass_activity_gene, .keep_all = TRUE) %>%
    slice_head(n = n)
  if (nrow(diverse) >= n) return(diverse)

  fill <- ranked %>%
    filter(!family_id %in% diverse$family_id) %>%
    slice_head(n = n - nrow(diverse))
  bind_rows(diverse, fill)
}

dot_wide <- wide %>%
  filter(moderate) %>%
  arrange(desc(abs(comp_r)), comp_p) %>%
  distinct(family_id, .keep_all = TRUE) %>%
  group_by(compass_context) %>%
  group_modify(~ bind_rows(
    select_top_direction(.x, "up", n = 5),
    select_top_direction(.x, "down", n = 5)
  )) %>%
  ungroup() %>%
  mutate(
    compass_context = factor(compass_context, levels = context_order),
    parent_context = factor(context_to_parent(as.character(compass_context)), levels = parent_order),
    direction_group = if_else(comp_r >= 0, "Positive", "Negative"),
    direction_order = if_else(direction_group == "Positive", 1L, 2L),
    direction_rank = if_else(direction_group == "Positive", -comp_r, comp_r)
  ) %>%
  arrange(parent_context, compass_context, direction_order, direction_rank, comp_p) %>%
  mutate(
    marker_label = rppa_phosphoprotein %>%
      str_remove("_Caution$") %>%
      str_replace_all("_p", " p") %>%
      str_replace_all("_", " "),
    display_label = paste0(compass_activity_gene, " (", marker_label, "), r=", sprintf("%.2f", comp_r)),
    display_key = paste(display_label, compass_context, direction_group, family_id, sep = "___")
  ) %>%
  mutate(
    display_key = factor(display_key, levels = rev(unique(display_key))),
    y_pos = rev(row_number())
  )

dot_long <- dot_wide %>%
  select(
    family_id, association_label, display_label, display_key, y_pos, parent_context, compass_context, direction_group,
    comp_abs_r, sigcom_abs_r, comp_r, sigcom_r, comp_p, sigcom_p
  ) %>%
  pivot_longer(
    cols = c(comp_r, sigcom_r),
    names_to = "source",
    values_to = "pearson_r"
  ) %>%
  mutate(
    source = recode(source, comp_r = "COMPASS", sigcom_r = "SigCom LINCS"),
    source = factor(source, levels = rev(source_levels))
  )

write_csv(dot_wide, file.path(out_dir, "paper_contexts_compass_positive_unique_families_for_signed_dotplot.csv"))
write_csv(dot_long, file.path(out_dir, "paper_contexts_compass_positive_unique_families_for_signed_dotplot_long.csv"))

plot_height <- max(10, min(28, 2.6 + nrow(dot_wide) * 0.17))
x_limits <- range(dot_long$pearson_r, na.rm = TRUE)
x_limits <- c(
  max(-1, x_limits[1] - 0.03),
  min(1, x_limits[2] + 0.03)
)
x_breaks <- c(-0.75, -0.5, -0.3, 0.3, 0.5, 0.75)
x_breaks <- x_breaks[x_breaks >= x_limits[1] & x_breaks <= x_limits[2]]

parent_boundary_contexts <- c(
  "nsclc_large_cell",
  "melanoma_cutaneous",
  "ovarian_hgsc",
  "pdac_adenocarcinoma",
  "crc_colon",
  "prostate_other",
  "breast_er_pr_positive_her2_negative"
)
parent_boundary_lines <- tibble(
  compass_context = factor(parent_boundary_contexts, levels = context_order),
  yintercept = Inf
)

p_dot <- ggplot(
  dot_long,
  aes(x = pearson_r, y = display_key, color = compass_context, shape = source)
) +
  geom_hline(
    data = parent_boundary_lines,
    aes(yintercept = yintercept),
    color = "grey30",
    linewidth = 0.45,
    linetype = "dotted"
  ) +
  geom_segment(
    data = dot_wide,
    aes(
      x = comp_r,
      xend = sigcom_r,
      y = display_key,
      yend = display_key,
      group = family_id
    ),
    inherit.aes = FALSE,
    color = "grey72",
    linewidth = 0.38
  ) +
  geom_point(size = 5.2, stroke = 0.68, alpha = 0.95) +
  geom_vline(xintercept = c(-0.5, -0.3, 0.3, 0.5), color = "grey48", linewidth = 0.25, linetype = "dotted") +
  facet_grid(
    compass_context ~ .,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(compass_context = context_labels)
  ) +
  scale_color_manual(values = context_palette, drop = FALSE, guide = "none") +
  scale_shape_manual(values = c("SigCom LINCS" = 1, "COMPASS" = 16)) +
  scale_x_continuous(limits = x_limits, breaks = x_breaks, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  labs(
    title = "COMPASS-positive RPPA associations",
    subtitle = "Top positive/negative families per context; matched SigCom LINCS shown as open circles",
    x = "Pearson r",
    y = "Protein activity (marker)",
    shape = "Signature source"
  ) +
  theme_classic(base_size = 11.4) +
  theme(
    axis.text.x = element_text(color = "black", size = 9.8),
    axis.text.y = element_text(color = "black", size = 7.6),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.background = element_rect(fill = "grey94", color = "grey70", linewidth = 0.25),
    strip.text.y = element_text(face = "bold", size = 9.2, angle = 0),
    legend.position = "bottom",
    legend.title = element_text(size = 10.8),
    legend.text = element_text(size = 10.2),
    panel.spacing.y = unit(0.04, "lines"),
    panel.border = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    plot.title = element_text(face = "bold", size = 13.4),
    plot.subtitle = element_text(size = 10.6)
  )

ggsave(
  file.path(out_dir, "03_compass_positive_matched_sigcom_by_subcontext_signed_r_ge_0_30_paper_contexts_only.pdf"),
  p_dot,
  width = 9.6,
  height = plot_height,
  units = "in",
  device = "pdf",
  useDingbats = FALSE,
  limitsize = FALSE
)
ggsave(
  file.path(out_dir, "03_compass_positive_matched_sigcom_by_subcontext_signed_r_ge_0_30_paper_contexts_only.png"),
  p_dot,
  width = 9.6,
  height = plot_height,
  units = "in",
  dpi = 350,
  limitsize = FALSE
)

message("Paper-context plots written to: ", out_dir)
message("Dotplot families with |COMPASS r| >= 0.30: ", nrow(dot_wide))
print(stats)
print(perm_summary)
