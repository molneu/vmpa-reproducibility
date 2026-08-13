#!/usr/bin/env Rscript
# ============================================================
# Precision–Recall Curves for VIPER / COMPASS / SigCom / Collectri
# using Z-scores from benchmark file with null models
# - Uses 'collection_PR' column as method identifier
# - Flips Z-scores so higher = stronger activation
# - Computes AUCs and F1 maxima
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(precrec)
  library(ggplot2)
})

# ---------- PATHS ----------
results_dir <- here::here("reproducibility", "figure3", "results")
in_file <- file.path(results_dir, "251111_BENCHMARKING_withNullModels_Zscores.csv")
out_dir <- results_dir

# ---------- READ DATA ----------
df <- read_csv(in_file, show_col_types = FALSE)

# ---------- PROCESS DATA ----------
df2 <- df %>%
  dplyr::mutate(
    Z_adj = dplyr::if_else(expected_nes_sign == "pos", -z_score, z_score),
    Z_adj = -Z_adj,
    Label = dplyr::if_else(random_flag == 1, 1, 0),
    Method = collection_PR
  ) %>%
  dplyr::filter(!is.na(Z_adj), !is.na(Label), !is.na(Method))

# ---------- PRECREC INPUT ----------
score_list <- split(df2$Z_adj, df2$Method)
label_list <- split(df2$Label, df2$Method)
method_names <- names(score_list)

# ---------- RUN PR ANALYSIS ----------
precrec_obj <- precrec::evalmod(
  scores   = score_list,
  labels   = label_list,
  modnames = method_names,
  mode     = "rocprc"
)

# ---------- DEFINE CUSTOM COLOR PALETTE ----------
# Red shades for COMPASS, blue shades for others

method_colors <- c(
  "COMPASS_250"      = "#d73027",  # red-orange for COMPASS
  "VIPER_default"    = "#1f78b4",  # deep blue
  "VIPER_noFilter"   = "#9ad0ec",  # light sky blue
  "SigCom250"        = "#2e7d32",  # deep forest green (darker)
  "Collectri_ULM"    = "#9ccc65"   # light yellow-green (brighter)
)

# ---------- PLOT: Full PR curves ----------
msdf <- fortify(precrec_obj)

p_full <- ggplot(msdf %>% dplyr::filter(curvetype == "PRC"),
                 aes(x = x, y = y, color = modname)) +
  geom_line(size = 1.2) +
  labs(
    title = "Precision–Recall Curves by Method",
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = method_colors, name = "Method") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

p_full

ggsave(
  file.path(out_dir, "251111_PRcurve_byMethod_full_colored.pdf"),
  p_full, width = 6.5, height = 4.5, bg = "white"
)

# ---------- PLOT: Zoomed PR curves (0.5–1 precision range) ----------
p_zoom <- ggplot(msdf %>% dplyr::filter(curvetype == "PRC"),
                 aes(x = x, y = y, color = modname)) +
  geom_line(size = 1.2) +
  labs(
    title = "Precision–Recall Curves by Method (Zoomed)",
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = method_colors, name = "Method") +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  file.path(out_dir, "251111_PRcurve_byMethod_zoomed_colored.pdf"),
  p_zoom, width = 6.5, height = 4.5, bg = "white"
)

# ---------- AUC SUMMARY ----------
auc_df <- precrec::auc(precrec_obj)
cat("\nAUCs for Precision–Recall curves:\n")
print(auc_df %>% dplyr::filter(curvetypes == "PRC"))

auc_df <- precrec::auc(precrec_obj)
cat("\nAUCs for ROC curves:\n")
print(auc_df %>% dplyr::filter(curvetypes == "ROC"))

best <- auc_df %>%
  dplyr::filter(curvetypes == "PRC") %>%
  dplyr::slice_max(aucs, n = 1)
cat("\nBest performing method:",
    best$modnames, "with AUC =", round(best$aucs, 3), "\n")

n_by_method <- sapply(label_list, length)
n_by_method

# View(df2)

# ---------- F1 scores ----------
f1_df <- msdf %>%
  dplyr::filter(curvetype == "PRC") %>%
  dplyr::mutate(F1 = 2 * (x * y) / (x + y)) %>%
  dplyr::group_by(modname) %>%
  dplyr::summarise(
    F1_max = max(F1, na.rm = TRUE),
    Recall_at_maxF1 = x[which.max(F1)],
    Precision_at_maxF1 = y[which.max(F1)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(F1_max))

cat("\nMax F1 scores by method:\n")
print(f1_df)


# ============================================================
# ⚙️ Compute PR and F1 scores per context × method
# ============================================================

# Make sure we have labels and adjusted Z-scores ready
df2_precrec <- df2 %>%
  dplyr::filter(!is.na(Z_adj), !is.na(Label), !is.na(collection_PR), !is.na(Context))

# Split by context × method
df_split <- split(df2_precrec, interaction(df2_precrec$Context, df2_precrec$collection_PR, drop = TRUE))

# Evaluate PR curves for each combination
pr_list <- lapply(df_split, function(subdf) {
  if (length(unique(subdf$Label)) < 2) return(NULL)  # skip if all 0 or all 1
  precrec::evalmod(
    scores = subdf$Z_adj,
    labels = subdf$Label,
    modnames = unique(subdf$collection_PR),
    mode = "rocprc"
  )
})

# Combine results
msdf_context <- purrr::map_dfr(pr_list, fortify, .id = "Context_Method_ID")

# Extract Context and Method from the split ID
msdf_context <- msdf_context %>%
  tidyr::separate(Context_Method_ID, into = c("Context", "Method"), sep = "\\.")

# Compute F1 per context × method
f1_context_df <- msdf_context %>%
  dplyr::filter(curvetype == "PRC") %>%
  dplyr::mutate(F1 = 2 * (x * y) / (x + y)) %>%
  dplyr::group_by(Context, Method) %>%
  dplyr::summarise(
    F1_max = max(F1, na.rm = TRUE),
    Recall_at_maxF1 = x[which.max(F1)],
    Precision_at_maxF1 = y[which.max(F1)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(Context, desc(F1_max))

# Add total n (real + null) used for each context × method
n_counts <- df2_precrec %>%
  dplyr::group_by(Context, collection_PR) %>%
  dplyr::summarise(n_total = dplyr::n(), .groups = "drop") %>%
  dplyr::rename(Method = collection_PR)

f1_context_df <- f1_context_df %>%
  dplyr::left_join(n_counts, by = c("Context", "Method"))

cat("\n📊 F1 scores with total sample counts:\n")
# View(f1_context_df)

cat("\n📊 Context-specific F1 scores:\n")
# View(f1_context_df)

# ============================================================
# 🧮 Compute AUC (PRC and ROC) per context × method
# ============================================================

# Reuse the same pr_list used for F1 computation
auc_context_df <- purrr::map_dfr(names(pr_list), function(id) {
  obj <- pr_list[[id]]
  if (is.null(obj)) return(NULL)
  
  # Extract AUCs for both PRC and ROC
  auc_tbl <- precrec::auc(obj) %>%
    dplyr::filter(curvetypes %in% c("PRC", "ROC")) %>%
    dplyr::select(curvetypes, aucs, modnames)
  
  # Extract Context and Method from the split ID (e.g. glioma.COMPASS_250)
  parts <- strsplit(id, "\\.")[[1]]
  data.frame(
    Context = parts[1],
    Method  = parts[2],
    curvetypes = auc_tbl$curvetypes,
    AUC = auc_tbl$aucs,
    stringsAsFactors = FALSE
  )
})

# Pivot to wide format (optional, easier to merge)
auc_context_wide <- auc_context_df %>%
  tidyr::pivot_wider(names_from = curvetypes, values_from = AUC, names_prefix = "AUC_")

cat("\n📊 Context-specific AUCs:\n")
print(auc_context_wide)

# Compute total N per context × method
n_context_df <- df2_precrec %>%
  dplyr::group_by(Context, collection_PR) %>%
  dplyr::summarise(N_total = dplyr::n(), .groups = "drop") %>%
  dplyr::rename(Method = collection_PR)

# Merge everything
context_summary <- f1_context_df %>%
  dplyr::left_join(auc_context_wide, by = c("Context", "Method")) %>%
  dplyr::left_join(n_context_df, by = c("Context", "Method")) %>%
  dplyr::select(Context, Method, N_total, F1_max, AUC_PRC, AUC_ROC,
                Recall_at_maxF1, Precision_at_maxF1)

# Sort for readability
context_summary <- context_summary %>%
  dplyr::arrange(Context, desc(F1_max))

cat("\n✅ Summary per context × method:\n")
# View(context_summary)


# ============================================================
# 📊 Bar plot: F1_max per Method × Context (grouped by Context)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
})

# Order methods by average F1 across contexts
method_order <- f1_context_df %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(mean_F1 = mean(F1_max, na.rm = TRUE)) %>%
  dplyr::arrange(desc(mean_F1)) %>%
  dplyr::pull(Method)

# Prepare data for plotting
f1_plot_df <- f1_context_df %>%
  dplyr::mutate(
    Method  = factor(Method, levels = method_order),
    Context = forcats::fct_reorder(Context, F1_max, .fun = max, .desc = TRUE)
  )

# Method color palette (same as PR curves)
method_colors <- c(
  "COMPASS_250"      = "#d73027",
  "VIPER_default"    = "#1f78b4",
  "VIPER_noFilter"   = "#9ad0ec",
  "SigCom250"        = "#2e7d32",
  "Collectri_ULM"    = "#66bb6a"
)

# Create the plot
p_f1_bycontext <- ggplot(f1_plot_df, aes(x = Context, y = F1_max, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = method_colors, name = "Method") +
  labs(
    title = "F1 Scores by Method across Contexts",
    x = "Context",
    y = "F1 Score"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


p_f1_bycontext


# Compute expected random baseline F1 per context × method
random_baseline <- df2_precrec %>%
  dplyr::group_by(Context, collection_PR) %>%
  dplyr::summarise(
    prevalence = mean(Label == 1, na.rm = TRUE),
    expected_random_F1 = prevalence,  # ≈ prevalence baseline
    .groups = "drop"
  )

random_baseline %>%
  dplyr::arrange(desc(expected_random_F1)) %>%
  head()
