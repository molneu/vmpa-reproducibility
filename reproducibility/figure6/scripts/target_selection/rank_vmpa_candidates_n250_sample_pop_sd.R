#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

out_dir <- here::here("reproducibility", "figure6", "results", "target_selection", "vmpa_wb_correlations")
infile <- file.path(out_dir, "vmpa_temsi_50_500_average_response_by_setting.csv")
annotation_file <- file.path(out_dir, "vmpa_unbiased_top5pct_comparison_annotated.csv")
if (!file.exists(infile)) stop("Missing VMPA response file: ", infile)

suppressPackageStartupMessages(library(ggplot2))

x <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)
x <- x[x$gene_set_size == 250 & x$scaling == "sample_pop_sd", , drop = FALSE]
x$unique_signatures <- as.logical(x$unique_signatures)
x$gene <- sub("\\.[0-9]+$", "", toupper(x$target_gene))
x$top5_either_dose <- x$up_percentile_50nM >= 0.95 | x$up_percentile_500nM >= 0.95
x$top5_both_doses <- x$up_percentile_50nM >= 0.95 & x$up_percentile_500nM >= 0.95
x$up_both_doses <- x$delta_50nM_vs_DMSO > 0 & x$delta_500nM_vs_DMSO > 0

target_mode_key <- interaction(x$activity_target, x$unique_signatures, drop = TRUE)
ranked <- do.call(rbind, lapply(split(x, target_mode_key), function(z) {
  data.frame(
    target_id = z$activity_target[1],
    gene = z$gene[1],
    unique_signatures = z$unique_signatures[1],
    n_cell_lines = length(unique(z$cell_line)),
    n_lines_top5_both_doses = sum(z$top5_both_doses),
    n_lines_top5_either_dose = sum(z$top5_either_dose),
    n_lines_up_both_doses = sum(z$up_both_doses),
    median_mean_delta = median(z$mean_temsi_delta_vs_DMSO),
    minimum_mean_delta = min(z$mean_temsi_delta_vs_DMSO),
    median_mean_up_percentile = median(z$mean_up_percentile),
    median_rank_mean = median(z$rank_mean_up),
    stringsAsFactors = FALSE
  )
}))
rownames(ranked) <- NULL
ranked$signature_mode <- ifelse(ranked$unique_signatures, "unique=TRUE", "unique=FALSE")

# Transparent data-only score. Top-5% recurrence is weighted first, then
# positive direction, then whether the target enters the top 5% at either dose.
ranked$evidence_score <-
  4 * ranked$n_lines_top5_both_doses +
  2 * ranked$n_lines_up_both_doses +
  ranked$n_lines_top5_either_dose +
  ranked$median_mean_up_percentile
ranked <- ranked[order(
  -ranked$evidence_score,
  -ranked$n_lines_top5_both_doses,
  -ranked$n_lines_up_both_doses,
  -ranked$median_mean_delta,
  ranked$median_rank_mean
), , drop = FALSE]
ranked$rank <- ave(
  seq_len(nrow(ranked)),
  ranked$unique_signatures,
  FUN = seq_along
)

# Add cell-line-specific values for the comparison table.
for (cl in c("BN118", "BN118R", "BN91", "BN91R")) {
  z <- x[x$cell_line == cl, , drop = FALSE]
  key_z <- interaction(z$activity_target, z$unique_signatures, drop = TRUE)
  key_r <- interaction(ranked$target_id, ranked$unique_signatures, drop = TRUE)
  m <- match(key_r, key_z)
  ranked[[paste0(cl, "_top5_both")]] <- z$top5_both_doses[m]
  ranked[[paste0(cl, "_up_both")]] <- z$up_both_doses[m]
  ranked[[paste0(cl, "_mean_delta")]] <- z$mean_temsi_delta_vs_DMSO[m]
}

if (file.exists(annotation_file)) {
  ann <- read.csv(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
  ann <- ann[!duplicated(ann$gene), c(
    "gene", "oncogenic_annotation", "annotation_confidence",
    "mechanism_note", "actionability"
  ), drop = FALSE]
  ranked <- merge(ranked, ann, by = "gene", all.x = TRUE, sort = FALSE)
} else {
  ranked$oncogenic_annotation <- NA_character_
  ranked$annotation_confidence <- NA_character_
  ranked$mechanism_note <- NA_character_
  ranked$actionability <- NA_character_
}
ranked$oncogenic_annotation[is.na(ranked$oncogenic_annotation)] <- "not annotated here"
ranked$annotation_confidence[is.na(ranked$annotation_confidence)] <- "none"
ranked$mechanism_note[is.na(ranked$mechanism_note)] <- ""
ranked$actionability[is.na(ranked$actionability)] <- "not annotated here"
ranked <- ranked[order(ranked$unique_signatures, ranked$rank), , drop = FALSE]

write.csv(
  ranked,
  file.path(out_dir, "vmpa_ranked_candidates_n250_sample_pop_sd_unique_TRUE_FALSE.csv"),
  row.names = FALSE
)

top_plot <- do.call(rbind, lapply(c(FALSE, TRUE), function(mode) {
  head(ranked[ranked$unique_signatures == mode, , drop = FALSE], 20)
}))
top_plot$plot_label <- paste0(top_plot$target_id, " (", top_plot$signature_mode, ")")
write.csv(
  top_plot,
  file.path(out_dir, "vmpa_top20_candidates_n250_sample_pop_sd_unique_TRUE_FALSE.csv"),
  row.names = FALSE
)

dot <- ggplot(
  top_plot,
  aes(
    x = median_mean_delta,
    y = reorder(plot_label, evidence_score),
    color = n_lines_top5_both_doses,
    size = n_lines_up_both_doses
  )
) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_point(alpha = 0.9) +
  facet_wrap(~ signature_mode, ncol = 2, scales = "free_y") +
  scale_color_gradient(low = "#9ecae1", high = "#b2182b", breaks = 0:4) +
  scale_size_continuous(range = c(2.5, 6), breaks = 0:4) +
  labs(
    x = "Median VMPA activity delta: mean Temsi 50/500 nM vs DMSO",
    y = NULL,
    color = "Cell lines in top 5% at both doses",
    size = "Cell lines positive at both doses",
    title = "Ranked VMPA candidates",
    subtitle = "Glioma context; gene-set size 250; sample population SD; top 20 per signature mode"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

bar <- ggplot(
  top_plot,
  aes(
    x = reorder(plot_label, evidence_score),
    y = evidence_score,
    fill = n_lines_top5_both_doses
  )
) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.1f", evidence_score)), hjust = -0.15, size = 2.5) +
  facet_wrap(~ signature_mode, ncol = 2, scales = "free_x") +
  scale_fill_gradient(low = "#9ecae1", high = "#b2182b", breaks = 0:4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Evidence score",
    fill = "Cell lines in top 5% at both doses",
    title = "VMPA candidate evidence score",
    subtitle = "Score = 4 x top-5% recurrence + 2 x positive direction + either-dose recurrence + percentile"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "vmpa_candidate_dotplot_n250_sample_pop_sd.pdf"), dot, width = 12, height = 10)
ggsave(file.path(out_dir, "vmpa_candidate_dotplot_n250_sample_pop_sd.png"), dot, width = 12, height = 10, dpi = 220)
ggsave(file.path(out_dir, "vmpa_candidate_barplot_n250_sample_pop_sd.pdf"), bar, width = 12, height = 10)
ggsave(file.path(out_dir, "vmpa_candidate_barplot_n250_sample_pop_sd.png"), bar, width = 12, height = 10, dpi = 220)

cat("Wrote ranked n=250 sample_pop_sd candidate tables and plots to ", normalizePath(out_dir), "\n", sep = "")
