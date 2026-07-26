#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
})

out_dir <- here::here("reproducibility", "figure6", "results", "target_selection", "vmpa_wb_correlations")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

parse_activity_filename <- function(path) {
  base <- basename(path)
  m <- regexec("^vmpa_glioma_n([0-9]+)_(none|sample_z|sample_pop_sd|signature_z)_unique_(TRUE|FALSE)_activity_matrix\\.csv$", base)
  hit <- regmatches(base, m)[[1]]
  if (length(hit) == 0) stop("Unexpected activity filename: ", base)
  data.frame(
    gene_set_size = as.integer(hit[2]),
    scaling = hit[3],
    unique_signatures = as.logical(hit[4]),
    stringsAsFactors = FALSE
  )
}

sample_metadata <- function(sample_id) {
  data.frame(
    sample_id = sample_id,
    cell_line = sub("_(DMSO|50nM|500nM)$", "", sample_id),
    treatment = sub("^.*_(DMSO|50nM|500nM)$", "\\1", sample_id),
    stringsAsFactors = FALSE
  )
}

target_gene <- function(x) {
  sub("__.*$", "", x)
}

make_response_table <- function(path) {
  params <- parse_activity_filename(path)
  mat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  targets <- mat[[1]]
  scores <- as.matrix(mat[, -1, drop = FALSE])
  storage.mode(scores) <- "numeric"
  rownames(scores) <- targets

  sm <- sample_metadata(colnames(scores))
  response <- data.frame()
  for (cl in unique(sm$cell_line)) {
    dmso <- sm$sample_id[sm$cell_line == cl & sm$treatment == "DMSO"]
    if (length(dmso) != 1) next
    for (dose in c("50nM", "500nM")) {
      treated <- sm$sample_id[sm$cell_line == cl & sm$treatment == dose]
      if (length(treated) != 1) next
      delta <- scores[, treated] - scores[, dmso]
      tmp <- data.frame(
        params,
        cell_line = cl,
        treatment = dose,
        activity_target = rownames(scores),
        target_gene = target_gene(rownames(scores)),
        dmso_score = scores[, dmso],
        treated_score = scores[, treated],
        delta_vs_DMSO = as.numeric(delta),
        stringsAsFactors = FALSE
      )
      tmp$rank_up <- rank(-tmp$delta_vs_DMSO, ties.method = "min")
      tmp$rank_down <- rank(tmp$delta_vs_DMSO, ties.method = "min")
      tmp$n_targets_in_setting <- nrow(tmp)
      tmp$up_percentile <- 1 - ((tmp$rank_up - 1) / pmax(tmp$n_targets_in_setting - 1, 1))
      response <- rbind(response, tmp)
    }
  }
  response
}

activity_files <- list.files(
  out_dir,
  pattern = "^vmpa_glioma_n[0-9]+_.*_activity_matrix\\.csv$",
  full.names = TRUE
)
activity_files <- sort(activity_files)
if (length(activity_files) == 0) stop("No VMPA activity matrices found in ", out_dir)

message("Reading ", length(activity_files), " activity matrices...")
response <- do.call(rbind, lapply(activity_files, make_response_table))
response <- response[order(response$gene_set_size, response$unique_signatures, response$scaling,
                           response$cell_line, response$treatment, response$rank_up), , drop = FALSE]

write.csv(response, file.path(out_dir, "vmpa_treatment_response_all_settings.csv"), row.names = FALSE)

top_up <- response[response$rank_up <= 5, , drop = FALSE]
top_up$direction <- "top_up"
top_down <- response[response$rank_down <= 5, , drop = FALSE]
top_down$direction <- "top_down"
top_bottom <- rbind(top_up, top_down)
top_bottom <- top_bottom[order(top_bottom$gene_set_size, top_bottom$unique_signatures,
                               top_bottom$scaling, top_bottom$cell_line,
                               top_bottom$treatment, top_bottom$direction,
                               top_bottom$rank_up, top_bottom$rank_down), , drop = FALSE]
write.csv(top_bottom, file.path(out_dir, "vmpa_top5_up_down_by_cellline_condition_all_settings.csv"), row.names = FALSE)

bcl2l1 <- response[toupper(response$target_gene) == "BCL2L1", , drop = FALSE]
bcl2l1 <- bcl2l1[order(bcl2l1$cell_line, bcl2l1$treatment, bcl2l1$gene_set_size,
                       bcl2l1$unique_signatures, bcl2l1$scaling,
                       bcl2l1$rank_up), , drop = FALSE]
write.csv(bcl2l1, file.path(out_dir, "BCL2L1_response_ranks_all_settings.csv"), row.names = FALSE)

bcl2l1_summary <- aggregate(
  cbind(delta_vs_DMSO, rank_up, rank_down, up_percentile) ~
    cell_line + treatment + gene_set_size + scaling + unique_signatures,
  data = bcl2l1,
  FUN = function(x) c(min = min(x), median = median(x), max = max(x))
)
write.csv(bcl2l1_summary, file.path(out_dir, "BCL2L1_response_summary_by_setting.csv"), row.names = FALSE)

selected <- response[
  response$gene_set_size == 250 &
    response$scaling == "signature_z" &
    response$unique_signatures == TRUE,
  ,
  drop = FALSE
]
selected_top <- rbind(
  selected[selected$rank_up <= 5, , drop = FALSE],
  selected[selected$rank_down <= 5, , drop = FALSE]
)
selected_genes <- unique(selected_top$target_gene)
selected_heat <- selected[selected$target_gene %in% selected_genes, , drop = FALSE]
selected_heat$condition <- paste(selected_heat$cell_line, selected_heat$treatment, sep = "_")

if (nrow(selected_heat) > 0) {
  p <- ggplot(selected_heat, aes(x = condition, y = target_gene, fill = delta_vs_DMSO)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b83232", midpoint = 0) +
    labs(x = NULL, y = NULL, fill = "Delta vs DMSO",
         title = "Top/bottom VMPA activity responses",
         subtitle = "Glioma, n=250, signature_z, unique=TRUE") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  ggsave(file.path(plot_dir, "top_bottom_response_heatmap_n250_signature_z_unique_TRUE.pdf"),
         p, width = 10, height = max(5, 0.18 * length(unique(selected_heat$target_gene)) + 2))
}

if (nrow(bcl2l1) > 0) {
  bplot <- bcl2l1
  bplot$setting <- paste0("n", bplot$gene_set_size, "_", bplot$scaling, "_unique", bplot$unique_signatures)
  bplot$condition <- paste(bplot$cell_line, bplot$treatment, sep = "_")
  p2 <- ggplot(bplot, aes(x = setting, y = delta_vs_DMSO, color = cell_line, shape = treatment)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(size = 1.8, alpha = 0.85) +
    facet_wrap(~ condition, ncol = 4, scales = "free_y") +
    labs(x = NULL, y = "BCL2L1 activity delta vs DMSO",
         title = "BCL2L1 response across VMPA settings") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.minor = element_blank())
  ggsave(file.path(plot_dir, "BCL2L1_delta_across_all_settings.pdf"), p2, width = 14, height = 8)
}

message("Wrote treatment-response outputs to: ", normalizePath(out_dir))
