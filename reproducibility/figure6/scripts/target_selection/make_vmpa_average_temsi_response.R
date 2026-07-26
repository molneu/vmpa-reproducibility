#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(here))

out_dir <- here::here("reproducibility", "figure6", "results", "target_selection", "vmpa_wb_correlations")
infile <- file.path(out_dir, "vmpa_treatment_response_all_settings.csv")
if (!file.exists(infile)) stop("Missing treatment-response table: ", infile)

response <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

key_cols <- c(
  "gene_set_size", "scaling", "unique_signatures",
  "cell_line", "activity_target", "target_gene"
)

wide_50 <- response[response$treatment == "50nM", c(key_cols, "delta_vs_DMSO", "rank_up", "up_percentile")]
wide_500 <- response[response$treatment == "500nM", c(key_cols, "delta_vs_DMSO", "rank_up", "up_percentile")]
names(wide_50)[names(wide_50) == "delta_vs_DMSO"] <- "delta_50nM_vs_DMSO"
names(wide_50)[names(wide_50) == "rank_up"] <- "rank_up_50nM"
names(wide_50)[names(wide_50) == "up_percentile"] <- "up_percentile_50nM"
names(wide_500)[names(wide_500) == "delta_vs_DMSO"] <- "delta_500nM_vs_DMSO"
names(wide_500)[names(wide_500) == "rank_up"] <- "rank_up_500nM"
names(wide_500)[names(wide_500) == "up_percentile"] <- "up_percentile_500nM"

avg <- merge(wide_50, wide_500, by = key_cols)
avg$mean_temsi_delta_vs_DMSO <- rowMeans(avg[, c("delta_50nM_vs_DMSO", "delta_500nM_vs_DMSO")])
avg$min_temsi_delta_vs_DMSO <- pmin(avg$delta_50nM_vs_DMSO, avg$delta_500nM_vs_DMSO)
avg$both_doses_up <- avg$delta_50nM_vs_DMSO > 0 & avg$delta_500nM_vs_DMSO > 0
avg$both_doses_top10pct <- avg$up_percentile_50nM >= 0.90 & avg$up_percentile_500nM >= 0.90
avg$both_doses_top5pct <- avg$up_percentile_50nM >= 0.95 & avg$up_percentile_500nM >= 0.95

rank_groups <- split(
  seq_len(nrow(avg)),
  paste(avg$cell_line, avg$gene_set_size, avg$scaling, avg$unique_signatures, sep = "\r")
)
avg$rank_mean_up <- NA_integer_
avg$rank_mean_down <- NA_integer_
avg$mean_up_percentile <- NA_real_
for (idx in rank_groups) {
  n <- length(idx)
  avg$rank_mean_up[idx] <- rank(-avg$mean_temsi_delta_vs_DMSO[idx], ties.method = "min")
  avg$rank_mean_down[idx] <- rank(avg$mean_temsi_delta_vs_DMSO[idx], ties.method = "min")
  avg$mean_up_percentile[idx] <- 1 - ((avg$rank_mean_up[idx] - 1) / pmax(n - 1, 1))
}

avg <- avg[order(avg$cell_line, avg$gene_set_size, avg$unique_signatures,
                 avg$scaling, avg$rank_mean_up), , drop = FALSE]
write.csv(avg, file.path(out_dir, "vmpa_temsi_50_500_average_response_by_setting.csv"), row.names = FALSE)

rank_cols <- c(
  "cell_line", "gene_set_size", "scaling", "unique_signatures",
  "rank_mean_up", "rank_mean_down", "mean_up_percentile",
  "activity_target", "target_gene",
  "mean_temsi_delta_vs_DMSO", "min_temsi_delta_vs_DMSO",
  "delta_50nM_vs_DMSO", "delta_500nM_vs_DMSO",
  "rank_up_50nM", "rank_up_500nM",
  "up_percentile_50nM", "up_percentile_500nM",
  "both_doses_up", "both_doses_top10pct", "both_doses_top5pct"
)
avg_rank_list <- avg[, intersect(rank_cols, names(avg)), drop = FALSE]
write.csv(avg_rank_list, file.path(out_dir, "vmpa_average_temsi_rank_list_all_conditions.csv"), row.names = FALSE)

avg_top_bottom <- rbind(
  transform(avg_rank_list[avg_rank_list$rank_mean_up <= 25, , drop = FALSE], direction = "top_up"),
  transform(avg_rank_list[avg_rank_list$rank_mean_down <= 25, , drop = FALSE], direction = "top_down")
)
avg_top_bottom <- avg_top_bottom[order(avg_top_bottom$cell_line, avg_top_bottom$gene_set_size,
                                       avg_top_bottom$unique_signatures, avg_top_bottom$scaling,
                                       avg_top_bottom$direction, avg_top_bottom$rank_mean_up,
                                       avg_top_bottom$rank_mean_down), , drop = FALSE]
write.csv(avg_top_bottom, file.path(out_dir, "vmpa_average_temsi_top25_up_down_all_conditions.csv"), row.names = FALSE)

summarise_group <- function(x) {
  data.frame(
    n_settings = nrow(x),
    n_both_doses_up = sum(x$both_doses_up),
    frac_both_doses_up = mean(x$both_doses_up),
    n_both_doses_top10pct = sum(x$both_doses_top10pct),
    frac_both_doses_top10pct = mean(x$both_doses_top10pct),
    n_both_doses_top5pct = sum(x$both_doses_top5pct),
    frac_both_doses_top5pct = mean(x$both_doses_top5pct),
    mean_temsi_delta_median = median(x$mean_temsi_delta_vs_DMSO),
    mean_temsi_delta_min = min(x$mean_temsi_delta_vs_DMSO),
    mean_temsi_delta_max = max(x$mean_temsi_delta_vs_DMSO),
    min_temsi_delta_median = median(x$min_temsi_delta_vs_DMSO),
    median_rank_50nM = median(x$rank_up_50nM),
    median_rank_500nM = median(x$rank_up_500nM),
    best_rank_50nM = min(x$rank_up_50nM),
    best_rank_500nM = min(x$rank_up_500nM),
    stringsAsFactors = FALSE
  )
}

split_key <- paste(avg$cell_line, avg$target_gene, sep = "\r")
summary_by_gene <- do.call(rbind, lapply(split(avg, split_key), function(x) {
  data.frame(
    cell_line = x$cell_line[1],
    target_gene = x$target_gene[1],
    summarise_group(x),
    check.names = FALSE
  )
}))
rownames(summary_by_gene) <- NULL
summary_by_gene <- summary_by_gene[order(
  summary_by_gene$cell_line,
  -summary_by_gene$frac_both_doses_up,
  -summary_by_gene$frac_both_doses_top10pct,
  -summary_by_gene$mean_temsi_delta_median,
  summary_by_gene$median_rank_50nM + summary_by_gene$median_rank_500nM
), , drop = FALSE]
write.csv(summary_by_gene, file.path(out_dir, "vmpa_temsi_consistently_up_by_cellline_gene.csv"), row.names = FALSE)

strict_up <- summary_by_gene[
  summary_by_gene$frac_both_doses_up == 1 &
    summary_by_gene$frac_both_doses_top10pct >= 0.5,
  ,
  drop = FALSE
]
write.csv(strict_up, file.path(out_dir, "vmpa_temsi_robustly_up_top10pct_by_cellline_gene.csv"), row.names = FALSE)

bcl2l1 <- summary_by_gene[toupper(summary_by_gene$target_gene) == "BCL2L1", , drop = FALSE]
write.csv(bcl2l1, file.path(out_dir, "BCL2L1_average_temsi_response_summary.csv"), row.names = FALSE)

for (cl in unique(summary_by_gene$cell_line)) {
  top <- summary_by_gene[summary_by_gene$cell_line == cl, , drop = FALSE]
  write.csv(head(top, 100), file.path(out_dir, paste0("top100_average_temsi_response_", cl, ".csv")), row.names = FALSE)
}

cat("Wrote average temsirolimus response summaries to", normalizePath(out_dir), "\n")
