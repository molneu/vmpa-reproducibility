#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biobase)
  library(GSVA)
  library(limma)
  library(cmapR)
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(patchwork)
  library(grid)
  library(here)
})

here::i_am("reproducibility/supplementary_figure8/run_supplementary_figure8.R")

plot_root_dir <- here("reproducibility", "supplementary_figure8", "results")
dir.create(plot_root_dir, recursive = TRUE, showWarnings = FALSE)
save_source_data <- FALSE
render_volcano_plots <- TRUE
render_gene_rank_plot <- TRUE

out_dir <- file.path(plot_root_dir, "large_font_30mm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(plot_root_dir, "plot_intermediates")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

plot_base_size <- 6
plot_title_size <- 5
plot_label_size <- 1.8
rank_base_size <- 5
plot_width <- 5.8
plot_height <- 1.85
rank_plot_width <- 150 / 25.4
rank_plot_height <- 60 / 25.4

eset <- readRDS(here("reproducibility", "figure6", "data", "ExpressionSet_GSE145128.rds"))
eset <- eset[, eset$relapse_TYPE != "eR"]
eset <- eset[, !(eset$patient_ID %in% c("B046", "B078"))]

expr <- exprs(eset)
meta <- pData(eset)

relapse_factor <- factor(meta$relapse_TYPE, levels = c("n", "cR"))
design <- model.matrix(~ relapse_factor)

format_p <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

clean_label <- function(x, max_chars = 30) {
  x <- gsub("^HALLMARK_|^REACTOME_|^BIOCARTA_|^PID_|^WP_|^KEGG_|^KEGG_MEDICUS_", "", x)
  x <- gsub("_", " ", x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

display_title <- function(x) {
  dplyr::case_when(
    x == "SigComLINCS_bottom250_gliomaDriverTargets" ~ "SigCom LINCS",
    TRUE ~ x
  )
}

volcano_colors <- c(
  Down = "#4575B4",
  NS = grDevices::adjustcolor("gray70", alpha.f = 0.15),
  Up = "#D73027"
)

make_vmpa_sets <- function(context, n = 250, min_conf = 1, driver_filter = TRUE) {
  subset_dir <- here("reproducibility", "figure3", "data", "subsets")
  gct <- readRDS(file.path(subset_dir, paste0(context, "_subset.rds")))

  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  if (driver_filter) {
    keep_idx <- keep_idx[gct@cdesc$cancer_driver_summary[keep_idx] != "None"]
  }

  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol

  ids <- gct@cdesc$id[keep_idx]
  conf <- gct@cdesc$cps_conf_total[keep_idx]
  targets <- gct@cdesc$cmap_name[keep_idx]
  names_out <- paste0(ids, "_c", conf)

  sets <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    ord <- suppressWarnings(order(vals, decreasing = FALSE, na.last = "keep"))
    unique(head(rownames(mat)[ord], n))
  })
  names(sets) <- names_out
  sets <- lapply(sets, function(x) intersect(x, rownames(expr)))
  keep <- lengths(sets) >= 10

  list(
    sets = sets[keep],
    metadata = tibble(signature = names_out[keep], target = targets[keep])
  )
}

make_sigcom_sets <- function(target_filter) {
  f <- here("reproducibility", "figure3", "data", "240918 LINC_gene sets_CRISPR_CONSENSUS_SigComLINCS_bottom_250.csv")
  raw <- read.csv(f, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  raw <- raw[, intersect(colnames(raw), target_filter), drop = FALSE]
  sets <- lapply(raw, function(x) unique(intersect(na.omit(x[x != ""]), rownames(expr))))
  sets <- sets[lengths(sets) >= 10]
  list(
    sets = sets,
    metadata = tibble(signature = names(sets), target = names(sets))
  )
}

run_signature_volcano <- function(collection_id, sets, metadata) {
  message("Processing ", collection_id, " with ", length(sets), " signatures")
  cache_file <- file.path(cache_dir, paste0(collection_id, "_signature_limma_tt.rds"))
  if (file.exists(cache_file)) {
    message("Using cached limma table: ", cache_file)
    tt <- readRDS(cache_file)
  } else {
    gsva_scores <- gsva(gsvaParam(expr, sets))
    n_features <- nrow(gsva_scores)
    scaled_scores <- scale(gsva_scores, scale = apply(gsva_scores, 2, sd) * sqrt((n_features - 1) / n_features))

    tt <- topTable(eBayes(lmFit(scaled_scores, design)), coef = 2, number = Inf, sort.by = "t") |>
      rownames_to_column("signature") |>
      left_join(metadata, by = "signature") |>
      mutate(
        target = ifelse(is.na(target), signature, target),
        logP = -log10(P.Value),
        color_group = case_when(
          logFC > 0 & P.Value < 0.05 ~ "Up",
          logFC < 0 & P.Value < 0.05 ~ "Down",
          TRUE ~ "NS"
        ),
        color_group = factor(color_group, levels = c("Down", "NS", "Up")),
        is_ctnnb1 = target == "CTNNB1" | grepl("^CTNNB1", signature)
      ) |>
      arrange(desc(t))
    saveRDS(tt, cache_file)
  }

  if (save_source_data) {
    write_csv(tt, file.path(out_dir, paste0(collection_id, "_GSVA_limma_results.csv")))
  }

  top_up <- tt |> filter(logFC > 0) |> arrange(desc(t)) |> slice_head(n = 5)
  top_down <- tt |> filter(logFC < 0) |> arrange(t) |> slice_head(n = 5)
  best_ctnnb1 <- tt |> filter(is_ctnnb1) |> arrange(desc(t)) |> slice_head(n = 1)

  label_rows <- bind_rows(
    top_up |> slice_head(n = 3),
    top_down |> slice_head(n = 3),
    best_ctnnb1
  ) |>
    distinct(signature, .keep_all = TRUE) |>
    mutate(label = clean_label(target))

  table_rows <- bind_rows(
    top_up |> transmute(Group = "Top up", Rank = as.character(row_number()), Target = clean_label(target), logFC = sprintf("%.3f", logFC), P = format_p(P.Value), FDR = format_p(adj.P.Val)),
    top_down |> transmute(Group = "Top down", Rank = as.character(row_number()), Target = clean_label(target), logFC = sprintf("%.3f", logFC), P = format_p(P.Value), FDR = format_p(adj.P.Val))
  ) |>
    distinct(Group, Target, .keep_all = TRUE)

  volcano <- ggplot(tt, aes(x = logFC, y = logP, color = color_group)) +
    geom_point(alpha = 0.45, size = 2) +
    geom_point(data = best_ctnnb1, aes(x = logFC, y = logP), inherit.aes = FALSE, shape = 21, fill = "gold", color = "black", stroke = 0.4, size = 2.8) +
    ggrepel::geom_text_repel(data = label_rows, aes(x = logFC, y = logP, label = label), inherit.aes = FALSE, size = plot_label_size, min.segment.length = 0, segment.color = "gray55", segment.size = 0.25, box.padding = 0.25, point.padding = 0.2, max.overlaps = Inf, seed = 1) +
    scale_color_manual(values = volcano_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    geom_hline(yintercept = 1.3, color = "gray50", linetype = "dashed") +
    labs(title = display_title(collection_id), x = "GSVA enrichment score difference", y = expression(-log[10]~~Raw~P-value), color = "") +
    theme_minimal(base_size = plot_base_size) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = plot_title_size, face = "plain"), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), legend.position = "none", aspect.ratio = 1)

  table_plot <- ggplot() +
    annotation_custom(
      gridExtra::tableGrob(
        table_rows,
        rows = NULL,
        theme = gridExtra::ttheme_minimal(
          base_size = 6,
          core = list(fg_params = list(hjust = 0, x = 0.02), padding = grid::unit(c(2, 2), "mm")),
          colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02), padding = grid::unit(c(2, 2), "mm"))
        )
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    labs(title = "Top signatures") +
    theme_void(base_size = 6) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 6), plot.margin = margin(2, 2, 2, 2))

  ggsave(file.path(out_dir, paste0(collection_id, "_GSVA_limma_volcano.pdf")), volcano + table_plot + patchwork::plot_layout(widths = c(1.0, 2.8)), width = plot_width, height = plot_height, dpi = 300)
}

run_context_plots <- function() {
  vmpa_glioma <- make_vmpa_sets("glioma")
  glioma_driver_targets <- sort(unique(vmpa_glioma$metadata$target))
  collections <- list(
    VMPA_glioma = vmpa_glioma,
    VMPA_melanoma = make_vmpa_sets("melanoma"),
    VMPA_nsclc = make_vmpa_sets("nsclc"),
    SigComLINCS_bottom250_gliomaDriverTargets = make_sigcom_sets(glioma_driver_targets)
  )
  invisible(lapply(names(collections), function(nm) {
    run_signature_volcano(nm, collections[[nm]]$sets, collections[[nm]]$metadata)
  }))
}

msig <- msigdbr(species = "Homo sapiens") |>
  distinct(gs_collection, gs_subcollection, gs_name, gene_symbol)

make_msig_sets <- function(gs_collection, gs_subcollection) {
  df <- msig |> filter(gs_collection == .env$gs_collection)
  if (!is.na(gs_subcollection) && grepl("[*]$", gs_subcollection)) {
    prefix <- sub("[*]$", "", gs_subcollection)
    df <- df |> filter(startsWith(gs_subcollection, prefix))
  } else if (!is.na(gs_subcollection)) {
    df <- df |> filter(gs_subcollection == .env$gs_subcollection)
  }
  sets <- split(df$gene_symbol, df$gs_name)
  sets <- lapply(sets, function(x) intersect(unique(x), rownames(expr)))
  sets[lengths(sets) >= 10 & lengths(sets) <= 500]
}

run_msig_volcano <- function(collection_id, gs_collection, gs_subcollection) {
  message("Processing ", collection_id)
  cache_file <- file.path(cache_dir, paste0(collection_id, "_msigdb_limma_tt.rds"))
  if (file.exists(cache_file)) {
    message("Using cached limma table: ", cache_file)
    tt <- readRDS(cache_file)
  } else {
    sets <- make_msig_sets(gs_collection, gs_subcollection)
    gsva_scores <- gsva(gsvaParam(expr, sets))
    n_features <- nrow(gsva_scores)
    scaled_scores <- scale(gsva_scores, scale = apply(gsva_scores, 2, sd) * sqrt((n_features - 1) / n_features))

    tt <- topTable(eBayes(lmFit(scaled_scores, design)), coef = 2, number = Inf, sort.by = "t") |>
      rownames_to_column("pathway") |>
      arrange(desc(t)) |>
      mutate(
        logP = -log10(P.Value),
        color_group = case_when(
          logFC > 0 & P.Value < 0.05 ~ "Up",
          logFC < 0 & P.Value < 0.05 ~ "Down",
          TRUE ~ "NS"
        ),
        color_group = factor(color_group, levels = c("Down", "NS", "Up")),
        is_wnt = grepl("WNT|BETA.?CATENIN|CTNNB|BCAT", pathway, ignore.case = TRUE),
        is_bcat = grepl("BCAT", pathway, ignore.case = TRUE)
      )
    saveRDS(tt, cache_file)
  }

  if (save_source_data) {
    write_csv(tt, file.path(out_dir, paste0(collection_id, "_GSVA_limma_results.csv")))
  }

  top_up <- tt |> filter(logFC > 0) |> arrange(desc(t)) |> slice_head(n = 5)
  top_down <- tt |> filter(logFC < 0) |> arrange(t) |> slice_head(n = 5)
  best_wnt <- tt |> filter(is_wnt) |> arrange(P.Value) |> slice_head(n = 1)
  highlight_rows <- if (collection_id == "C6") tt |> filter(is_bcat) else best_wnt

  label_rows <- bind_rows(
    top_up |> slice_head(n = 3),
    top_down |> slice_head(n = 3),
    highlight_rows
  ) |>
    distinct(pathway, .keep_all = TRUE) |>
    mutate(label = clean_label(pathway, max_chars = 14))

  table_rows <- bind_rows(
    top_up |> transmute(Group = "Top up", Rank = as.character(row_number()), Pathway = clean_label(pathway), logFC = sprintf("%.3f", logFC), P = format_p(P.Value), FDR = format_p(adj.P.Val)),
    top_down |> transmute(Group = "Top down", Rank = as.character(row_number()), Pathway = clean_label(pathway), logFC = sprintf("%.3f", logFC), P = format_p(P.Value), FDR = format_p(adj.P.Val))
  )

  volcano <- ggplot(tt, aes(x = logFC, y = logP, color = color_group)) +
    geom_point(alpha = 0.45, size = 2) +
    geom_point(data = highlight_rows, aes(x = logFC, y = logP), inherit.aes = FALSE, shape = 21, fill = "gold", color = "black", stroke = 0.4, size = 2.8) +
    ggrepel::geom_text_repel(data = label_rows, aes(x = logFC, y = logP, label = label), inherit.aes = FALSE, size = plot_label_size, min.segment.length = 0, segment.color = "gray55", segment.size = 0.25, box.padding = 0.25, point.padding = 0.2, max.overlaps = Inf, seed = 1) +
    scale_color_manual(values = volcano_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    geom_hline(yintercept = 1.3, color = "gray50", linetype = "dashed") +
    labs(title = display_title(collection_id), x = "GSVA enrichment score difference", y = expression(-log[10]~~Raw~P-value), color = "") +
    theme_minimal(base_size = plot_base_size) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = plot_title_size, face = "plain"), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), legend.position = "none", aspect.ratio = 1)

  table_plot <- ggplot() +
    annotation_custom(
      gridExtra::tableGrob(
        table_rows,
        rows = NULL,
        theme = gridExtra::ttheme_minimal(
          base_size = 6,
          core = list(fg_params = list(hjust = 0, x = 0.02), padding = grid::unit(c(2, 2), "mm")),
          colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02), padding = grid::unit(c(2, 2), "mm"))
        )
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    labs(title = "Top pathways") +
    theme_void(base_size = 6) +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 6), plot.margin = margin(2, 2, 2, 2))

  ggsave(file.path(out_dir, paste0(collection_id, "_GSVA_limma_volcano.pdf")), volcano + table_plot + patchwork::plot_layout(widths = c(1.0, 2.8)), width = plot_width, height = plot_height, dpi = 300)
}

run_msig_plots <- function() {
  specs <- tibble::tribble(
    ~collection_id, ~gs_collection, ~gs_subcollection,
    "HALLMARK", "H", NA_character_,
    "C2_CGP", "C2", "CGP",
    "C2_CP", "C2", "CP*",
    "C6", "C6", NA_character_
  )
  invisible(lapply(seq_len(nrow(specs)), function(i) {
    run_msig_volcano(specs$collection_id[i], specs$gs_collection[i], specs$gs_subcollection[i])
  }))
}

run_gene_rank_plot <- function() {
  cache_file <- file.path(cache_dir, "gene_rank_limma_tt.rds")
  if (file.exists(cache_file)) {
    message("Using cached gene-rank limma table: ", cache_file)
    gene_tt <- readRDS(cache_file)
  } else {
    meta_gene <- pData(eset)
    meta_gene$relapse_TYPE <- factor(meta_gene$relapse_TYPE, levels = c("n", "cR"))
    meta_gene$patient_ID <- factor(meta_gene$patient_ID)
    gene_design <- model.matrix(~ patient_ID + relapse_TYPE, data = meta_gene)

    gene_tt <- topTable(eBayes(lmFit(expr, gene_design)), coef = "relapse_TYPEcR", number = Inf, sort.by = "t") |>
      rownames_to_column("gene") |>
      mutate(
        rank_by_p = rank(P.Value, ties.method = "min"),
        direction = case_when(
          logFC > 0 ~ "Up in recurrent",
          logFC < 0 ~ "Down in recurrent",
          TRUE ~ "No change"
        )
      ) |>
      arrange(desc(t)) |>
      mutate(rank_up_to_down = row_number(), signed_score = t)
    saveRDS(gene_tt, cache_file)
  }

  wnt_genes <- c("CTNNB1", "CCND1", "WNT5B", "WNT5A", "LRP5", "LRP6", "FZD2", "FZD8", "SFRP2", "LGR5", "BCL9L", "ZNRF3", "DKK1", "LEF1", "MYC")
  highlight <- gene_tt |>
    filter(gene %in% wnt_genes) |>
    mutate(
      label = paste0(gene, " (", rank_up_to_down, ")"),
      label_x = case_when(
        gene == "SFRP2" ~ 350,
        gene == "LRP5" ~ 1550,
        gene == "LRP6" ~ 2750,
        gene == "ZNRF3" ~ 4100,
        gene == "FZD8" ~ 450,
        gene == "DKK1" ~ 1850,
        gene == "CTNNB1" ~ 3250,
        gene == "LEF1" ~ 4750,
        gene == "WNT5A" ~ 13250,
        gene == "MYC" ~ 14850,
        gene == "FZD2" ~ 16450,
        gene == "BCL9L" ~ 13900,
        gene == "WNT5B" ~ 15750,
        gene == "LGR5" ~ 17600,
        gene == "CCND1" ~ 16000,
        TRUE ~ rank_up_to_down
      ),
      label_y = case_when(
        gene == "SFRP2" ~ -0.85,
        gene == "LRP5" ~ -0.85,
        gene == "LRP6" ~ -0.85,
        gene == "ZNRF3" ~ -0.85,
        gene == "FZD8" ~ -2.0,
        gene == "DKK1" ~ -2.0,
        gene == "CTNNB1" ~ -2.0,
        gene == "LEF1" ~ -2.0,
        gene == "WNT5A" ~ 0.95,
        gene == "MYC" ~ 0.95,
        gene == "FZD2" ~ 0.95,
        gene == "BCL9L" ~ 2.1,
        gene == "WNT5B" ~ 2.1,
        gene == "LGR5" ~ 2.1,
        gene == "CCND1" ~ 3.25,
        TRUE ~ signed_score
      ),
      hjust = 0.5
    )

  if (save_source_data) {
    write_csv(highlight, file.path(out_dir, "selected_WNT_beta_catenin_gene_ranks.csv"))
  }

  p <- ggplot(gene_tt, aes(x = rank_up_to_down, y = signed_score)) +
    geom_col(aes(fill = signed_score > 0), width = 1, color = NA, alpha = 0.8) +
    annotate("rect", xmin = 0.5, xmax = 250.5, ymin = -Inf, ymax = Inf, fill = "#D73027", alpha = 0.06) +
    annotate("rect", xmin = nrow(gene_tt) - 250 + 0.5, xmax = nrow(gene_tt) + 0.5, ymin = -Inf, ymax = Inf, fill = "#4575B4", alpha = 0.06) +
    geom_hline(yintercept = 0, color = "gray35", linewidth = 0.18) +
    geom_vline(xintercept = c(250.5, nrow(gene_tt) - 250 + 0.5), color = "gray35", linetype = "dashed", linewidth = 0.18) +
    geom_point(data = highlight, aes(x = rank_up_to_down, y = 0), inherit.aes = FALSE, shape = 21, fill = "gold", color = "black", stroke = 0.15, size = 0.9) +
    geom_segment(data = highlight, aes(x = rank_up_to_down, xend = label_x, y = 0, yend = label_y), inherit.aes = FALSE, color = "gray35", linewidth = 0.15) +
    geom_text(
      data = highlight,
      aes(x = label_x, y = label_y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size = plot_label_size,
      color = "black"
    ) +
    scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"), labels = c("TRUE" = "Up in recurrent", "FALSE" = "Down in recurrent")) +
    labs(
      title = "Expression-ranked DEG list, all genes",
      subtitle = "Genes ranked from strongest up in recurrent to strongest down in recurrent; labels show rank in parentheses.",
      x = "Gene rank from top upregulated to top downregulated",
      y = "limma moderated t statistic",
      fill = ""
    ) +
    coord_cartesian(ylim = c(-13, 12), clip = "off") +
    theme_minimal(base_size = rank_base_size) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 5, face = "plain"), plot.subtitle = element_text(size = 5), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), legend.position = "bottom", plot.margin = margin(2, 5, 2, 5))

  ggsave(file.path(out_dir, "WNT_beta_catenin_gene_rank_barplot_all_genes.pdf"), p, width = rank_plot_width, height = rank_plot_height, dpi = 300)
  ggsave(file.path(out_dir, "WNT_beta_catenin_gene_rank_barplot_all_genes.svg"), p, width = rank_plot_width, height = rank_plot_height, dpi = 300)
}

if (isTRUE(render_volcano_plots)) {
  run_context_plots()
  run_msig_plots()
}
if (isTRUE(render_gene_rank_plot)) {
  run_gene_rank_plot()
}
