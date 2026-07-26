#!/usr/bin/env Rscript

# ===== 0) LIBRARIES =====
library(here)        # determine project root
library(readr)       # fast CSV I/O
library(dplyr)       # data manipulation (explicit dplyr:: prefixes)
library(precrec)     # PRC/ROC calculations (evalmod, auc)
library(ggplot2)     # plotting
library(patchwork)  # combine multiple ggplot2 plots

# ===== 1) SETUP =====

here::i_am("reproducibility/figure3/scripts/figure3e.R")

base_dir    <- here("reproducibility", "figure3")      # project figure3 folder
data_dir    <- file.path(base_dir, "data")              # data input path
results_dir <- file.path(base_dir, "results")           # output path
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)  # ensure exists

# ===== 2) READ DATA =====
nes_file  <- file.path(data_dir, "250721 Benchmarking_COMPASS_across_contexts_fgsea.csv")
null_file <- file.path(data_dir, "250725_FIG3E_COMPARISON_fGSEA_100NULLmodels_SUMMARY.csv")
nes_df       <- readr::read_csv(nes_file,  col_types = cols())  # real FGSEA results
null_summary <- readr::read_csv(null_file, col_types = cols())  # null-model summary

# ===== 3) MERGE & COMPUTE Z-SCORE =====
combined <- nes_df %>%
  dplyr::left_join(
    null_summary %>%
      dplyr::select(pathway, GeneSet_Context, GEO, grand_mean_NES, grand_sd_NES),
    by = c("pathway", "GSE" = "GEO")  # match gene-set ID and GSE/GEO
  ) %>%
  dplyr::mutate(
    z_score = (NES - grand_mean_NES) / grand_sd_NES  # standardize enrichment scores
  )

# ===== 4) SELECT HIGHEST-CONFIDENCE GENE-SETS & COMPUTE MEAN =====
# 4.1 extract confidence level from suffix "_cN"
combined_best_conf <- combined %>%
  dplyr::mutate(
    conf = as.integer(sub(".*_c(\\d+)$", "\\1", pathway))  # parse N from "_cN"
  )

# 4.2 keep only top confidence per GENE × Context (ties all kept)
best_conf_sets <- combined_best_conf %>%
  dplyr::group_by(GENE, Context, GeneSet_Context) %>%
  dplyr::filter(conf == max(conf, na.rm = TRUE)) %>%
  dplyr::ungroup()

# 4.3 reduce to high-confidence only
combined_highconf <- combined %>%
  dplyr::semi_join(
    best_conf_sets %>% dplyr::select(pathway),
    by = "pathway"
  )

# 4.4 compute mean z_score across replicates per GeneSet
combined_highconf_mean <- combined_highconf %>%
  dplyr::group_by(GSE, GENE, random_flag, Context, GeneSet_Context, Perturbation) %>%
  dplyr::summarise(
    z_score = mean(z_score, na.rm = TRUE),  # average across replicates
    .groups = "drop"
  )

# ===== 5) PRC IN/OUT-OF-CONTEXT & AUC COLLECTION =====
contexts  <- unique(combined_highconf_mean$GeneSet_Context)            # contexts to iterate
auc_table <- tibble::tibble(Context=character(), AUC_in=numeric(), AUC_out=numeric())
plots     <- vector("list", length(contexts)); names(plots) <- contexts

for (ctx in contexts) {
  # 5.1 subset to context and adjust z-scores
  df_ctx <- combined_highconf_mean %>%
    dplyr::filter(GeneSet_Context == ctx) %>%
    dplyr::mutate(
      z_score = -z_score,                              # flip sign for PR curves
      Label   = as.integer(random_flag == 1),           # 1=real, 0=null
      # invert again for activation perturbations
      z_score = if_else(Perturbation == "ACT", -z_score, z_score) #harmonize signs
    )
  
  # 5.2 split into in vs out of context
  in_ctx  <- df_ctx %>% dplyr::filter(Context == ctx)
  out_ctx <- df_ctx %>% dplyr::filter(Context != ctx)
  
  # 5.3 compute precision–recall curves
  pr_in  <- precrec::evalmod(scores=in_ctx$z_score,  labels=in_ctx$Label, modnames=ctx, mode="rocprc")
  pr_out <- precrec::evalmod(scores=out_ctx$z_score, labels=out_ctx$Label, modnames=ctx, mode="rocprc")
  
  # 5.4 extract AUC and sample counts
  auc_in  <- round(attr(pr_in,  "auc")[2,4], 3)
  auc_out <- round(attr(pr_out, "auc")[2,4], 3)
  n_in     <-       attr(pr_in,  "data_info")[4]
  n_out    <-       attr(pr_out, "data_info")[4]
  auc_table <- dplyr::bind_rows(
    auc_table,
    tibble::tibble(Context=ctx, AUC_in=auc_in, AUC_out=auc_out)
  )
  
  # 5.5 prepare data for plotting
  df1        <- fortify(pr_in)  %>% dplyr::mutate(Source="In context")
  df2        <- fortify(pr_out) %>% dplyr::mutate(Source="Out of context")
  plot_data  <- dplyr::bind_rows(df1, df2)
  
  # 5.6 build ggplot
  p <- ggplot2::ggplot(
    plot_data %>% dplyr::filter(curvetype == "PRC"),
    ggplot2::aes(x=x, y=y, color=Source)
  ) +
    ggplot2::geom_line(size=1) +
    ggplot2::labs(
      title = paste0("AUC in/out for ", ctx, ": ", auc_in, "/", auc_out),
      x     = "Recall",
      y     = "Precision",
      color = paste0(ctx, " n=", n_in, ",", n_out)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values=c(
      "In context"     = "#0072B2",
      "Out of context" = "#D55E00"
    )) +
    ggplot2::theme(
      legend.title      = ggplot2::element_text(size=10),
      legend.key.size   = grid::unit(0.5, "cm"),
      legend.text       = ggplot2::element_text(size=12),
      panel.border      = ggplot2::element_rect(colour="black", fill=NA, size=0.5),
      axis.text.x       = ggplot2::element_text(angle=45, hjust=1),
      axis.title        = ggplot2::element_text(size=14),
      axis.text         = ggplot2::element_text(size=14),
      axis.ticks        = ggplot2::element_line(color="black", size=0.5),
      axis.ticks.length = grid::unit(0.25, "cm")
    )
  
  plots[[ctx]] <- p  # store plot
}

# ===== 6) ASSEMBLE & SAVE ALL CONTEXT PLOTS =====
plots_clean   <- purrr::keep(plots, ~inherits(.x, "ggplot"))  # drop any NULLs
combined_plot <- patchwork::wrap_plots(plots_clean, ncol=4)     # 4 cols x 2 rows

ggplot2::ggsave(
  filename = file.path(results_dir, "Fig 3e_PRC_in_out_contexts.pdf"),
  plot     = combined_plot,
  width    = 20,   # 4 panels × 5 units width each
  height   = 8,    # 2 rows  × 4 units height each
  bg       = "white"
)

# ===== 7) PRINT AUC SUMMARY =====
auc_table <- auc_table %>%
  dplyr::mutate(Delta = AUC_in - AUC_out) %>%  # difference in performance
  dplyr::arrange(dplyr::desc(Delta))           # sort descending by Delta

message("AUCs for Precision-recall curves (in-context, out-of-context):")
print(auc_table)  # show in console
