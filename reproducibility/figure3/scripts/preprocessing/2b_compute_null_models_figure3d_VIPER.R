#!/usr/bin/env Rscript
# ============================================================
# VIPER Benchmark with 100 Null Models per Gene and Context
# - Runs VIPER in three modes on real data
# - Generates 100 permuted null models per mode, gene, and context
# - Uses ARACNe context-specific regulons (only TF-target subset)
# ============================================================

suppressPackageStartupMessages({
  library(viper)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(data.table)
  library(tidyverse)
  library(here)
  library(aracne.networks)
})

# ====== STEP 1: Directory setup ======
data_dir    <- here("reproducibility", "figure3", "data")
results_dir <- here("reproducibility", "figure3", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# ====== STEP 2: Load ARACNe regulons per context ======
regulon_map <- list(
  glioma   = regulongbm,
  nsclc    = regulonluad,
  ovarian  = regulonov,
  pdac     = regulonpaad,
  crc      = reguloncoad,
  prostate = regulonprad,
  breast   = regulonbrca
)

# ====== STEP 3: TF → SYMBOL map for regulons ======
regulon_tf_symbol_map <- lapply(regulon_map, function(regul_ctx) {
  tf_ids <- names(regul_ctx)
  tf_syms <- suppressMessages(mapIds(
    org.Hs.eg.db,
    keys     = tf_ids,
    column   = "SYMBOL",
    keytype  = "ENTREZID",
    multiVals = "first"
  ))
  tf_syms[!is.na(tf_syms)]
})

# ====== STEP 4: Parameters and annotation ======
params_file <- file.path(data_dir, "BENCHMARKING PARAMETERS ALL CONTEXTS_Collectri.csv")
params <- read.csv(params_file, stringsAsFactors = FALSE)

annot_file <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv.gz")
if (!file.exists(annot_file)) {
  stop("Missing gene annotation file: ", annot_file)
}
annot <- fread(annot_file) %>% column_to_rownames("GeneID")

# ====== STEP 5: Helper to extract DESeq2 stats ======
extract_stats <- function(ds_obj) {
  res    <- as.data.frame(results(ds_obj))
  merged <- merge(res, annot, by = "row.names", sort = FALSE)
  merged <- merged[is.finite(merged$stat), ]
  lst    <- split(merged$stat, merged$Symbol)
  sapply(lst, function(x) x[which.max(abs(x))])
}

# ====== STEP 6: Helper for VIPER output ======
collect_viper_wide <- function(msv, mode, flag, null_id = NA, gene) {
  if (inherits(msv, "error") || is.null(msv$es$nes)) return(NULL)
  df <- as.data.frame(msv$es$nes)
  df$TF_Entrez <- rownames(df)
  colnames(df)[1] <- "NES"
  df$pval <- msv$es$p.value
  df$Symbol <- suppressMessages(mapIds(
    org.Hs.eg.db,
    keys     = df$TF_Entrez,
    column   = "SYMBOL",
    keytype  = "ENTREZID",
    multiVals = "first"
  ))
  df <- df %>% filter(Symbol == gene)
  if (!nrow(df)) return(NULL)
  df %>%
    transmute(pathway = Symbol, pval = pval, NES = NES,
              collection = mode, random_flag = flag, null_id = null_id)
}

# ====== STEP 7: Main loop over datasets ======
final_results <- vector("list", nrow(params))

for (i in seq_len(nrow(params))) {
  GENE         <- params$GENE[i]
  GSE          <- params$GEO[i]
  gsms         <- params$gsms[i]
  context      <- tolower(params$Context[i])
  perturbation <- params$Perturbation[i]
  
  message("\n\n>>> ", GENE, " @ ", GSE, " (", context, ")")
  
  # ===== Load GEO data =====
  count_file <- file.path(data_dir, paste0(GSE, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
  if (!file.exists(count_file)) {
    stop("Missing raw count file: ", count_file)
  }
  tbl <- fread(count_file) %>% column_to_rownames("GeneID") %>% as.matrix()
  
  sml <- strsplit(gsms, "")[[1]]
  sel <- sml != "X"
  sml <- sml[sel]
  tbl <- tbl[, sel, drop = FALSE]
  
  groups    <- factor(ifelse(sml == "0", "ctrl", "treat"))
  min_samps <- min(table(groups))
  keep_genes <- rowSums(tbl >= 1) >= min_samps
  tbl <- tbl[keep_genes, , drop = FALSE]
  
  # ===== DESeq2 differential analysis =====
  coldata <- data.frame(Group = groups, row.names = colnames(tbl))
  ds      <- DESeqDataSetFromMatrix(tbl, colData = coldata, design = ~ Group)
  ds      <- DESeq(ds, test = "Wald", sfType = "poscount", fitType = "mean")
  stat_o  <- extract_stats(ds)
  
  # ===== SYMBOL → ENTREZ mapping =====
  sig_entrez <- suppressMessages(mapIds(
    org.Hs.eg.db,
    keys     = names(stat_o),
    column   = "ENTREZID",
    keytype  = "SYMBOL",
    multiVals = "first"
  ))
  sig_entrez <- sig_entrez[!is.na(sig_entrez)]
  signature_real <- stat_o[names(stat_o) %in% names(sig_entrez)]
  names(signature_real) <- sig_entrez[names(signature_real)]
  
  # ===== Regulon and TF subset =====
  regul <- regulon_map[[context]]
  if (is.null(regul)) {
    message("  ⚠️ No regulon for context: ", context, " — skipping VIPER.")
    next
  }
  
  tf_map_ctx <- regulon_tf_symbol_map[[context]]
  tf_ids_match <- names(tf_map_ctx)[tf_map_ctx == GENE]
  if (length(tf_ids_match) == 0) {
    warning("No TF match for ", GENE, " in ", context, "; skipping VIPER.")
    next
  }
  regul_target <- regul[tf_ids_match]
  message("  ✔ Using TF-specific regulon for target gene: ", GENE, " in context: ", context)
  
  # ===== Define VIPER modes =====
  viper_modes <- list(
    VIPER_default  = function(sig) msviper(sig, regul_target, verbose = FALSE),
    VIPER_noFilter = function(sig) msviper(sig, regul_target, ges.filter = FALSE, verbose = FALSE),
    VIPER_adaptive = function(sig) msviper(sig, regul_target, adaptive.size = TRUE, verbose = FALSE)
  )
  
  # ===== Real VIPER model =====
  message("  ▶ Running VIPER (real model)...")
  viper_real_wide <- purrr::map2_dfr(
    viper_modes, names(viper_modes),
    ~ tryCatch(collect_viper_wide(.x(signature_real), .y, 1, NA, GENE),
               error = function(e) NULL)
  )
  
  # ===== 100 Null VIPER models =====
  message("  ▶ Generating 100 VIPER null models...")
  null_models <- vector("list", 100)
  for (n in seq_len(100)) {
    permuted_signature <- setNames(sample(signature_real), names(signature_real))
    null_df <- purrr::map2_dfr(
      viper_modes, names(viper_modes),
      ~ tryCatch(
        collect_viper_wide(.x(permuted_signature), .y, 0, null_id = n, GENE),
        error = function(e) NULL
      )
    )
    null_models[[n]] <- null_df
    if (n %% 10 == 0) message("    ...completed ", n, "/100 nulls")
  }
  viper_null_wide <- dplyr::bind_rows(null_models)
  
  # ===== Combine & store =====
  combined <- bind_rows(viper_real_wide, viper_null_wide)
  combined <- combined %>%
    mutate(GSE = GSE, GENE = GENE, Context = params$Context[i],
           Perturbation = perturbation, Samples = gsms,
           n_ctrl = sum(sml == "0"), n_treat = sum(sml == "1"))
  
  final_results[[i]] <- combined
}

# ====== STEP 8: Save combined results ======
final_df <- bind_rows(final_results)
output_file <- file.path(results_dir, "251109_VIPER_100NullModels_ALL_CONTEXTS.csv")
write_csv(final_df, output_file)
message("✅ Done: wrote ", nrow(final_df), " rows to ", output_file)

# ====== STEP 9: Compute grand mean and SD of NES for VIPER null models ======
message("\n\n📊 Computing grand mean and SD of NES across 100 VIPER null models...")

# Filter only null models (random_flag == 0) and VIPER collections
viper_null_summary <- final_df %>%
  filter(random_flag == 0 & grepl("^VIPER_", collection)) %>%
  group_by(GSE, GENE, Context, collection) %>%
  summarise(
    grand_mean_NES = mean(NES, na.rm = TRUE),
    grand_sd_NES   = sd(NES, na.rm = TRUE),
    n_nulls        = n(),
    .groups = "drop"
  )

# Save summary table
summary_file <- file.path(results_dir, "251109_VIPER_100NullModels_NES_summary.csv")
write_csv(viper_null_summary, summary_file)

message("✅ Summary table saved: ", summary_file)
message("   → ", nrow(viper_null_summary), " rows summarised across all contexts.")



View(final_df %>%
  filter(random_flag == 0 & grepl("^VIPER_", collection)) %>%
  group_by(GSE, GENE, Context, collection) %>%
  summarise(
    n_nulls = n_distinct(null_id),
    nes_min = min(NES, na.rm = TRUE),
    nes_max = max(NES, na.rm = TRUE),
    unique_nes = n_distinct(round(NES, 4)),
    all_nas = all(is.na(NES))
  ))
