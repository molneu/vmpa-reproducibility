#!/usr/bin/env Rscript
# ============================================================
# Benchmarking COMPASS (b250), SigCom, Collectri, and 3-mode VIPER
# across multiple contexts using context-specific regulons
# - VIPER uses ARACNe regulons per context
# - COMPASS derives signatures from LINCS CMAP subsets
# - SigCom and Collectri gene sets are included for comparison
# - Includes both real (observed) and null (randomized) models
# - If a regulon is missing, VIPER is skipped but all other analyses run
# ============================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(DESeq2)
  library(viper)
  library(BiocParallel)
  library(data.table)
  library(decoupleR)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(here)
  library(aracne.networks)
  library(GSEABase)
  library(cmapR)
})

# ====== STEP 1: COMPASS helper function ======
# Extracts top-N up/down-regulated genes from CMAP subsets to form gene sets

compass_gsc <- function(context,
                        subset_dir,
                        n        = 200,
                        min_conf = 1,
                        targets  = NULL,
                        output   = c("list","df","gsc")) {
  output <- match.arg(output)
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) stop("No subset file found for context: ", context)
  gct <- readRDS(subset_file)
  
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  if (!length(keep_idx)) return(if (output=="gsc") GeneSetCollection(list()) else list())
  
  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol
  ids      <- gct@cdesc$id[keep_idx]
  conf_tot <- gct@cdesc$cps_conf_total[keep_idx]
  out_names <- paste0(ids, "_c", conf_tot)
  
  # Optional: filter by target genes
  if (!is.null(targets)) {
    sel <- gct@cdesc$cmap_name[keep_idx] %in% targets
    mat       <- mat[, sel, drop = FALSE]
    out_names <- out_names[sel]
    if (!ncol(mat)) {
      warning("No signatures found for targets: ", paste(targets, collapse=","))
      return(if (output=="gsc") GeneSetCollection(list()) else list())
    }
  }
  
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    ord <- order(mat[, j], decreasing = FALSE, na.last = "keep")
    head(rownames(mat)[ord], n)
  })
  names(res_list) <- out_names
  
  # Return GeneSetCollection, data.frame, or list
  if (output == "gsc") {
    gs <- lapply(out_names, function(nm)
      GeneSet(geneIds = res_list[[nm]], setName = nm))
    return(GeneSetCollection(gs))
  }
  if (output == "df") {
    max_len <- max(lengths(res_list))
    df <- as.data.frame(
      do.call(cbind, lapply(res_list, function(v) { length(v) <- max_len; v })),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(df)
  }
  res_list
}

# ====== STEP 2: Directory setup and configuration ======
data_dir   <- here("reproducibility","figure3","data")
subset_dir <- file.path(data_dir, "subsets")
results_dir <- here("reproducibility","figure3","results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)
nproc   <- BiocParallel::multicoreWorkers() - 1
BPPARAM <- MulticoreParam(workers = nproc)

# ====== STEP 3: Load ARACNe regulons per context ======
regulon_map <- list(
  glioma   = regulongbm,
  nsclc    = regulonluad,
  ovarian  = regulonov,
  pdac     = regulonpaad,
  crc      = reguloncoad,
  prostate = regulonprad,
  breast   = regulonbrca
)

# Precompute TF symbol map (ENTREZ → SYMBOL)
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

# ====== STEP 4: External gene sets (SigCom and Collectri) ======
sigcom_file    <- file.path(data_dir, "240918 LINC_gene sets_CRISPR_CONSENSUS_SigComLINCS_bottom_250.csv")
sigcom_all     <- read.csv(sigcom_file, row.names = 1, stringsAsFactors = FALSE)

collectri_file <- file.path(data_dir, "Collectri_genesets_forfGSEA.csv")
collectri_all  <- read.csv(collectri_file, row.names = 1, stringsAsFactors = FALSE)

net <- get_collectri(organism = "human", split_complexes = FALSE)

# ====== STEP 5: Gene annotation and parameters ======
annot_file <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv.gz")
if (!file.exists(annot_file)) {
  stop("Missing gene annotation file: ", annot_file)
}
annot <- fread(annot_file) %>% column_to_rownames("GeneID")

params_file <- file.path(data_dir, "BENCHMARKING PARAMETERS ALL CONTEXTS_Collectri.csv")
params <- read.csv(params_file, stringsAsFactors = FALSE)
final_results <- vector("list", nrow(params))

# ====== STEP 6: Helper for extracting DESeq2 stats ======
extract_stats <- function(ds_obj) {
  res    <- as.data.frame(results(ds_obj))
  merged <- merge(res, annot, by = "row.names", sort = FALSE)
  merged <- merged[is.finite(merged$stat), ]
  lst    <- split(merged$stat, merged$Symbol)
  sapply(lst, function(x) x[which.max(abs(x))])
}

# ====== STEP 7: Main loop over datasets ======
for (i in seq_len(nrow(params))) {
  GENE         <- params$GENE[i]
  GSE          <- params$GEO[i]
  gsms         <- params$gsms[i]
  context      <- tolower(params$Context[i])
  perturbation <- params$Perturbation[i]
  
  message("\n\n\n\n>>> ", GENE, " @ ", GSE, " (", context, ")")
  
  # ===== Load GEO data and DESeq2 =====
  raw_url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts",
    "&acc=", GSE,
    "&file=", GSE, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"
  )
  tbl <- fread(raw_url) %>% column_to_rownames("GeneID") %>% as.matrix()
  
  sml <- strsplit(gsms, "")[[1]]
  sel <- sml != "X"
  sml <- sml[sel]
  tbl <- tbl[, sel, drop = FALSE]
  
  groups    <- factor(ifelse(sml == "0", "ctrl", "treat"))
  min_samps <- min(table(groups))
  keep_genes <- rowSums(tbl >= 1) >= min_samps
  tbl <- tbl[keep_genes, , drop = FALSE]
  
  tbl_r <- tbl
  rownames(tbl_r) <- sample(rownames(tbl_r))
  
  coldata <- data.frame(Group = groups, row.names = colnames(tbl))
  ds      <- DESeqDataSetFromMatrix(tbl,   colData = coldata, design = ~ Group)
  ds_r    <- DESeqDataSetFromMatrix(tbl_r, colData = coldata, design = ~ Group)
  ds      <- DESeq(ds,   test = "Wald", sfType = "poscount", fitType = "mean")
  ds_r    <- DESeq(ds_r, test = "Wald", sfType = "poscount", fitType = "mean")
  
  stat_o   <- extract_stats(ds)
  stat_r_o <- extract_stats(ds_r)
  
  # ===== Map SYMBOL → ENTREZ =====
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
  signature_null <- stat_r_o[names(stat_r_o) %in% names(sig_entrez)]
  names(signature_null) <- sig_entrez[names(signature_null)]
  
  # ===== Context-specific regulon for VIPER =====
  regul <- regulon_map[[context]]
  if (!is.null(regul)) {
    message("  ✔ Using regulon for context: ", context)
    tf_map_ctx <- regulon_tf_symbol_map[[context]]
    tf_ids_match <- names(tf_map_ctx)[tf_map_ctx == GENE]
    if (length(tf_ids_match) == 0) {
      warning("No TF match for ", GENE, " in ", context, "; skipping VIPER.")
      regul_target <- NULL
    } else {
      regul_target <- regul[tf_ids_match]
      message("  ⚙️  Using TF-specific regulon for target gene: ", GENE)
    }
  } else {
    message("  ⚠️ No regulon for context: ", context, " — skipping VIPER.")
    regul_target <- NULL
  }
  
  # ===== VIPER (only if regul_target available) =====
  viper_real_wide <- list()
  viper_null_wide <- list()
  if (!is.null(regul_target) && length(signature_real) > 0) {
    message("  ▶ Running VIPER (real model)...")
    viper_modes <- list(
      VIPER_default  = function(sig) msviper(sig, regul_target, verbose = FALSE),
      VIPER_noFilter = function(sig) msviper(sig, regul_target, ges.filter = FALSE, verbose = FALSE),
      VIPER_adaptive = function(sig) msviper(sig, regul_target, adaptive.size = TRUE, verbose = FALSE)
    )
    viper_real_wide <- purrr::map2_dfr(
      viper_modes, names(viper_modes),
      ~ tryCatch(collect_viper_wide(.x(signature_real), .y, 1), error = function(e) NULL)
    )
  }
  if (!is.null(regul_target) && length(signature_null) > 0) {
    message("  ▶ Running VIPER (null model)...")
    viper_null_wide <- purrr::map2_dfr(
      viper_modes, names(viper_modes),
      ~ tryCatch(collect_viper_wide(.x(signature_null), .y, 0), error = function(e) NULL)
    )
  }
  
  # ===== FGSEA: COMPASS, SigCom, Collectri =====
  message("  ▶ Running FGSEA for context: ", context)
  sizes <- c(250)
  compass_sets <- purrr::map(sizes, function(n) {
    gs <- compass_gsc(context, subset_dir, min_conf = 1, n = n, output = "list")
    gs[grep(paste0("^", GENE, "_"), names(gs))]
  }) %>% set_names(paste0("b", sizes))
  
  sigcom    <- sigcom_all[, grepl(GENE, colnames(sigcom_all)), drop = FALSE]
  collectri <- collectri_all[, grepl(GENE, colnames(collectri_all)), drop = FALSE]
  all_sets  <- c(compass_sets, list(SigCom = sigcom, Collectri = collectri))
  
  fg_real <- purrr::map(all_sets, ~ fgseaMultilevel(
    pathways = .x, stats = stat_o, minSize = 10, maxSize = 500, BPPARAM = BPPARAM))
  fg_null <- purrr::map(all_sets, ~ fgseaMultilevel(
    pathways = .x, stats = stat_r_o, minSize = 10, maxSize = 500, BPPARAM = BPPARAM))
  
  extract_fgsea <- function(fg_list, flag) {
    purrr::map2_dfr(fg_list, names(fg_list), function(df, nm) {
      sub <- df %>% filter(str_detect(pathway, paste0("^", GENE, "_|\\b", GENE, "\\b"))) %>%
        dplyr::select(pathway, pval, NES)
      if (!nrow(sub)) return(NULL)
      sub %>% mutate(collection = nm, random_flag = flag)
    })
  }
  
  real_df <- extract_fgsea(fg_real, 1)
  null_df <- extract_fgsea(fg_null, 0)
  
  # ===== Collectri ULM =====
  collectri_ulm   <- run_ulm(
    mat = stat_o, net = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)
  collectri_ulm_r <- run_ulm(
    mat = stat_r_o, net = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)
  
  extract_ulm <- function(ulm_df, flag) {
    sub <- ulm_df %>% filter(str_detect(source, paste0("\\b", GENE, "\\b"))) %>%
      dplyr::select(source, p_value, score)
    if (!nrow(sub)) return(NULL)
    sub %>% rename(pathway = source, pval = p_value, NES = score) %>%
      mutate(collection = "Collectri_ULM", random_flag = flag)
  }
  
  ulm_real <- extract_ulm(collectri_ulm,   1)
  ulm_null <- extract_ulm(collectri_ulm_r, 0)
  
  # ===== Combine all results =====
  combined <- bind_rows(real_df, null_df, ulm_real, ulm_null, viper_real_wide, viper_null_wide)
  if (!nrow(combined)) next
  
  combined <- combined %>%
    mutate(GSE = GSE, GENE = GENE, Context = params$Context[i],
           Perturbation = perturbation, Samples = gsms,
           n_ctrl = sum(sml == "0"), n_treat = sum(sml == "1"))
  
  final_results[[i]] <- combined
}

# ====== STEP 8: Combine and save results ======
final_df <- bind_rows(final_results)
long_results <- final_df %>%
  mutate(expected_nes_sign = case_when(
    Perturbation %in% c("ACT","OE") ~ "pos",
    Perturbation %in% c("INH","KD","KO") ~ "neg",
    TRUE ~ NA_character_
  ))

output_file <- file.path(results_dir, "251109_Benchmarking_VIPER_COMPASSb250_ALL_CONTEXTS.csv")
write_csv(long_results, output_file)
message("✅ Done: wrote ", nrow(long_results), " rows to ", output_file)
