#!/usr/bin/env Rscript
# Benchmarking COMPASS Gene Sets with fgsea (with null model)
# Author: Igor Cima
# Created: 2025-07-09    Updated: 2025-07-27

# ====== STEP 1: Load libraries ======
library(fgsea)            # Fast GSEA
library(DESeq2)           # Differential expression
library(BiocParallel)     # Parallel computing
library(data.table)       # Fast I/O
library(decoupleR)         # Network analysis
library(tidyverse)        # dplyr, purrr, tidyr, ggplot2, etc.
library(org.Hs.eg.db)     # Gene annotation
library(here)             # Project‐root relative paths

# ====== STEP 2: compass_gsc() ======

compass_gsc <- function(context,
                        subset_dir = here("reproducibility","figure3","data","subsets"),
                        n          = 200,
                        min_conf   = 1,
                        targets    = NULL,
                        output     = c("list","df","gsc")) {
  output <- match.arg(output)
  require(cmapR)
  require(GSEABase)
  
  # 2.1 Load the context‐specific subset
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) {
    stop("No subset file found for context: ", context)
  }
  gct <- readRDS(subset_file)
  
  # 2.2 Filter by confidence
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  if (length(keep_idx) == 0) {
    return(if (output=="gsc") GeneSetCollection(list()) else if (output=="df") data.frame() else list())
  }
  
  # 2.3 Extract matrix & set gene symbols as rownames
  mat <- gct@mat[, keep_idx, drop=FALSE]
  rownames(mat) <- gct@rdesc$symbol
  
  # 2.4 Build output names: "<id>_c<TOTAL_CONF>"
  ids       <- gct@cdesc$id[keep_idx]
  conf_tot  <- gct@cdesc$cps_conf_total[keep_idx]
  out_names <- paste0(ids, "_c", conf_tot)
  
  # 2.5 (Optional) filter to specified targets
  if (!is.null(targets)) {
    sel <- gct@cdesc$cmap_name[keep_idx] %in% targets
    mat <- mat[, sel, drop=FALSE]
    out_names <- out_names[sel]
    if (ncol(mat)==0) {
      warning("No signatures found for targets: ", paste(targets, collapse=","))
      return(if (output=="gsc") GeneSetCollection(list()) else if (output=="df") data.frame() else list())
    }
  }
  
  # 2.6 Extract bottom‐n genes per signature
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    ord <- order(mat[,j], decreasing=FALSE, na.last="keep")
    head(rownames(mat)[ord], n)
  })
  names(res_list) <- out_names
  
  # 2.7 Return in requested format
  if (output == "gsc") {
    gs <- lapply(out_names, function(nm) GeneSet(geneIds=res_list[[nm]], setName=nm))
    return(GeneSetCollection(gs))
  }
  if (output == "df") {
    max_len <- max(lengths(res_list))
    df <- as.data.frame(
      do.call(cbind, lapply(res_list, function(v) { length(v) <- max_len; v })),
      stringsAsFactors=FALSE,
      check.names=FALSE
    )
    return(df)
  }
  # default: list
  res_list
}

# ====== STEP 3: Setup directories & parameters ======
# parallel workers
nproc   <- BiocParallel::multicoreWorkers() - 1
BPPARAM <- MulticoreParam(workers = nproc)

# paths
data_dir         <- here("reproducibility","figure3","data")
subset_dir       <- file.path(data_dir, "subsets")
fig3_results_dir <- here("reproducibility","figure3","results")
dir.create(fig3_results_dir, recursive=TRUE, showWarnings=FALSE)

set.seed(123)

# ====== STEP 4: Load external gene sets ======
# 4.1 SigCom LINCS consensus LINC sets
sigcom_file   <- file.path(data_dir, "240918 LINC_gene sets_CRISPR_CONSENSUS_SigComLINCS_bottom_200.csv")
sigcom_all    <- read.csv(sigcom_file, row.names=1, stringsAsFactors=FALSE)

# 4.2 Collectri gene sets
collectri_file <- file.path(data_dir, "Collectri_genesets_forfGSEA.csv")
collectri_all  <- read.csv(collectri_file, row.names=1, stringsAsFactors=FALSE)

# 4.3 Collectri network
net <- get_collectri(organism="human", split_complexes=FALSE)

# ====== STEP 5: Load annotation data ======
annot_file <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv.gz")
if (!file.exists(annot_file)) {
  stop("Missing gene annotation file: ", annot_file)
}
annot <- fread(annot_file) %>% column_to_rownames("GeneID")

# ====== STEP 6: Load benchmarking parameters ======
params_file <- file.path(data_dir, "BENCHMARKING PARAMETERS ALL CONTEXTS_Collectri.csv")
params <- read.csv(params_file, stringsAsFactors=FALSE)

final_results <- vector("list", nrow(params))

# ====== STEP 7: Main benchmarking loop ======
for (i in seq_len(nrow(params))) {
  
  GENE         <- params$GENE[i]
  GSE          <- params$GEO[i]
  gsms         <- params$gsms[i]
  context      <- tolower(params$Context[i])
  perturbation <- params$Perturbation[i]
  
  message(">>> ", GENE, " @ ", GSE, " (", context, ")")
  
  # 7.1 Load & subset counts
  raw_url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts",
    "&acc=", GSE,
    "&file=", GSE, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"
  )
  tbl <- fread(raw_url) %>% column_to_rownames("GeneID") %>% as.matrix()
  sml <- strsplit(gsms, "")[[1]]
  sel <- sml != "X"
  sml <- sml[sel]
  tbl <- tbl[, sel, drop=FALSE]
  
  # 7.2 Filter low‐count genes
  groups    <- factor(ifelse(sml=="0","ctrl","treat"))
  min_samps <- min(table(groups))
  keep_genes <- rowSums(tbl >= 1) >= min_samps
  message("  dropping ", sum(!keep_genes), " low‐count genes")
  tbl <- tbl[keep_genes, , drop=FALSE]
  
  # 7.3 Build null counts by shuffling genes
  tbl_r <- tbl
  rownames(tbl_r) <- sample(rownames(tbl_r))
  
  # 7.4 DESeq2 on real & null
  coldata <- data.frame(Group=groups, row.names=colnames(tbl))
  ds      <- DESeqDataSetFromMatrix(tbl,   colData=coldata, design=~Group)
  ds_r    <- DESeqDataSetFromMatrix(tbl_r, colData=coldata, design=~Group)
  ds      <- DESeq(ds,   test="Wald", sfType="poscount", fitType="mean")
  ds_r    <- DESeq(ds_r, test="Wald", sfType="poscount", fitType="mean")
  
  # 7.5 Extract & annotate stats (unique per SYMBOL)
  extract_stats <- function(ds_obj) {
    res    <- as.data.frame(results(ds_obj))
    merged <- merge(res, annot, by="row.names", sort=FALSE)
    merged <- merged[is.finite(merged$stat), ]
    lst    <- split(merged$stat, merged$Symbol)
    sapply(lst, function(x) x[which.max(abs(x))])
  }
  stat_o   <- extract_stats(ds)
  stat_r_o <- extract_stats(ds_r)
  
  # 7.6 Prepare FGSEA collections filtered to GENE
  sizes        <- c(25,50,100,150,200,250,300)
  compass_sets <- map(sizes, function(n) {
    gs <- compass_gsc(context, subset_dir, min_conf=1, n=n, output="list")
    gs[grep(paste0("^",GENE,"_"), names(gs))]
  }) %>% set_names(paste0("b", sizes))
  
  sigcom    <- sigcom_all[, grepl(GENE, colnames(sigcom_all)), drop=FALSE]
  collectri  <- collectri_all[, grepl(GENE, colnames(collectri_all)), drop=FALSE]
  all_sets  <- c(compass_sets, list(SigCom=sigcom, Collectri=collectri))
  
  total_sets <- sum(map_int(all_sets, length))
  if (total_sets == 0) {
    warning("  no FGSEA sets for ", GENE); next
  }
  
  # 7.7 Run FGSEA (real & null)
  fg_real <- map(all_sets, ~ fgseaMultilevel(
    pathways   = .x,
    stats      = stat_o,
    minSize    = 10,
    maxSize    = 500,
    BPPARAM    = BPPARAM
  ))
  fg_null <- map(all_sets, ~ fgseaMultilevel(
    pathways   = .x,
    stats      = stat_r_o,
    minSize    = 10,
    maxSize    = 500,
    BPPARAM    = BPPARAM
  ))
  
  # 7.8 Extract FGSEA results
  extract_fgsea <- function(fg_list, flag) {
    map2_dfr(fg_list, names(fg_list), function(df, nm) {
      sub <- df %>%
        filter(str_detect(pathway,
                          paste0("^",GENE,"_|\\b",GENE,"\\b"))) %>%
        dplyr::select(pathway, pval, NES)
      if (nrow(sub) == 0) return(NULL)
      set_names(sub, paste0(c("pathway","pval","NES"), "_", nm)) %>%
        mutate(random_flag = flag)
    })
  }
  real_df <- extract_fgsea(fg_real, 1)
  null_df <- extract_fgsea(fg_null, 0)
  
  # 7.9 Run Collectri ULM (real & null)
  collectri_ulm   <- run_ulm(mat = stat_o,    net = net,
                             .source = 'source', .target = 'target',
                             .mor    = 'mor',
                             minsize = 5)
  collectri_ulm_r <- run_ulm(mat = stat_r_o,  net = net,
                             .source = 'source', .target = 'target',
                             .mor    = 'mor',
                             minsize = 5)
  
  # 7.10 Extract ULM results
  extract_ulm <- function(ulm_df, flag) {
    sub <- ulm_df %>%
      filter(str_detect(source, paste0("\\b",GENE,"\\b"))) %>%
      dplyr::select(source, p_value, score)
    if (nrow(sub) == 0) return(NULL)
    set_names(sub, paste0(c("pathway","pval","NES"), "_Collectri_ULM")) %>%
      mutate(random_flag = flag)
  }
  ulm_real <- extract_ulm(collectri_ulm,   1)
  ulm_null <- extract_ulm(collectri_ulm_r, 0)
  
  # 7.11 Combine all and attach metadata
  combined <- bind_rows(real_df, null_df, ulm_real, ulm_null)
  if (is.null(combined)) {
    warning("  no hits at all for ", GENE); next
  }
  combined <- combined %>%
    mutate(
      GSE          = GSE,
      GENE         = GENE,
      Context      = params$Context[i],
      Perturbation = perturbation,
      Samples      = gsms,
      n_ctrl       = sum(sml=="0"),
      n_treat      = sum(sml=="1")
    )
  final_results[[i]] <- combined
}

# ====== STEP 8: Wrap‐up & save results ======
final_df   <- bind_rows(final_results)

long_results <- final_df %>%
  pivot_longer(
    cols         = matches("^(pathway|pval|NES)_"),
    names_to     = c(".value", "collection"),
    names_pattern= "([^_]+)_(.+)"
  ) %>%
  filter(!is.na(pathway)) %>%
  mutate(
    expected_nes_sign = case_when(
      Perturbation %in% c("ACT","OE")      ~ "pos",
      Perturbation %in% c("INH","KD","KO") ~ "neg",
      TRUE                                  ~ NA_character_
    )
  ) %>%
  dplyr::select(
    GSE, GENE, Context, Perturbation, collection,
    pathway, pval, NES, random_flag, expected_nes_sign,
    Samples, n_ctrl, n_treat
  )

# write CSV to results directory
output_file <- file.path(data_dir, "250710_Benchmarking_fgsea_results_ALL_CONTEXTS.csv")
write_csv(long_results, output_file)
message("Done: wrote ", nrow(long_results), " rows to ", output_file)
