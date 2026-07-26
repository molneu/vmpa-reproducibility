#!/usr/bin/env Rscript
# Benchmarking COMPASS gene sets across contexts with fgsea
# Author: Igor Cima (adapted)
# Date: 2025-07-21

# ====== LIBRARIES ======
library(fgsea)
library(DESeq2)
library(BiocParallel)
library(data.table)
library(tidyverse)
library(org.Hs.eg.db)
library(tidyr)

# ====== compass_gsc function ======
# (unchanged from your version)
compass_gsc <- function(context,
                        subset_dir = "subsets",
                        n          = 200,
                        min_conf   = 1,
                        targets    = NULL,
                        output     = c("list","df","gsc")) {
  output <- match.arg(output)
  require(cmapR); require(GSEABase)
  subset_file <- file.path(subset_dir, paste0(context,"_subset.rds"))
  if (!file.exists(subset_file)) stop("No subset for ", context)
  gct <- readRDS(subset_file)
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  mat      <- gct@mat[, keep_idx, drop=FALSE]
  rownames(mat) <- gct@rdesc$symbol
  ids      <- gct@cdesc$id[keep_idx]
  conf_tot <- gct@cdesc$cps_conf_total[keep_idx]
  out_names<- paste0(ids, "_c", conf_tot)
  if (!is.null(targets)) {
    keep_tgt <- gct@cdesc$cmap_name[keep_idx] %in% targets
    mat       <- mat[, keep_tgt, drop=FALSE]
    out_names <- out_names[keep_tgt]
  }
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[,j])
    ord  <- suppressWarnings(order(vals, decreasing=FALSE, na.last="keep"))
    head(rownames(mat)[ord], n)
  })
  names(res_list) <- out_names
  if (output=="gsc") return(GeneSetCollection(
    lapply(out_names, function(nm) GeneSet(res_list[[nm]], setName=nm))
  ))
  if (output=="df") {
    max_len <- max(lengths(res_list))
    df <- as.data.frame(do.call(cbind, lapply(res_list, function(x){
      length(x) <- max_len; x
    })), stringsAsFactors=FALSE, check.names=FALSE)
    return(df)
  }
  res_list
}

# ====== SETUP ======
# parallel FGSEA
nproc   <- BiocParallel::multicoreWorkers() - 1
BPPARAM <- MulticoreParam(workers = nproc)
DIR     <- here::here("reproducibility", "figure3", "data", "subsets")
set.seed(123)

# ====== LOAD EXTERNAL ANNOTATIONS & PARAMETERS ======
annot_path <- here::here("reproducibility", "figure3", "data", "Human.GRCh38.p13.annot.tsv.gz")
if (!file.exists(annot_path)) {
  stop("Missing gene annotation file: ", annot_path)
}
annot <- fread(annot_path) %>% column_to_rownames("GeneID")

params <- read.csv(here::here("reproducibility", "figure3", "data", "BENCHMARKING PARAMETERS ALL CONTEXTS_Collectri.csv"),
  stringsAsFactors = FALSE
)

# collect all possible contexts
all_contexts <- unique(tolower(params$Context))

# ====== MAIN LOOP ======
final_results <- vector("list", nrow(params))

for (i in seq_len(nrow(params))) {
  GENE    <- params$GENE[i]
  GSE     <- params$GEO[i]
  gsms    <- params$gsms[i]
  context <- tolower(params$Context[i])
  perturbation <- params$Perturbation[i]
  
  message(">>> ", GENE, "@", GSE, " (", context, ")")
  # 1) Load & subset counts
  url_base <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
  path_raw <- paste0(url_base, "&acc=", GSE,
                     "&file=", GSE, "_raw_counts_GRCh38.p13_NCBI.tsv.gz")
  tbl <- fread(path_raw) %>%
    column_to_rownames("GeneID") %>%
    as.matrix()
  
  sml <- strsplit(gsms, "")[[1]]
  sel <- sml != "X"
  sml <- sml[sel]
  tbl <- tbl[, sel, drop=FALSE]
  
  # 2) Filter low‐count genes
  groups    <- factor(ifelse(sml=="0","ctrl","treat"))
  min_samps <- min(table(groups))
  keep_g    <- rowSums(tbl >= 1) >= min_samps
  tbl       <- tbl[keep_g, , drop=FALSE]
  
  # 3) Scrambled null (shuffle rownames)
  tbl_r <- tbl
  rownames(tbl_r) <- sample(rownames(tbl_r))
  
  # 4) DESeq2 on real & null
  coldata <- data.frame(Group=groups, row.names=colnames(tbl))
  ds       <- DESeqDataSetFromMatrix(tbl,   colData=coldata, design=~Group)
  ds_r     <- DESeqDataSetFromMatrix(tbl_r, colData=coldata, design=~Group)
  ds   <- DESeq(ds,   test="Wald", sfType="poscount", fitType="mean")
  ds_r <- DESeq(ds_r, test="Wald", sfType="poscount", fitType="mean")
  
  # 5) Extract & annotate stats (unique SYMBOLS)
  extract_stats <- function(ds_obj) {
    res    <- as.data.frame(results(ds_obj))
    merged <- merge(res, annot, by.x="row.names", by.y="row.names", sort=FALSE)
    merged <- merged[is.finite(merged$stat), ]
    lst    <- split(merged$stat, merged$Symbol)
    sapply(lst, function(x) x[which.max(abs(x))])
  }
  stat_o   <- extract_stats(ds)
  stat_r_o <- extract_stats(ds_r)
  
  # 6) Build 250-gene COMPASS buckets for this GENE across every context
  compass_sets <- list()
  for (ctx in all_contexts) {
    tmp <- compass_gsc(ctx, DIR, n = 250, output = "list")
    hits <- grep(paste0("^", GENE, "_"), names(tmp), value = TRUE)
    if (length(hits)) {
      compass_sets <- c(compass_sets, tmp[hits])
    }
  }
  
  # 7) Run FGSEA (real & null)
  fg_real <- fgseaMultilevel(
    pathways = compass_sets,
    stats    = stat_o,
    minSize  = 10,
    maxSize  = 500,
    BPPARAM  = BPPARAM
  )
  fg_null <- fgseaMultilevel(
    pathways = compass_sets,
    stats    = stat_r_o,
    minSize  = 10,
    maxSize  = 500,
    BPPARAM  = BPPARAM
  )
  
  # 8) Label & combine
  real_df <- fg_real %>%
    dplyr::select(pathway, pval, NES) %>%
    mutate(random_flag = 1)
  null_df <- fg_null %>%
    dplyr::select(pathway, pval, NES) %>%
    mutate(random_flag = 0)
  
  combined <- bind_rows(real_df, null_df) %>%
    mutate(
      GSE           = GSE,
      GENE          = GENE,
      Context       = params$Context[i],
      Perturbation  = perturbation,
      Samples       = gsms,
      n_ctrl        = sum(sml=="0"),
      n_treat       = sum(sml=="1")
    )
  
  final_results[[i]] <- combined
}

# ====== WRAP-UP ======
final_df    <- bind_rows(final_results)

head(final_df)

long_results <- final_df %>%
  dplyr::select(
    GSE, GENE, Context, Perturbation, pathway,
    pval, NES, random_flag, Samples, n_ctrl, n_treat
  )

# write out
write_csv(long_results,
          here::here("reproducibility", "figure3", "data", "250721 Benchmarking_COMPASS_across_contexts_fgsea.csv"))

message("Done: ", nrow(long_results), " total rows.")
