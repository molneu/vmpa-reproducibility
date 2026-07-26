# Benchmarking null‐model FGSEA & Collectri‐ULM for Figure 3d

# IMPORTANT!...needs to run STEP7 to STEP9 manually for each scrambled dataset

# Author: Igor Cima
# Created: 2025-07-27  Updated: 2025-07-28

# ===== STEP 1: Load libraries =====
library(readr)       # fast CSV I/O
library(dplyr)       # data manipulation (explicit dplyr::select())
library(purrr)       # functional programming helpers
library(fgsea)       # fast GSEA
library(decoupleR)    # network-based ULM
library(doParallel)  # parallel backend for foreach
library(foreach)     # foreach loops
library(cmapR)       # for compass_gsc()
library(GSEABase)    # GeneSetCollection
library(here)        # project-root paths

# ===== STEP 2: Define directories =====
base_dir    <- here("reproducibility","figure3")
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
subset_dir  <- file.path(data_dir, "subsets")

# ensure results directory exists
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

# ===== STEP 3: Load real-data FGSEA results =====
nes_file <- file.path(results_dir,
                      "250710 Benchmarking_fgsea results_ALL_CONTEXTS.csv")
nes_df <- readr::read_csv(nes_file, col_types = cols())

tested_sets <- nes_df %>%
  dplyr::select(Context, GENE, collection, new_id) %>%
  distinct()

# ===== STEP 4: Setup COMPASS buckets =====

# Compass function

compass_gsc <- function(context,
                        subset_dir   = subset_dir,
                        n            = 200,
                        min_conf     = 1,
                        targets      = NULL,
                        driver_filter = FALSE,
                        output       = c("list", "df", "gsc")) {
  output <- match.arg(output)
  require(cmapR)
  require(GSEABase)
  
  # 1) Load the context‐specific subsetted GCT object
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) {
    stop("No subset file found for context: ", context)
  }
  gct <- readRDS(subset_file)
  
  # 2) Select high‐confidence signatures by total CPS score
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  
  # 2b) (optional) filter by cancer‐driver annotation
  if (driver_filter) {
    cds <- gct@cdesc$cancer_driver_summary
    keep_idx <- keep_idx[ cds[keep_idx] != "None" ]
  }
  
  if (length(keep_idx) == 0) {
    if (output == "gsc") return(GeneSetCollection(list()))
    if (output == "df")  return(data.frame())
    return(list())
  }
  
  # 3) Extract matrix & reassign rownames to gene symbols
  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol
  
  # 4) Build output names: "<id>_c<TOTAL_CONF>"
  ids      <- gct@cdesc$id[keep_idx]
  conf_tot <- gct@cdesc$cps_conf_total[keep_idx]
  out_names <- paste0(ids, "_c", conf_tot)
  
  # 5) (Optional) filter to only specified targets
  if (!is.null(targets)) {
    keep_tgt <- gct@cdesc$cmap_name[keep_idx] %in% targets
    if (!any(keep_tgt)) {
      warning("No signatures found for targets: ", paste(targets, collapse = ", "))
      if (output == "gsc") return(GeneSetCollection(list()))
      if (output == "df")  return(data.frame())
      return(list())
    }
    mat       <- mat[, keep_tgt,       drop = FALSE]
    out_names <- out_names[keep_tgt]
  }
  
  # 6) Extract bottom N genes per signature
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    ord  <- suppressWarnings(order(vals, decreasing = FALSE, na.last = "keep"))
    head(rownames(mat)[ord], n)
  })
  names(res_list) <- out_names
  
  # 7) Return as GeneSetCollection
  if (output == "gsc") {
    gs <- lapply(out_names, function(setname) {
      GeneSet(geneIds = res_list[[setname]], setName = setname)
    })
    return(GeneSetCollection(gs))
  }
  
  # 8) Return as data.frame
  if (output == "df") {
    max_len <- max(lengths(res_list))
    df <- as.data.frame(
      do.call(cbind, lapply(res_list, function(vec) {
        length(vec) <- max_len; vec
      })),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(df)
  }
  
  # 9) Default: return as list
  res_list
}

# Setup buckets

compass_sizes <- c(25,50,100,150,200,250,300)
filtered_gs <- tested_sets %>%
  filter(collection %in% paste0("b", compass_sizes)) %>%
  group_split(Context, collection) %>%
  set_names(map_chr(., ~ paste0(.x$Context[1], "_", .x$collection[1]))) %>%
  map(function(df) {
    ctx <- df$Context[1]
    n   <- as.integer(sub("^b", "", df$collection[1]))
    all_sets <- compass_gsc(
      context    = ctx,
      subset_dir = subset_dir,
      n          = n,
      min_conf   = 1,
      output     = "list"
    )
    names(all_sets) <- sub("_c\\d+$", "", names(all_sets))
    all_sets[df$new_id]
  })

# ===== STEP 5: Pull & filter SigComLINCS =====
sigcom_files <- list(
  SigCom150 = file.path(data_dir, "240918 LINC_gene sets_CRISPR_CONSENSUS_SigComLINCS_bottom_150.csv"),
  SigCom200 = file.path(data_dir, "240918 LINC_gene sets_CRISPR_CONSENSUS_SigComLINCS_bottom_200.csv"),
  SigCom250 = file.path(data_dir, "240918 LINC_gene sets_CRISPR_CONSENSUS_SigComLINCS_bottom_250.csv")
)
sigcom_lists <- imap(sigcom_files, function(path, nm) {
  df <- readr::read_csv(path, col_types = cols())
  keep_cols <- intersect(names(df), tested_sets$GENE)
  df %>%
    dplyr::select(all_of(keep_cols)) %>%
    map(~ na.omit(as.character(.x))) %>%
    set_names(keep_cols)
})

# ===== STEP 6: Pull & filter Collectri regulons =====
collectri_file <- file.path(data_dir, "Collectri_genesets_forfGSEA.csv")
collectri_df   <- readr::read_csv(collectri_file, col_types = cols())
keep_coll <- tested_sets %>%
  filter(collection == "Collectri") %>%
  pull(GENE) %>%
  intersect(names(collectri_df), .)
collectri_list <- collectri_df %>%
  dplyr::select(all_of(keep_coll)) %>%
  map(~ na.omit(as.character(.x))) %>%
  set_names(keep_coll)

# Collectri network for ULM
t(net) <- get_collectri(organism="human", split_complexes=FALSE)
null_collections <- c(filtered_gs, sigcom_lists, list(Collectri = collectri_list))

# ===== STEP 7: Load scrambled T‐statistics =====

#run this for each of the following datasets:

#"240924 Benchmarking_R_list_25NULL MODELS FOR EACH PARAM.rds"
#"240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_2.rds"                                         
#"240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_3.rds"                                         
#"240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_4.rds" 

scrambled_file <- file.path(data_dir,
"240924 Benchmarking_R_list_25NULL MODELS FOR EACH PARAM.rds")

scrambled_T <- readRDS(scrambled_file)
gc()

# ===== STEP 8: Parallel FGSEA + Collectri‐ULM =====
# setup parallel backend
num_cores <- parallel::detectCores() - 1
cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)

# progress bar
total <- length(scrambled_T)
pb <- txtProgressBar(min=0, max=total, style=3)
count <- 0L
combine_fn <- function(x, y) {
  count <<- count + 1L
  setTxtProgressBar(pb, count)
  rbind(x, y)
}

# run
all_fgsea <- foreach(i = seq_along(scrambled_T), .combine=combine_fn,
                     .packages=c("fgsea","decoupleR","dplyr","purrr")) %dopar% {
                       scr <- scrambled_T[[i]]
                       g   <- scr$GENE[1]
                       GSE <- scr$GEO[1]
                       ctx <- scr$Context[1]
                       sid <- scr$ScrambleID[1]
                       
                       rv <- setNames(scr$stat, scr$Symbol)
                       rv <- rv[is.finite(rv) & !duplicated(names(rv))]
                       rv <- sort(rv, decreasing=TRUE)
                       
                       # build pathways list
                       gs_pathways <- list()
                       for(col in names(null_collections)) {
                         pathways <- null_collections[[col]]
                         hits <- grep(paste0("^",g,"$|^",g,"_"), names(pathways), value=TRUE)
                         for(h in hits) gs_pathways[[paste0(col, "::", h)]] <- pathways[[h]]
                       }
                       if (length(gs_pathways)==0) return(NULL)
                       gs_pathways <- purrr::map(gs_pathways, ~ .x[!is.na(.x)])
                       
                       # FGSEA real & null
                       fg <- fgsea(pathways = gs_pathways, stats=rv, minSize=10, maxSize=500, nperm=1e4)
                       fg$ScrambleID <- sid
                       
                       # Collectri ULM
                       net_sub <- net[net$source==g, ]
                       ulm <- tryCatch(
                         run_ulm(mat=rv, net=net_sub, .source="source", .target="target", .mor="mor", minsize=5),
                         error=function(e) data.frame(source=g, score=NA_real_, p_value=NA_real_, stringsAsFactors=FALSE)
                       )
                       ulm_df <- ulm %>%
                         select(pathway = source, NES = score) %>%
                         mutate(Collection="Collectri_ULM", GeneSetName=pathway,
                                TargetGene=g, GEO=GSE, Context=ctx, ScrambleID=sid)
                       
                       # annotate FGSEA output
                       parts <- strsplit(fg$pathway, "::", fixed=TRUE)
                       fg$Collection  <- vapply(parts, `[`, 1, FUN.VALUE="")
                       fg$GeneSetName <- vapply(parts, `[`, 2, FUN.VALUE="")
                       fg$TargetGene  <- g
                       fg$GEO         <- GSE
                       fg$Context     <- ctx
                       
                       fg_df <- fg %>%
                         select(pathway, NES, Collection, GeneSetName, TargetGene, GEO, Context, ScrambleID)
                       bind_rows(fg_df, ulm_df)
                     }

# teardown
close(pb)
stopCluster(cl)

# ===== STEP 9: Save full null FGSEA dataset =====

output_all <- file.path(results_dir, "250715_fGSEA_25NULLmodels_dataset1.csv")
readr::write_csv(all_fgsea, output_all)
message("Wrote full null FGSEA results to: ", output_all)

# ===== STEP 10: Summarize null models =====

# here we have one combined, so summarize directly

# 1) List your four null‐model CSVs
files <- paste0("250715_fGSEA_25NULLmodels_dataset", 1:4, ".csv")

# 2) Read and bind into one data.frame
all_null <- files %>%
  map_df(read_csv, col_types = cols())

# 3) Summarize: mean ± sd of NES, grouped by all the metadata
summary_null <- all_null %>%
  group_by(
    pathway,
    Collection,
    GeneSetName,
    TargetGene,
    GEO,
    Context
  ) %>%
  summarise(
    grand_mean_NES = mean(NES, na.rm = TRUE),
    grand_sd_NES   = sd(  NES, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  summary_null,
  here("reproducibility", "figure3", "data", "250715_fGSEA_100NULLmodels_SUMMARY.csv")
)
