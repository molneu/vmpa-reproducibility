#!/usr/bin/env Rscript

# ===== LIBRARIES =====
library(readr)
library(dplyr)
library(purrr)
library(fgsea)
library(decoupleR)
library(doParallel)
library(cmapR)    # for compass_gsc()
library(GSEABase) # if you need GeneSetCollection
library(here)

# ===== 1) Load real‐data FGSEA results =====
nes_df <- read_csv(
  here("reproducibility", "figure3", "data", "250721 Benchmarking_COMPASS_across_contexts_fgsea.csv"),
  col_types = cols()
)

library(dplyr)
library(purrr)

# 0) Start from your nes_df
#    pull out a distinct list of (Context, new_id=pathway)
tested_sets <- nes_df %>%
  dplyr::select(Context, new_id = pathway) %>%
  distinct()

# 1) Only b250
compass_size <- 250
subset_dir    <- here("reproducibility", "figure3", "data", "subsets")

# 2) For each context, fetch only the b250 bucket and keep exactly those gene‐sets you tested
library(purrr)

filtered_gs <- split(tested_sets, tested_sets$Context) %>% 
  imap(function(df, ctx) {
    gs_all <- compass_gsc(
      context    = ctx,
      subset_dir = subset_dir,
      n          = 250,
      min_conf   = 1,
      output     = "list"
    )
    # pick only the ones you tested, may yield NULLs
    out <- gs_all[df$new_id]
    # drop any NULL elements (e.g. missing in gs_all)
    compact(out)
  }) %>% 
  # now drop any contexts that ended up with zero gene‐sets
  keep(~ length(.x) > 0)

# you’ll now have a list named like “glioma”, “melanoma”, etc.,
# each element is a named list of your b250 pathways (no empties).
# 3) That is now your null_collections (COMPASS‐only, size‐250)
null_collections <- filtered_gs

cat("Built COMPASS‐only b250 collections:\n")
print(names(null_collections))

# 4) Save for later
saveRDS(
  null_collections,
  here("reproducibility", "figure3", "data", "250725_null_collections_COMPASS_b250.rds")
)

# ===== 7) Load filtered null‐collections =====
null_collections <- readRDS(
  here("reproducibility", "figure3", "data", "250725_null_collections_COMPASS_b250.rds")
)

# ===== 8) Load scrambled T‐statistics =====

#"240924 Benchmarking_R_list_25NULL MODELS FOR EACH PARAM.rds"
#"240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_2.rds"                                         
#"240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_3.rds"                                         
#"240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_4.rds"    

scrambled_T <- readRDS(here("reproducibility", "figure3", "data", "240925_Benchmarking_R_list_25NULL_MODELS_FOR_EACH_PARAM_4.rds"))

gc()

# ===== 9) Parallel FGSEA  =====
# —— 0) flatten all size-250 COMPASS collections into one list ——  
all_compass250 <- purrr::flatten(filtered_gs)

# —— 1) parallel setup ——  
num_cores <- 4  
rscript   <- file.path(R.home("bin"), "Rscript")  
cl <- makePSOCKcluster(rep("localhost", num_cores), rscript = rscript)  
doParallel::registerDoParallel(cl)

total <- length(scrambled_T)  
pb    <- txtProgressBar(min = 0, max = total, style = 3)  
count <- 0L  
comb_fun2 <- function(a, b) {  
  count <<- count + 1L  
  setTxtProgressBar(pb, count)  
  bind_rows(a, b)  
}

start_time <- Sys.time()

all_fgsea <- foreach(i = seq_along(scrambled_T),
                     .combine  = comb_fun2,
                     .packages = c("fgsea","dplyr","purrr")) %dopar% {
                       
                       scr        <- scrambled_T[[i]]
                       g          <- scr[["GENE"]][1]
                       GSE        <- scr[["GEO"]][1]
                       ctx        <- scr[["Context"]][1]
                       scrambleID <- scr[["ScrambleID"]][1]
                       
                       # 1) ranking vector
                       rv <- setNames(scr[["stat"]], scr[["Symbol"]])
                       rv <- rv[is.finite(rv) & !duplicated(names(rv))]
                       rv <- sort(rv, decreasing = TRUE)
                       
                       # 2) filter to only those size-250 COMPASS sets for this gene
                       sel <- grep(paste0("^", g, "_"), names(all_compass250))
                       if (length(sel) == 0) return(NULL)
                       my_pathways <- all_compass250[sel]
                       
                       # 3) strip NAs & attrs
                       my_pathways <- map(my_pathways, function(x) {
                         x <- x[!is.na(x)]
                         attributes(x) <- NULL
                         x
                       })
                       
                       # 4) fgsea
                       set.seed(123)
                       fg <- fgsea(
                         pathways = my_pathways,
                         stats    = rv,
                         minSize  = 10,
                         maxSize  = 500,
                         nperm    = 5000
                       )
                       
                       # 5) annotate
                       fg %>%
                         mutate(
                           Collection   = "b250",
                           GeneSetName  = pathway,
                           TargetGene   = g,
                           GEO          = GSE,
                           Context      = ctx,
                           ScrambleID   = scrambleID
                         ) %>%
                         dplyr::select(
                           pathway, NES, pval, padj,
                           Collection, GeneSetName,
                           TargetGene, GEO, Context, ScrambleID
                         )
                     }

# —— teardown ——  
end_time <- Sys.time()  
message("Elapsed time: ", signif(end_time - start_time, 3))  
close(pb)  
stopCluster(cl)

# all_fgsea now has only those b250 gene-sets whose prefix matches each scramble’s GENE  
# `all_fgsea` now holds FGSEA results (NES, pval, padj) for every size-250 COMPASS gene‐set  

### ===== result 
# `all_fgsea` is now a data.frame with columns:
#   pathway, NES, Collection, GeneSetName, TargetGene, GEO, Context, ScrambleID

library(dplyr)
library(stringr)

all_fgsea <- all_fgsea %>%
  # 1) pull out the cell‐line code from anywhere in "pathway"
  mutate(
    GeneSetCell = str_extract(
      pathway,
      "U251MG|A375|A549|ES2|HT29|MCF7|PC3|YAPC"
    ),
    # 2) map that to your disease context
    GeneSetContext = recode(
      GeneSetCell,
      U251MG = "glioma",
      A375   = "melanoma",
      A549   = "nsclc",
      ES2    = "ovarian",
      HT29   = "crc",
      MCF7   = "breast",
      PC3    = "prostate",
      YAPC   = "pdac",
      .default = NA_character_
    )
  ) %>%
  # 3) keep only the columns you care about
  dplyr::select(
    pathway, NES, pval, padj,
    Collection, GeneSetName,
    TargetGene, GEO, Context, ScrambleID,
    GeneSetCell, GeneSetContext
  )

write.csv(
  all_fgsea,
  here("reproducibility", "figure3", "data", "250725_FIG3E_COMPARISON_fGSEA_25NULLmodels_dataset4.csv"),
  row.names = FALSE
)

# ===== 10) Combine and sumarize =====

# 1) List your four null‐model CSVs
files <- here("reproducibility", "figure3", "data", paste0("250725_FIG3E_COMPARISON_fGSEA_25NULLmodels_dataset", 1:4, ".csv"))

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
    GEO_Context = Context,
    GeneSet_Context = GeneSetContext
  ) %>%
  summarise(
    grand_mean_NES = mean(NES, na.rm = TRUE),
    grand_sd_NES   = sd(  NES, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  summary_null,
  here("reproducibility", "figure3", "data", "250725_FIG3E_COMPARISON_fGSEA_100NULLmodels_SUMMARY.csv")
)



