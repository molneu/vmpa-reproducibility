# Supplementary metadata table for all contexts
library(dplyr)
library(here)

# contexts you want to cover
contexts <- c("glioma","melanoma","nsclc","gastric","ovarian",
              "crc","breast","prostate","pdac","headneck")

# directory where subset .rds live
subset_dir <- here("reproducibility", "figure3", "data", "subsets")

# collect metadata
metadata_list <- lapply(contexts, function(ctx) {
  gct <- readRDS(file.path(subset_dir, paste0(ctx, "_subset.rds")))
  
  # gct@cdesc is the metadata per gene set
  meta <- gct@cdesc %>%
    as.data.frame() %>%
    mutate(Context = ctx)
  
  return(meta)
})

# combine all contexts
metadata_df <- bind_rows(metadata_list)

# Select, rename, and filter metadata table
metadata_filtered <- metadata_df %>%
  dplyr::filter(cps_conf_total >= 0) %>%   
  dplyr::select(
    cmap_name,
    context,
    id,
    sig_id,
    cancer_driver_summary,
    tas,
    cps_conf_total
  ) %>%
  dplyr::rename(
    Symbol               = cmap_name,
    Context              = context,
    Compass_id           = id,
    Tas              = tas,
    CMap_id               = sig_id,
    Cancer_driver_status = cancer_driver_summary,
    Confidence_score     = cps_conf_total
  )

length(unique(metadata_filtered$Symbol))

# Save as supplementary source data
results_dir <- here("reproducibility", "figure3", "results")
write.csv(metadata_filtered,
          file = file.path(results_dir, "Suppl. Table 1_compass_catalogue_metadata.csv"),
          row.names = FALSE)
