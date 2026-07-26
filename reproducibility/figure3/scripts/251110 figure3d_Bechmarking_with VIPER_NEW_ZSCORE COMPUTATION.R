#!/usr/bin/env Rscript
# ============================================================
# Merge experimental NES data with VIPER + fGSEA null model summaries
# Compute Z-scores and add standardized collection_PR column for PR analysis
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ====== File paths ======
results_dir <- here("reproducibility", "figure3", "results")
data_dir    <- here("reproducibility", "figure3", "data")

exp_file   <- file.path(results_dir, "251110_CLEAN_Benchmarking_VIPER_COMPASSb250_ALL_CONTEXTS_intersectionVIPER_2modes.csv")
viper_file <- file.path(results_dir, "251109_VIPER_100NullModels_NES_summary.csv")
fgsea_file <- file.path(data_dir, "250715_fGSEA_100NULLmodels_SUMMARY.csv")

# ====== Load data ======
exp_data    <- read_csv(exp_file)
viper_nulls <- read_csv(viper_file)
fgsea_nulls <- read_csv(fgsea_file)

# ============================================================
# STEP 1: Merge VIPER null statistics
# ============================================================

merged_viper <- exp_data %>%
  dplyr::filter(grepl("^VIPER_", collection)) %>%
  dplyr::left_join(
    viper_nulls %>%
      dplyr::select(GSE, GENE, Context, collection, grand_mean_NES, grand_sd_NES),
    by = c("GSE", "GENE", "Context", "collection")
  )

# ============================================================
# STEP 2: Normalize context names for fGSEA/COMPASS/SigCom/Collectri
# ============================================================

normalize_context <- function(x) {
  x <- trimws(x)
  dplyr::recode(x,
                "CRC" = "crc",
                "Glioblastoma" = "glioma",
                "PDAC" = "pdac",
                "Breast Cancer" = "breast",
                "Ovarian Cancer" = "ovarian",
                "Melanoma" = "melanoma",
                "Prostate Cancer" = "prostate",
                "NSCLC" = "nsclc",
                .default = tolower(x)
  )
}

# Clean both datasets
exp_data_clean <- exp_data %>%
  dplyr::mutate(
    Context = normalize_context(Context),
    GENE = toupper(GENE)
  )

fgsea_nulls_clean <- fgsea_nulls %>%
  dplyr::rename(GSE = GEO, collection = Collection, GENE = TargetGene) %>%
  dplyr::mutate(
    Context = normalize_context(Context),
    GENE = toupper(GENE)
  )

# ============================================================
# STEP 3: Handle COMPASS "_cX" suffix in GeneSetName
# ============================================================

exp_data_clean <- exp_data_clean %>%
  dplyr::mutate(
    GeneSetName = dplyr::case_when(
      grepl("_b[0-9]+$", collection) ~ stringr::str_remove(pathway, "_c[0-9]+$"),
      TRUE ~ NA_character_
    )
  )

# ============================================================
# STEP 4: Merge fGSEA/COMPASS/SigCom/Collectri null statistics
# ============================================================

# Separate joins by method type
merged_compass <- exp_data_clean %>%
  dplyr::filter(grepl("_b[0-9]+$", collection)) %>%  # COMPASS collections only
  dplyr::left_join(
    fgsea_nulls_clean %>%
      dplyr::select(GSE, GENE, Context, collection, GeneSetName,
                    grand_mean_NES, grand_sd_NES),
    by = c("GSE", "GENE", "Context", "collection", "GeneSetName")
  )

merged_sigcollect <- exp_data_clean %>%
  dplyr::filter(collection %in% c("SigCom250", "Collectri_ULM")) %>%
  dplyr::left_join(
    fgsea_nulls_clean %>%
      dplyr::select(GSE, GENE, Context, collection,
                    grand_mean_NES, grand_sd_NES),
    by = c("GSE", "GENE", "Context", "collection")
  )

# Combine the two parts for non-VIPER methods
merged_fgsea <- dplyr::bind_rows(merged_compass, merged_sigcollect)

# ============================================================
# STEP 5: Combine all merged results and compute Z-scores
# ============================================================

merged_all <- dplyr::bind_rows(merged_viper, merged_fgsea) %>%
  dplyr::mutate(
    z_score = dplyr::if_else(
      !is.na(grand_sd_NES) & grand_sd_NES > 0,
      (NES - grand_mean_NES) / grand_sd_NES,
      NA_real_
    )
  )

# ============================================================
# STEP 6: Add unified collection_PR for PR/ROC analyses
# ============================================================

merged_all <- merged_all %>%
  dplyr::mutate(
    collection_PR = dplyr::case_when(
      grepl("_b[0-9]+$", collection) ~ paste0("COMPASS_", stringr::str_extract(collection, "b[0-9]+$") %>%
                                                stringr::str_remove("^b")),
      TRUE ~ collection
    )
  )

# ============================================================
# STEP 7: Save combined result
# ============================================================

output_file <- file.path(results_dir, "251111_BENCHMARKING_withNullModels_Zscores.csv")
readr::write_csv(merged_all, output_file)

message("✅ Done: wrote ", nrow(merged_all), " rows to ", output_file)
# ============================================================
# STEP 8: Quick summary check
# ============================================================

summary_check <- merged_all %>%
  dplyr::group_by(collection_PR) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_zscore = sum(!is.na(z_score)),
    prop_zscore = n_zscore / n_total,
    .groups = "drop"
  )

print(summary_check)

# ============================================================
# STEP 9: Verify duplicate structure (random_flag)
# ============================================================

check_pairs <- exp_data_clean %>%
  dplyr::count(GSE, GENE, Context, collection, random_flag) %>%
  dplyr::count(GSE, GENE, Context, collection) %>%
  dplyr::filter(n > 2)

if (nrow(check_pairs) == 0) {
  message("✅ All duplicates explained by random_flag (expected).")
} else {
  message("⚠️ Found unexpected many-to-many duplicates beyond random_flag:")
  print(check_pairs)
}

