#!/usr/bin/env Rscript
# Precision–Recall Curves for COMPASS benchmarking (Figure 3d)
# Author: Igor Cima
# Created: 2025-07-28  Updated: 2025-07-29

# ===== STEP 1: Load libraries =====
library(readr)       # fast CSV I/O
library(dplyr)       # data manipulation (explicit dplyr::select())
library(tidyr)       # pivot_wider, pivot_longer
library(stringr)     # string operations
library(precrec)     # evalmod, auc
library(ggplot2)     # plotting
library(here)        # project-root paths

# ===== STEP 2: Define directories & input files =====
here::i_am("reproducibility/figure3/scripts/figure3d.R")

base_dir    <- here("reproducibility", "figure3")
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
# ensure results directory exists
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

nes_file  <- file.path(data_dir, "250710 Benchmarking_fgsea results_ALL_CONTEXTS.csv")  # real FGSEA results
null_file <- file.path(data_dir, "250715_fGSEA_100NULLmodels_SUMMARY.csv")                # null-model summary

# ===== STEP 3: Read data =====
nes_df       <- readr::read_csv(nes_file, col_types = cols())       # real FGSEA outputs
null_summary <- readr::read_csv(null_file, col_types = cols())      # null-model means & SDs

# ===== STEP 4: Match & combine with null summary =====
nes2 <- nes_df %>%
  mutate(
    join_collection = case_when(
      str_detect(collection, "^b\\d+") ~ paste0(Context, "_", collection),  # COMPASS buckets
      collection == "SigCom"            ~ "SigCom150",                       # SigCom key
      TRUE                                ~ collection                         # Collectri or ULM
    ),
    GeneSetName = case_when(
      str_detect(collection, "^b\\d+") ~ new_id,                             # COMPASS uses new_id
      TRUE                                 ~ GENE                                 # others use gene name
    )
  )

combined <- nes2 %>%
  left_join(
    null_summary %>%
      dplyr::rename(
        join_collection = Collection,
        TargetGene      = TargetGene,
        GEO             = GEO
      ),
    by = c(
      "join_collection",            # key in nes2
      "GSE"         = "GEO",      # align GSE <-> GEO
      "GENE"        = "TargetGene",
      "GeneSetName" = "GeneSetName"
    )
  )%>%

# compute z-score, guard against zero SD\ ncombined_z <- combined %>%
mutate(
  z_score = if_else(
    grand_sd_NES > 0,
    (NES - grand_mean_NES) / grand_sd_NES,
    NA_real_
  )
)

# ===== STEP 5: Annotate each gene‐set with its confidence level =====
suppressWarnings(
  combined_conf <- combined %>%
  dplyr::mutate(
    conf = as.integer(
      sub(".*_c(\\d+)$", "\\1", pathway.x)   # extract the N from “…_cN”
    )
  )
)

# ===== STEP 6: Select only the highest‐confidence sets per GENE & Context =====
combined_best <- combined_conf %>%
  dplyr::group_by(GENE, Context.x) %>%
  dplyr::mutate(
    max_conf = max(conf, na.rm = TRUE)  # returns -Inf if all conf are NA
  ) %>%
  dplyr::filter(
    is.na(conf) | conf == max_conf     # keep NAs or highest‐conf rows
  ) %>%
  dplyr::select(-max_conf) %>%         # clean up
  dplyr::ungroup()

# optional: write merged data for record
readr::write_csv(
  combined_best,
  file.path(results_dir, "Source data_Fig 3d_zscores.csv")
)

# ===== STEP 5: Prepare data for PR curves =====
df <- combined_best %>%
  mutate(
    z_adj = if_else(expected_nes_sign == "pos", -z_score, z_score)  # adjust for direction
  ) %>%
  filter(! collection %in% c("SigCom", "Collectri", "Collectri_ULM"))  # drop non-COMPASS contexts

# pivot wide for COMPASS buckets
wide_tbl <- df %>%
  dplyr::select(GSE, random_flag, pathway.x, collection, z_adj) %>%
  distinct() %>%
  tidyr::pivot_wider(
    id_cols      = c(GSE, random_flag, pathway.x),
    names_from   = collection,
    values_from  = z_adj,
    names_prefix = "zscore_"
  )

# ===== STEP 6: Build PR curve inputs for model =====

# COMPASS bucket sizes to include
compass_sizes <- c(25,50,100,150,200,250,300)
keep_cols     <- paste0("zscore_b", compass_sizes)
modnames_orig <- keep_cols

pr_input_orig <- wide_tbl %>%
  mutate(across(dplyr::starts_with("zscore_"), ~ - .x)) %>%  # flip sign for PRC
  dplyr::mutate(Label = as.numeric(random_flag == 1))        # 1 = real, 0 = null

# split into lists for evalmod
score_list_orig <- pr_input_orig %>% dplyr::select(all_of(keep_cols)) %>% as.list()
label_list_orig <- rep(list(pr_input_orig$Label), length(score_list_orig))

# ===== STEP 7: Run PRC for model =====

precrec_obj_orig <- evalmod(
scores   = score_list_orig,
labels   = label_list_orig,
modnames = modnames_orig,
mode     = "rocprc"    # compute both ROC and PRC
)

# ===== STEP 8: Plot PR curves  =====

msdf_orig <- fortify(precrec_obj_orig)

p1 <- ggplot(msdf_orig %>% filter(curvetype == "PRC"), aes(x = x, y = y, color = modname)) +
  geom_line(size = 1) +
  labs(
    title = "Unselected COMPASS gene sets",
    x     = "Recall", y = "Precision"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    zscore_b25  = "#0072B2",
    zscore_b50  = "#D55E00",
    zscore_b100 = "#009E73",
    zscore_b150 = "#D62728",
    zscore_b200 = "#CC79A7",
    zscore_b250 = "#56B4E9",
    zscore_b300 = "#E69F00"
  )) +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(results_dir, "Fig 3d_compass_PR_gene_set_size.pdf"),
  p1, width=5, height=4, bg="white"
)

# ===== STEP 9: Print AUC  =====

auc_df_orig <- auc(precrec_obj_orig)
message("AUCs for Precision-recall curves (COMPASS model):")
auc_df_orig %>% filter(curvetypes == "PRC") %>% print()

# pick the row with the max AUC
best_row <- auc_df_orig %>%
  dplyr::slice_max(aucs, n = 1)
                  
message("Best gene set size: ", best_row$modnames)

# ===== STEP 10: Prepare b250 benchmarking dataset =====

# rebuild wide with GENE to isolate b250

df <- combined_best %>%
  mutate(
    z_adj = if_else(expected_nes_sign == "pos", -z_score, z_score)  # adjust for direction
  )

# 2) Rebuild wide to isolate b250
wide_b250 <- df %>%
  dplyr::select(GSE, GENE, random_flag, pathway.x, collection, z_adj) %>%
  distinct() %>%
  tidyr::pivot_wider(
    id_cols      = c(GSE, GENE, random_flag, pathway.x),
    names_from   = collection,
    values_from  = z_adj,
    names_prefix = "zscore_"
  )

b250_long <- wide_b250 %>%
  dplyr::select(GSE, GENE, random_flag, pathway.x, zscore_b250) %>%
  dplyr::rename(pathway = pathway.x, z = zscore_b250)

# 3) Build best_b250_rows using pre‐computed high‐confidence selections
best_b250_rows <- df %>%
  # only b250 rows
  dplyr::filter(collection == "b250") %>%
  # keep only those pathways you marked high‐confidence
  dplyr::semi_join(
    combined_best %>% dplyr::select(GSE, GENE, pathway.x),
    by = c("GSE", "GENE", "pathway.x")
  ) %>%
  # join SigCom
  dplyr::left_join(
    df %>% dplyr::filter(collection=="SigCom") %>%
      dplyr::select(GSE, GENE, random_flag, z_adj) %>%
      dplyr::rename(z_SigCom = z_adj),
    by = c("GSE","GENE","random_flag")
  ) %>%
  # join Collectri
  dplyr::left_join(
    df %>% dplyr::filter(collection=="Collectri") %>%
      dplyr::select(GSE, GENE, random_flag, z_adj) %>%
      dplyr::rename(z_Collectri = z_adj),
    by = c("GSE","GENE","random_flag")
  ) %>%
  # join Collectri_ULM
  dplyr::left_join(
    df %>% dplyr::filter(collection=="Collectri_ULM") %>%
      dplyr::select(GSE, GENE, random_flag, z_adj) %>%
      dplyr::rename(z_Collectri_ULM = z_adj),
    by = c("GSE","GENE","random_flag")
  ) %>%
  # compute per‐GSE&GENE mean z for b250
  dplyr::left_join(
    b250_long %>%
      dplyr::group_by(GSE, GENE, random_flag) %>%
      dplyr::summarise(z_b250mean = mean(z, na.rm=TRUE), .groups="drop"),
    by = c("GSE","GENE","random_flag")
  ) %>%
  # flip all for PRC
  dplyr::mutate(across(dplyr::starts_with("z_"), ~ - .x))

# 4) Now run PRC exactly as before
keep_cols2 <- c("z_b250mean","z_SigCom","z_Collectri","z_Collectri_ULM")
modnames2  <- c("b250_mean","SigCom","Collectri","Collectri_ULM")
raw_scores2 <- best_b250_rows %>% dplyr::select(all_of(keep_cols2)) %>% as.list()
raw_labels2 <- best_b250_rows$random_flag

precrec_obj2 <- precrec::evalmod(
  scores   = raw_scores2,
  labels   = raw_labels2,
  modnames = modnames2,
  mode     = "rocprc"
)

# ===== STEP 12: Plot PR curves for best-b250 =====

msdf2 <- fortify(precrec_obj2)

p2 <- ggplot(msdf2 %>% filter(curvetype == "PRC"), aes(x=x,y=y,color=modname)) +
  geom_line(size=1) +
  labs(title="PR: Best b250 vs others", x="Recall", y="Precision") +
  theme_minimal() +
  scale_color_manual(values=c(
    b250_best      ="#0072B2",
    b250_mean      ="#56B4E9",
    SigCom         ="#D55E00",
    Collectri      ="#CC79A7",
    Collectri_ULM  ="#E69F00"
  )) +
  theme(legend.title=element_blank(), panel.border=element_rect(color="black",fill=NA,size=0.5), axis.text.x=element_text(angle=45,hjust=1))

ggsave(
  file.path(results_dir, "Fig 3d_compass_PR_benchmark.pdf"),
  p2, width=5, height=4, bg="white"
)


# ===== STEP 13: Print AUC for best-b250 =====

auc_df2 <- auc(precrec_obj2)
message("AUCs for Precision-recall curves (benchmarking):")
auc_df2 %>% filter(curvetypes=="PRC") %>% print()

