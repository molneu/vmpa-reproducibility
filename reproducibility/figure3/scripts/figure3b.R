# Figure 3 Pipeline: CMap LV5 ➔ Context‑specific LINC gene sets
# Author: Igor Cima
# Created: 2025-06-24  
# Updated: 2025-07-27

here::i_am("reproducibility/figure3/scripts/figure3b.R")

# Load required libraries -----------------------------------------------------
library(cmapR)    
library(dplyr)     
library(purrr)     
library(here) 
library(ggplot2)

# ====== STEP 1: Setup directories and paths ======
# Define output directories for figure 3 reproducibility
fig3_data_dir    <- here("reproducibility", "figure3", "data")
fig3_results_dir <- here("reproducibility", "figure3", "results")
shared_cmap_dir  <- here("reproducibility", "figure1", "data")

# Create directories if they do not exist (recursive, suppress warnings)
dir.create(fig3_data_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(fig3_results_dir, recursive = TRUE, showWarnings = FALSE)

# ====== STEP 2: Load & annotate column metadata ======
# 2.1 Define input files
gctx_file    <- file.path(shared_cmap_dir, "level5_beta_trt_xpr_n142901x12328.gctx")
siginfo_file <- file.path(shared_cmap_dir, "siginfo_beta.txt")

# 2.2 Parse GCTX and load signature info
trt_xpr   <- parse_gctx(gctx_file)                      # expression matrix
meta_full <- read.delim(siginfo_file, stringsAsFactors = FALSE)  # metadata table

# 2.3 Filter for transcriptional perturbations (pert_type == "trt_xpr")
meta_xpr  <- meta_full %>% filter(pert_type == "trt_xpr")

# 2.4 Annotate GCT object columns by matching on sig_id
annotated <- annotate_gct(trt_xpr, meta_xpr, dim="column", keyfield="sig_id")

# Clean up large intermediates
tmp_remove <- c("trt_xpr", "meta_full")
rm(list = tmp_remove); rm(tmp_remove)

# ====== STEP 3: Annotate row metadata (gene symbols & inference) ======
# 3.1 Load LINC gene annotation map
gene_map <- read.csv(
  file.path(shared_cmap_dir, "LINC gene annotations.csv"),
  stringsAsFactors = FALSE
)

# 3.2 Build lookup vectors
entrez2sym <- setNames(gene_map$gene.name, gene_map$entrez)
infer_map  <- setNames(gene_map$inference, gene_map$gene.name)

# 3.3 Annotate each row in the GCT object
rids <- annotated@rid  # row IDs (Entrez IDs)
annotated@rdesc <- annotated@rdesc %>%
  mutate(
    symbol    = unname(entrez2sym[rids]),
    inference = unname(infer_map[symbol])
  )

# ====== STEP 4: Annotate cancer‐driver status by perturbation gene ======

# Read in annotation CSVs
dressler_role_df   <- read.csv(file.path(fig3_data_dir, "250717_original_annotations_Dressler_role.csv"),   stringsAsFactors = FALSE)
dressler_status_df <- read.csv(file.path(fig3_data_dir, "250717_original_annotations_Dressler_status.csv"), stringsAsFactors = FALSE)
kinn_role_df       <- read.csv(file.path(fig3_data_dir, "250717_original_annotations_Kinnersley_role.csv"), stringsAsFactors = FALSE)

# Build named lookup vectors by symbol
map_dressler_role   <- setNames(dressler_role_df$Dressler_Role,       dressler_role_df$SYMBOL)
map_dressler_status <- setNames(dressler_status_df$Dressler_annotation, dressler_status_df$SYMBOL)
map_kinn_role       <- setNames(kinn_role_df$Kinnersley_Role,      kinn_role_df$SYMBOL)

# 4.2 Initialize new raw columns in column metadata
cd <- annotated@cdesc$cmap_name
annotated@cdesc <- annotated@cdesc %>%
  mutate(
    D_status_raw = map_dressler_status[cd],  # Dressler status
    D_role_raw   = map_dressler_role[cd],    # Dressler role
    K_role_raw   = map_kinn_role[cd]         # Kinnersley role
  )

# 4.3 Replace missing annotations with "None"
for (col in c("D_status_raw","D_role_raw","K_role_raw")) {
  annotated@cdesc[[col]][is.na(annotated@cdesc[[col]])] <- "None"
}

# 4.4 Recode raw annotations into concise labels
annotated@cdesc <- annotated@cdesc %>%
  mutate(
    D_status = case_when(
      D_status_raw %in% c("canonical cancer","canonical cancer and healthy") ~ "canonical",
      D_status_raw %in% c("candidate cancer","candidate cancer and healthy") ~ "candidate",
      TRUE                                                                    ~ "None"
    ),
    D_role = case_when(
      D_role_raw   == "Oncogene"              ~ "O",
      D_role_raw   == "Dual Role"             ~ "DR",
      D_role_raw   == "Tumour Suppressor"     ~ "TS",
      D_role_raw   %in% c("Unclassified","Conflicting annotations") ~ "Uncl",
      TRUE                                      ~ ""
    ),
    K_role = case_when(
      K_role_raw   == "Oncogene"        ~ "O",
      K_role_raw   == "Dual role"       ~ "DR",
      K_role_raw   == "Tumor Suppressor"~ "TS",
      TRUE                              ~ ""
    )
  )

# 4.5 Create a summary string combining all driver codes per signature
annotated@cdesc$cancer_driver_summary <- apply(
  annotated@cdesc[ , c("D_status","D_role","K_role")],
  1,
  function(codes) {
    parts <- character()
    if (codes["D_status"] != "None") parts <- c(parts, paste0("D_status::", codes["D_status"]))
    if (codes["D_role"]   != "")     parts <- c(parts, paste0("D_role::",   codes["D_role"]))
    if (codes["K_role"]   != "")     parts <- c(parts, paste0("K_role::",   codes["K_role"]))
    if (length(parts) == 0) "None" else paste(parts, collapse = "|")
  }
)

# ====== STEP 5: Filter signatures & annotate IDs ======
# 5.1 Filter for Activity Score (TAS) > 0.2
a02 <- subset_gct(annotated, cid = which(annotated@cdesc$tas > 0.2))
rm(annotated); gc()

# 5.2 Remove NA or blank cmap_name
a02 <- subset_gct(a02, cid = which(!is.na(a02@cdesc$cmap_name) & a02@cdesc$cmap_name != ""))

# 5.3 Drop control or border signatures (UnTrt, BRDN0...)
remove_idx <- which(
  a02@cdesc$cmap_name == "UnTrt" |
    grepl("^BRDN0", a02@cdesc$cmap_name)
)
if (length(remove_idx) > 0) {
  keep_idx <- setdiff(seq_along(a02@cdesc$cmap_name), remove_idx)
  a02       <- subset_gct(a02, cid = keep_idx)
}

# 5.4 Assign unique identifier column: "<cmap_name>_<sig_id>"
a02@cdesc$id <- paste0(a02@cdesc$cmap_name, "_", a02@cdesc$sig_id)

# ====== STEP 6: Split into context-specific subsets ======
# 6.1 Identify unique cell lines (contexts)
cell_types <- unique(meta_xpr$cell_iname)

# 6.2 Create a named list of GCT subsets per cell type
subsets <- map(cell_types, function(ct) {
  idx  <- which(a02@cdesc$cell_iname == ct)
  gct2 <- subset_gct(a02, cid = idx)
  gct2@cdesc <- a02@cdesc[idx, , drop = FALSE]  # preserve metadata
  gct2
})
names(subsets) <- cell_types

# 6.3 Keep only contexts with > 1000 signatures
subsets <- keep(subsets, ~ length(slot(.x, "cid")) > 1000)

# ====== STEP 7: Inject organ-of-origin context labels ======
# 7.1 Define mapping from cell line to organ context
context_map <- c(
  U251MG = "glioma", A375 = "melanoma", A549 = "nsclc",
  AGS    = "gastric", ES2  = "ovarian", HT29 = "crc",
  MCF7   = "breast", PC3  = "prostate", YAPC = "pdac",
  BICR6  = "headneck"
)

# 7.2 Apply context label to each subset's metadata
for (ct in names(subsets)) {
  cdesc <- subsets[[ct]]@cdesc
  cdesc$context <- context_map[cdesc$cell_iname]
  subsets[[ct]]@cdesc <- cdesc
}

# ====== STEP 8: Compute signature confidence scores per context ======
for(ct in names(subsets)){
  gct   <- subsets[[ct]]
  mat   <- as.data.frame(gct@mat)
  cdesc <- gct@cdesc
  
  # 8.1 Exemplar score
  exemplar <- as.integer(cdesc$is_exemplar_sig==1)
  names(exemplar) <- colnames(mat)
  
  # 8.2 Inference score
  inf <- setNames(numeric(ncol(mat)),colnames(mat))
  for(sig in colnames(mat)){
    ri <- match(sig,cdesc$sig_id)
    cn <- trimws(cdesc$cmap_name[ri])
    if(!is.na(cn)&&nzchar(cn)){
      gi <- match(cn,trimws(gct@rdesc$symbol))
      if(!is.na(gi)){
        val      <- mat[gi,sig]
        inf_type <- gct@rdesc$inference[gi]
        inf[sig] <- case_when(
          val>0 & inf_type=="landmark"     ~ -2,
          val>0 & inf_type=="best inferred"~ -1,
          val>0 & inf_type=="inferred"     ~ 0,
          val<0                              ~ 1,
          TRUE                               ~ 0
        )
      }
    }
  }
  
  # 8.3 Correlation score
  groups     <- split(colnames(mat),cdesc$cmap_name)
  corr_score <- setNames(rep(0,ncol(mat)),colnames(mat))
  for(grp in names(groups)){
    cols      <- groups[[grp]]
    if(length(cols)<2) next
    submat    <- cor(mat[,cols],use="pairwise.complete.obs")
    any_strong <- any(submat[upper.tri(submat)]>0.2,na.rm=TRUE)
    for(col_i in cols){
      others         <- setdiff(cols,col_i)
      vals_i         <- submat[col_i,others]
      score1         <- as.integer(any(vals_i>0.2,na.rm=TRUE))
      score2         <- if(score1==0&&any_strong) -1 else 0
      corr_score[col_i] <- score1+score2
    }
  }
  
  # 8.4 Total confidence
  total <- exemplar + inf + corr_score
  
  # 8.5 Attach scores to metadata by matching on sig_id
  cdesc <- cdesc %>% mutate(
    cps_conf_exemplar    = exemplar[sig_id],
    cps_conf_inference   = inf[sig_id],
    cps_conf_correlation = corr_score[sig_id],
    cps_conf_total       = total[sig_id]
  )
  subsets[[ct]]@cdesc <- cdesc
}

# ====== STEP 9: Save context subsets ======
#for (ct in names(subsets)) {
#  saveRDS(
#    subsets[[ct]],
#    file = file.path(fig3_data_dir, paste0(ct, "_subset.rds"))
#  )
#}

#read subsets

#fig3_data_dir <- here::here("reproducibility","figure3","data","subsets")
# subset_files <- list.files(
#  path      = fig3_data_dir,
#  pattern   = "_subset\\.rds$",
#  full.names= TRUE)
#subsets <- lapply(subset_files, readRDS)
#names(subsets) <- basename(subset_files) %>%
#  sub("_subset\\.rds$", "", .)


# ====== STEP 10: Generate Fig 3b barplot ======
# 10.1 Drop contexts without sufficient samples
keep_contexts <- setdiff(names(subsets), c("AGS", "BICR6"))
subs2 <- subsets[keep_contexts]

# 10.2 Count high-confidence signatures and unique perturbed genes per context
sig_counts    <- map_int(subs2, ~ sum(.x@cdesc$cps_conf_total >= 0, na.rm = TRUE))
unique_counts <- map_int(subs2, ~ .x@cdesc$cmap_name[.x@cdesc$cps_conf_total >= 0] %>% unique() %>% length())

# 10.3 Assemble summary tibble and reshape for plotting
gene_sets <- tibble(
  Context           = keep_contexts,
  n_gene_sets       = sig_counts,
  n_unique_proteins = unique_counts
) %>% arrange(desc(n_gene_sets))

melted_data <- gene_sets %>%
  tidyr::pivot_longer(
    cols      = -Context,         # every column except Context
    names_to  = "Signature",      # name of the column that will hold former column names
    values_to = "Gene"            # name of the column that will hold the values
  )

# 10.4 Set factor levels for consistent ordering
melted_data$Context  <- factor(melted_data$Context, levels = rev(gene_sets$Context))
melted_data$Signature <- factor(melted_data$Signature, levels = c("n_unique_proteins","n_gene_sets"))

write.csv(
  melted_data,
  file = file.path(fig3_results_dir, "Source Data Fig3b.csv"),
  row.names = FALSE
)


# 10.5 Plot horizontal bar chart of counts per context
final_plot <- ggplot(melted_data, aes(x = Context, y = Gene, fill = Signature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("n_gene_sets" = "#3E4A89",
                               "n_unique_proteins" = "#44AA99")) +
  coord_flip() +
  labs(x = "Context", y = "Count", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position   = "top",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks        = element_line(size = 0.3),
    axis.text.x       = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig3_results_dir, "Fig 3b_compass_gene_sets.pdf"),
  final_plot, width=3, height=5, bg="white"
)


# End of Figure 3b pipeline
