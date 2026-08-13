#!/usr/bin/env Rscript
# Figure 5 Reproducibility Script

# ====== LIBRARIES ======
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)
library(stringr)
library(caret)
library(pbapply)
library(tidyr)
library(tibble)
library(purrr)
library(DescTools)
library(pheatmap)
library(M3C)
library(GSVA)
library(GSEABase)
library(Biobase)
library(corrplot)
library(readr)
library(RColorBrewer)
library(limma)
library(here)
library(httr)

# ====== SETUP ======
fig5_dir <- here("reproducibility", "figure5", "data")
results_dir <- here("reproducibility", "figure5", "results")
dir.create(fig5_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ====== FUNCTION: compass_gsc ======
compass_gsc <- function(context,
                        subset_dir   = here("reproducibility","figure3","data","subsets"),
                        n            = 200,
                        min_conf     = 1,
                        targets      = NULL,
                        driver_filter = FALSE,
                        output       = c("list", "df", "gsc")) {
  output <- match.arg(output)
  require(cmapR)
  require(GSEABase)
  
  # 1) Load context-specific subset
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) stop("No subset file found for context: ", context)
  gct <- readRDS(subset_file)
  
  # 2) Filter by confidence
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  if (driver_filter) {
    cds <- gct@cdesc$cancer_driver_summary
    keep_idx <- keep_idx[ cds[keep_idx] != "None" ]
  }
  if (length(keep_idx) == 0) {
    if (output == "gsc") return(GeneSetCollection(list()))
    if (output == "df")  return(data.frame())
    return(list())
  }
  
  # 3) Extract matrix & assign gene symbols
  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol
  
  # 4) Build names with confidence suffix
  ids <- gct@cdesc$id[keep_idx]
  conf_tot <- gct@cdesc$cps_conf_total[keep_idx]
  out_names <- paste0(ids, "_c", conf_tot)
  
  # 5) Optional target restriction
  if (!is.null(targets)) {
    keep_tgt <- gct@cdesc$cmap_name[keep_idx] %in% targets
    if (!any(keep_tgt)) {
      warning("No signatures found for targets: ", paste(targets, collapse = ", "))
      if (output == "gsc") return(GeneSetCollection(list()))
      if (output == "df")  return(data.frame())
      return(list())
    }
    mat <- mat[, keep_tgt, drop = FALSE]
    out_names <- out_names[keep_tgt]
  }
  
  # 6) Extract bottom n genes per signature
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    ord  <- suppressWarnings(order(vals, decreasing = FALSE, na.last = "keep"))
    head(rownames(mat)[ord], n)
  })
  names(res_list) <- out_names
  
  # 7) Format output
  if (output == "gsc") {
    gs <- lapply(out_names, function(setname) {
      GeneSet(geneIds = res_list[[setname]], setName = setname)
    })
    return(GeneSetCollection(gs))
  }
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
  res_list
}

# ====== STEP 1: Check data folder ======
file_names <- c(
  "250428 GSE174554_seurat_malignant.harmony_CLEAN_LAYERS JOINT_IDHwt_all.rds",
  "250626 GSE174554_pseudobulked transcriptomes_SCALED and feature selected.rds",
  "Neftel_TableS2_MES&NPCjoined2.csv"
)

missing_files <- file_names[!file.exists(file.path(fig5_dir, file_names))]
if (length(missing_files) > 0) {
  stop(
    "Missing Figure 5 input files in ", fig5_dir, ": ",
    paste(missing_files, collapse = ", "),
    ". Download the Figshare data and place the Figure 5 inputs in this data folder."
  )
}

# ====== STEP 2: Load & preprocess ======

# Load malignant single-nuclei Seurat object
seurat_malignant.clean.jl <- readRDS(file.path(fig5_dir,
                                               "250428 GSE174554_seurat_malignant.harmony_CLEAN_LAYERS JOINT_IDHwt_all.rds"))

# Load transcriptome-selected features (precomputed)
selected_data_list_tr <- readRDS(file.path(fig5_dir,
                                           "250626 GSE174554_pseudobulked transcriptomes_SCALED and feature selected.rds"))

# Load Neftel metaprograms
neftel_joined2 <- read.csv(file.path(fig5_dir,
                                     "Neftel_TableS2_MES&NPCjoined2.csv"))

# ====== STEP 2: Process single-nucleus data ======

# 1) Extract count matrix and metadata
counts <- seurat_malignant.clean.jl[["RNA"]]$counts
meta <- seurat_malignant.clean.jl@meta.data

# 2) Build Seurat object and filter samples with too few cells
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
seurat_clean <- CreateSeuratObject(counts = counts, meta.data = meta)
cell_counts <- table(seurat_clean$orig.ident)
valid_samples <- names(cell_counts[cell_counts >= 100])
seurat_clean_filtered <- subset(seurat_clean, subset = orig.ident %in% valid_samples)

rm(seurat_clean)
rm(seurat_malignant.clean.jl)

# 3) Split by sample
sample_list <- SplitObject(seurat_clean_filtered, split.by = "orig.ident")
gc()

# 4) Define processing function per sample
process_sample <- function(seurat_obj, normalization = "log") {
  tryCatch({
    if (normalization == "sct") {
      seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
      reduction_name <- "pca"
    } else {
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000, verbose = FALSE)
      seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
      seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
      reduction_name <- "pca"
    }
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    resolutions <- 1.5
    result_list <- list()
    for (res in resolutions) {
      temp_obj <- seurat_obj
      temp_obj <- FindNeighbors(temp_obj, reduction = reduction_name, dims = 1:30, verbose = FALSE)
      temp_obj <- FindClusters(temp_obj, resolution = res, verbose = FALSE)
      temp_obj <- RunUMAP(temp_obj, reduction = reduction_name, dims = 1:30, verbose = FALSE)
      temp_obj$normalization <- normalization
      temp_obj$resolution <- res
      result_list[[paste0(normalization, "_res", res)]] <- temp_obj
    }
    return(result_list)
  }, error = function(e) {
    warning("Failed to process sample: ", e$message)
    return(NULL)
  })
}

# 5) Process all samples (only log normalization here)
results <- pbapply::pblapply(sample_list, function(obj) {
  list(
    log = process_sample(obj, normalization = "log")
    # sct = process_sample(obj, normalization = "sct")
  )
})

# 6) Build cluster summary metadata frame
cluster_summary <- purrr::map_dfr(names(results), function(sample_name) {
  sample_result <- results[[sample_name]]
  purrr::map_dfr(names(sample_result), function(norm_method) {
    normalization_result <- sample_result[[norm_method]]
    purrr::map_dfr(names(normalization_result), function(res_name) {
      seurat_obj <- normalization_result[[res_name]]
      meta <- seurat_obj@meta.data %>%
        dplyr::mutate(
          sample = sample_name,
          normalization = norm_method,
          resolution = seurat_obj$resolution,
          cluster = seurat_clusters,
          cell_id = rownames(seurat_obj@meta.data)
        ) %>%
        dplyr::select(sample, normalization, resolution, cluster, cell_id)
      return(meta)
    })
  })
})

# 7) Summarize cell counts per grouping
cluster_summary_count <- cluster_summary %>%
  dplyr::count(sample, normalization, resolution, cluster, name = "n_cells") %>%
  dplyr::arrange(sample, normalization, resolution, cluster) %>%
  dplyr::mutate(unique_cluster_id = paste(sample, normalization, resolution, cluster, sep = "_"))

# ====== STEP 3: PSEUDOBULK ======

pseudobulk_with_aggregate <- function(seurat_obj, sample_id, normalization, resolution) {
  Seurat::Idents(seurat_obj) <- "seurat_clusters"
  cluster_counts <- table(Seurat::Idents(seurat_obj))
  if (length(cluster_counts) <= 1) {
    message("⚠️ Skipping pseudobulk: only one cluster found.")
    return(NULL)
  }
  pseudobulk <- AggregateExpression(
    seurat_obj,
    assay = "RNA",
    return.seurat = TRUE,
    group.by = "seurat_clusters",
    slot = "counts"
  )
  pseudobulk@meta.data$sample_id <- sample_id
  pseudobulk@meta.data$normalization <- normalization
  pseudobulk@meta.data$resolution <- resolution
  pseudobulk@meta.data$cluster <- rownames(pseudobulk@meta.data)
  pseudobulk@meta.data$sample_cluster <- paste0(sample_id, "_", pseudobulk@meta.data$cluster)
  return(pseudobulk)
}

# Build pseudobulked list
pseudobulked_list <- list()
for (sample_id in names(results)) {
  for (norm_method in names(results[[sample_id]])) {
    norm_list <- results[[sample_id]][[norm_method]]
    for (res_name in names(norm_list)) {
      seurat_obj <- norm_list[[res_name]]
      resolution <- seurat_obj$resolution[1]
      message("Processing: ", sample_id, " | ", norm_method, " | ", resolution)
      pb_obj <- pseudobulk_with_aggregate(seurat_obj, sample_id, norm_method, resolution)
      key <- paste(sample_id, norm_method, res_name, sep = "_")
      pseudobulked_list[[key]] <- pb_obj
    }
  }
}

# Filter valid clusters and drop original results to save memory
valid_clusters <- cluster_summary_count %>%
  dplyr::filter(n_cells >= 20) %>%
  dplyr::mutate(sample_cluster = paste0(sample, "_", cluster))
rm(results)

# Merge pseudobulks
merged_pseudobulk <- merge(pseudobulked_list[[1]], pseudobulked_list[-1])
rm(pseudobulked_list)

# Subset for resolution=1.5 & log normalization
sub_seu <- subset(merged_pseudobulk, subset = resolution == 1.5 & normalization == "log")

# JoinLayers for later expression extraction
merged_pseudobulk.jl <- JoinLayers(merged_pseudobulk)

# Add unique_cluster_id to metadata
meta <- merged_pseudobulk.jl@meta.data
meta$unique_cluster_id <- paste(
  meta$sample_id,
  meta$normalization,
  meta$resolution,
  meta$seurat_clusters,
  sep = "_"
)
merged_pseudobulk.jl <- AddMetaData(merged_pseudobulk.jl, metadata = meta$unique_cluster_id, col.name = "unique_cluster_id")

# Merge cluster counts into metadata
meta <- merged_pseudobulk.jl@meta.data
meta$unique_cluster_id_clean <- gsub("_g", "_", meta$unique_cluster_id)
cluster_summary_count$unique_cluster_id <- as.character(cluster_summary_count$unique_cluster_id)
meta <- meta %>%
  dplyr::left_join(cluster_summary_count %>% dplyr::select(unique_cluster_id, n_cells),
                   by = c("unique_cluster_id_clean" = "unique_cluster_id")) %>%
  dplyr::rename(total_cells = n_cells)
merged_pseudobulk.jl@meta.data <- meta

# Build final expression dataframe (normalized data)
count_matrix <- as.data.frame(GetAssayData(merged_pseudobulk.jl, slot = "data"))
meta <- merged_pseudobulk.jl@meta.data
new_colnames <- meta %>%
  dplyr::transmute(col_id = paste(sample_id, normalization, resolution, cluster, sep = "_")) %>%
  dplyr::pull()
colnames(count_matrix) <- new_colnames
final_gene_expr_df <- as.data.frame(count_matrix)
write.csv(final_gene_expr_df, file.path(fig5_dir, "final_gene_expr_df.csv"))

# ====== STEP 4: COMPASS ======

# Load final expression (pseudobulk) matrix
final_gene_expr_df <- read.csv(file.path(fig5_dir, "final_gene_expr_df.csv"), row.names = 1)

# Filter metadata for clusters with >=20 cells
filtered_meta <- meta %>% dplyr::filter(total_cells >= 20)

# Subset expression matrix accordingly
final_gene_expr_filtered <- final_gene_expr_df[, colnames(final_gene_expr_df) %in% filtered_meta$unique_cluster_id]

# Define normalization-resolution tag(s)
combinations <- "log_1.5"
subset_expr_list <- lapply(combinations, function(tag) {
  cols <- grep(tag, colnames(final_gene_expr_filtered), value = TRUE)
  final_gene_expr_filtered[, cols]
})
names(subset_expr_list) <- combinations

# Load COMPASS gene sets for glioma
DIR <- here("reproducibility", "figure3", "data", "subsets")
compass_sets_gs <- compass_gsc("glioma", subset_dir = DIR, n = 250, min_conf = 1, output = "gsc")

# Run GSVA per subset
gsva_results_list <- lapply(subset_expr_list, function(df) {
  data.m <- as.matrix(df)
  data.es <- ExpressionSet(assayData = data.m)
  gsvapar <- gsvaParam(data.es, compass_sets_gs)
  gene_sets <- geneSets(gsvapar)
  gsvapar <- gsvaParam(data.es, gene_sets)
  gsva_result <- gsva(gsvapar)
  as.data.frame(exprs(gsva_result))
})

saveRDS(gsva_results_list, file.path(results_dir, "STEP 4 output_gsva_results_list.rds"))

# ====== STEP 5: CONSENSUS CLUSTERING ======

# Prepare selected data from COMPASS (scaling & variance selection)
selected_data_list <- list()
for (dataset_name in names(gsva_results_list)) {
  dataset <- gsva_results_list[[dataset_name]]
  N <- nrow(dataset)
  scaled_data <- scale(dataset, scale = apply(dataset, 2, sd) * sqrt((N - 1) / N))
  cor_matrix <- cor(t(scaled_data))
  high_cor <- findCorrelation(cor_matrix, cutoff = 0.8)
  scaled_data_corr <- scaled_data[-high_cor, ]
  num_genes <- nrow(scaled_data_corr)
  top_10_percent_genes <- ceiling(num_genes * 0.10)
  top_variance_genes <- head(order(apply(scaled_data_corr, 1, var), decreasing = TRUE), top_10_percent_genes)
  selected_data <- scaled_data_corr[top_variance_genes, ]
  selected_data_list[[dataset_name]] <- selected_data
}

# Consensus clustering via M3C
perform_consensus_clustering <- function(expr_data) {
  M3C(
    mydata     = expr_data,
    cores      = 8,
    clusteralg = "spectral",
    method     = 1,
    seed       = 1,
    pacx1      = 0.05,
    pacx2      = 0.95,
    maxK       = 15,
    removeplots = TRUE,
    silent      = FALSE
  )
}

m3c_results_list.1 <- perform_consensus_clustering(selected_data_list[["log_1.5"]])
m3c_results_list_tr.1 <- perform_consensus_clustering(selected_data_list_tr[["log_1.5"]])

combined_m3c_results <- list(
  compass_log_1.5 = m3c_results_list.1,
  transcriptome_log_1.5 = m3c_results_list_tr.1
)
saveRDS(combined_m3c_results, file.path(results_dir, "STEP 5 output_M3C object_combined_m3c_results.rds"))

# (optional) Load Step 5 output
#combined_m3c_results <- readRDS(file.path(results_dir, "STEP 5 output_M3C object_combined_m3c_results.rds"))
#m3c_results_list.1 <- combined_m3c_results[1]
#m3c_results_list_tr.1 <- combined_m3c_results[2]

# ====== STEP 6: ANNOTATION TABLE ======



n_clu <- length(combined_m3c_results[["compass_log_1.5"]][["assignments"]])
annotation_col <- data.frame(Cluster = character(n_clu), stringsAsFactors = FALSE)

# Retrieve ordered data and assign cluster names
k <- 11
ordered_data <- m3c_results_list.1[["realdataresults"]][[k]][["ordered_data"]]
cluster_names <- colnames(ordered_data)
annotation_col$Cluster <- cluster_names

# M3C cluster info (COMPASS)
annotation_col$M3C_cluster <- factor(
  combined_m3c_results[["compass_log_1.5"]][["realdataresults"]][[k]][["ordered_annotation"]]$consensuscluster
)
annotation_col$M3C_name <- colnames(m3c_results_list.1[["realdataresults"]][[k]][["consensus_matrix"]])
rownames(annotation_col) <- annotation_col$Cluster

# Transcriptome M3C cluster
ordered_data.tr <- m3c_results_list_tr.1[["realdataresults"]][[11]][["ordered_data"]]
cluster_names.tr <- colnames(ordered_data.tr)
annotation_col$Cluster_TR <- cluster_names.tr
annotation_col$M3C_cluster_TR <- factor(
  combined_m3c_results[["transcriptome_log_1.5"]][["realdataresults"]][[11]][["ordered_annotation"]]$consensuscluster
)

# Patient ID / pair_or_id annotation
seurat_malignant.clean.jl <- readRDS(file.path(fig5_dir,
                                               "250428 GSE174554_seurat_malignant.harmony_CLEAN_LAYERS JOINT_IDHwt_all.rds"))
seurat_metadata <- seurat_malignant.clean.jl@meta.data
patient_ids <- sapply(strsplit(annotation_col$Cluster, "_"), `[`, 1)
patient_ids_clean <- sub("v[0-9]+$", "", patient_ids)
pair_lookup <- aggregate(pair ~ orig.ident, data = seurat_metadata, FUN = function(x) x[!is.na(x)][1])
matched_pairs <- pair_lookup$pair[match(patient_ids_clean, pair_lookup$orig.ident)]
pair_or_id <- ifelse(is.na(matched_pairs), patient_ids_clean, matched_pairs)
annotation_col$PatientID <- factor(pair_or_id)

# Sex and progression status
seurat_meta <- seurat_malignant.clean.jl@meta.data[, c("orig.ident", "sex", "progression_status")]
metadata_lookup <- seurat_meta %>%
  dplyr::distinct(orig.ident, sex, progression_status)
patient_ids <- sapply(strsplit(annotation_col$Cluster, "_"), `[`, 1)
annotation_col$sex <- metadata_lookup$sex[match(patient_ids, metadata_lookup$orig.ident)]
annotation_col$sex <- stringr::str_to_title(annotation_col$sex) %>% trimws()
annotation_col$progression_status <- metadata_lookup$progression_status[match(patient_ids, metadata_lookup$orig.ident)] %>% trimws()

# Neftel metaprograms processing
neftel_joined2 <- neftel_joined2[, -c(5,6)]  # remove G1S/G2M
neftel_joined2_l <- tidyr::pivot_longer(
  neftel_joined2,
  cols = everything(),
  names_to = "Metaprogram",
  values_to = "Gene",
  values_drop_na = TRUE
)
neftel.broad2 <- split(x = neftel_joined2_l$Gene, f = neftel_joined2_l$Metaprogram)

# GSVA Neftel on final expression (use saved final_gene_expr_df)
data.m <- as.matrix(final_gene_expr_df)
data.es <- ExpressionSet(assayData = data.m)
gsvapar.n1 <- gsvaParam(data.es, neftel.broad2)
gene_sets.n1 <- geneSets(gsvapar.n1)
gsvapar.n1 <- gsvaParam(data.es, gene_sets.n1)
gsva_results_neftel.broad2 <- gsva(gsvapar.n1)
gsva_neftel.broad2_scores <- exprs(gsva_results_neftel.broad2)

# Define subtype by max score
subtypes.n.broad2.gsva <- c("AC", "MES", "NPC", "OPC")
max_indices <- apply(gsva_neftel.broad2_scores, 2, which.max)
defined_subtype_n.broad2.gsva <- subtypes.n.broad2.gsva[max_indices]
gsva_neftel.broad2_scores <- rbind(gsva_neftel.broad2_scores, defined_subtype_n.broad2.gsva)
rownames(gsva_neftel.broad2_scores)[nrow(gsva_neftel.broad2_scores)] <- "Neftel_Metaprogram"

# Map Neftel scores into annotation_col
sample_clusters <- as.character(annotation_col$Cluster)
neftel_scores_by_cluster <- gsva_neftel.broad2_scores["Neftel_Metaprogram", , drop = TRUE]
neftel_values2 <- neftel_scores_by_cluster[sample_clusters]
annotation_col$Neftel_Metaprogram2 <- neftel_values2

# Save annotation table
write.csv(annotation_col,
          file.path(results_dir, "STEP 6 output_Suppl. Table_Annotation_TABLE_consensus_clustering.csv"),
          row.names = FALSE)

# ====== STEP 7: VISUALIZATIONS & ANALYSES ======

# a) Cluster / patient counts ----
cluster_per_patient_count <- annotation_col %>%
  as.data.frame() %>%
  dplyr::group_by(PatientID) %>%
  dplyr::summarise(n_M3C_clusters = n_distinct(M3C_cluster))
patient_per_cluster_count <- annotation_col %>%
  as.data.frame() %>%
  dplyr::group_by(M3C_cluster) %>%
  dplyr::summarise(n_patients = n_distinct(PatientID))

# Display median performance table
cluster_pat_tbl <- ggtexttable(patient_per_cluster_count, rows = NULL,
            theme = ttheme("minimal", base_size = 14,
                           colnames.style = list(face = "bold")))

ggplot2::ggsave(
  filename = here("reproducibility","figure5","results","Suppl. Fig., Patient per cluster count.pdf"),
  plot = cluster_pat_tbl,
  width = 4, height = 4, dpi = 300, device = "pdf"
)

# b) Fig 5b: RCSI dotplots ----
rcsi_values <- m3c_results_list.1[["scores"]][["RCSI"]]
p_values_log <- m3c_results_list.1[["scores"]][["P_SCORE"]]
k_values <- seq(2, length(rcsi_values) + 1)
plot_data <- data.frame(Clusters = factor(k_values), RCSI = rcsi_values, P_value = p_values_log)

rcsi_values_t <- m3c_results_list_tr.1[["scores"]][["RCSI"]]
p_values_log_t <- m3c_results_list_tr.1[["scores"]][["P_SCORE"]]
k_values_t <- seq(2, length(rcsi_values_t) + 1)
plot_data_t <- data.frame(Clusters = factor(k_values_t), RCSI = rcsi_values_t, P_value = p_values_log_t)

custom_theme <- theme_classic() +
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.title = element_blank())

rcsi_compass <- ggplot(plot_data, aes(x = RCSI, y = P_value, color = Clusters, label = Clusters)) +
  geom_point(size = 3) +
  geom_text(vjust = -1) +
  labs(x = "Relative Cluster Stability Index (RCSI)", y = "-log10(p-value)", title = "Consensus clustering") +
  scale_x_continuous(limits = c(-0.8, 0.6), breaks = seq(-0.6, 1, by = 0.2), labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, by = 2)) +
  custom_theme +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "none")

ggplot2::ggsave(
  filename = here("reproducibility","figure5","results","Fig 5b_RIGHT panel_rcsi_dotplot_COMPASS.pdf"),
  plot = rcsi_compass,
  width = 4, height = 5, dpi = 300, device = "pdf"
)

rcsi_transcriptome <- ggplot(plot_data_t, aes(x = RCSI, y = P_value, color = Clusters, label = Clusters)) +
  geom_point(size = 3) +
  geom_text(vjust = -1) +
  labs(x = "Relative Cluster Stability Index (RCSI)", y = "-log10(p-value)", title = "Consensus clustering") +
  scale_x_continuous(limits = c(-0.6, 1.2), breaks = seq(-0.6, 1.2, by = 0.2), labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, by = 2)) +
  custom_theme +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.position = "none")

ggplot2::ggsave(
  filename = here("reproducibility","figure5","results","Fig 5b_LEFT panel_rcsi_dotplot_TRANSCRIPTOME.pdf"),
  plot = rcsi_transcriptome,
  width = 4, height = 5, dpi = 300, device = "pdf"
)

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()

# c) Fig 5c: Consensus heatmaps ----
# Prepare consensus matrices
cons_mat <- m3c_results_list.1[["realdataresults"]][[k]][["consensus_matrix"]]
cons_mat_TR <- m3c_results_list_tr.1[["realdataresults"]][[k]][["consensus_matrix"]]

# Rename consensus matrix columns to user cluster labels
id_map <- setNames(annotation_col$Cluster, annotation_col$M3C_name)
colnames(cons_mat) <- id_map[colnames(cons_mat)]
rownames(cons_mat) <- colnames(cons_mat)
colnames(cons_mat_TR) <- id_map[colnames(cons_mat_TR)]
rownames(cons_mat_TR) <- colnames(cons_mat_TR)

# Annotation for heatmaps
anno_pheatmap <- annotation_col[, rev(c("M3C_cluster", "Neftel_Metaprogram2"))]
anno_pheatmap <- anno_pheatmap[colnames(cons_mat), , drop = FALSE]
rownames(anno_pheatmap) <- colnames(cons_mat)
cluster_colors2 <- c(
  "1" = "#2ca02c", "2" = "#FFDD44", "3" = "#FF7F0E", "4" = "#7f7f7f",
  "5" = "#1f77b4", "6" = "#D62728", "7" = "#9467bd", "8" = "#8c564b",
  "9" = "#17becf", "10" = "#e377c2", "11" = "black"
)
neftel_colors2 <- c("AC" = "#66c2a5", "MES" = "#a63a79", "NPC" = "#377eb8", "OPC" = "#ffae42")
ann_colors <- list(M3C_cluster = cluster_colors2, Neftel_Metaprogram2 = neftel_colors2)

# Save COMPASS heatmap
png(filename = file.path(results_dir, "Fig 5c_RIGHT panel_Consensus_heatmap_COMPASS.png"),
    width = 2000, height = 1600, res = 300)
pheatmap::pheatmap(cons_mat,
                   annotation_col = anno_pheatmap,
                   annotation_colors = ann_colors,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   legend = FALSE,
                   border_color = NA)

dev.off()

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()

# Left panel (transcriptome)
anno_pheatmap_tr <- annotation_col[, rev(c("M3C_cluster_TR", "Neftel_Metaprogram2"))]
anno_pheatmap_tr <- anno_pheatmap_tr[colnames(cons_mat_TR), , drop = FALSE]
rownames(anno_pheatmap_tr) <- colnames(cons_mat_TR)

png(filename = file.path(results_dir, "Fig 5c_LEFT panel_Consensus_heatmap_TRANSCRIPTOMES.png"),
    width = 2000, height = 1600, res = 300)
pheatmap::pheatmap(cons_mat_TR,
                   annotation_col = anno_pheatmap_tr,
                   annotation_colors = ann_colors,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   legend = FALSE,
                   border_color = NA)
dev.off()

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()


# d) Fig 5d: PAC / Gini / global consensus across K ----
k_values_vec <- 2:length(m3c_results_list.1[["realdataresults"]])
global_COMPASS <- numeric(); pac_COMPASS <- numeric(); gini_COMPASS <- numeric()
global_TR <- numeric(); pac_TR <- numeric(); gini_TR <- numeric()

for (kk in k_values_vec) {
  mat_COMPASS <- m3c_results_list.1[["realdataresults"]][[kk]][["consensus_matrix"]]
  upper_vals_COMPASS <- mat_COMPASS[upper.tri(mat_COMPASS)]
  global_COMPASS <- c(global_COMPASS, mean(upper_vals_COMPASS))
  pac_COMPASS <- c(pac_COMPASS, ecdf(upper_vals_COMPASS)(0.9) - ecdf(upper_vals_COMPASS)(0.1))
  gini_COMPASS <- c(gini_COMPASS, DescTools::Gini(as.vector(mat_COMPASS)))
  
  mat_TR <- m3c_results_list_tr.1[["realdataresults"]][[kk]][["consensus_matrix"]]
  upper_vals_TR <- mat_TR[upper.tri(mat_TR)]
  global_TR <- c(global_TR, mean(upper_vals_TR))
  pac_TR <- c(pac_TR, ecdf(upper_vals_TR)(0.9) - ecdf(upper_vals_TR)(0.1))
  gini_TR <- c(gini_TR, DescTools::Gini(as.vector(mat_TR)))
}

plot_df_indices <- data.frame(
  k = rep(k_values_vec, times = 2),
  GlobalConsensus = c(global_COMPASS, global_TR),
  PAC = c(pac_COMPASS, pac_TR),
  Gini = c(gini_COMPASS, gini_TR),
  Source = rep(c("COMPASS", "TRANSCRIPTOMES"), each = length(k_values_vec))
)

custom_theme2 <- theme_classic() +
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.title = element_blank())

p_global <- ggplot(plot_df_indices, aes(x = k, y = GlobalConsensus, color = Source)) +
  geom_line() + geom_point() +
  labs(title = "Global Consensus vs. Number of Clusters", x = "k", y = "Global Consensus (Average co-clustering)") +
  custom_theme2
p_pac <- ggplot(plot_df_indices, aes(x = k, y = PAC, color = Source)) +
  geom_line() + geom_point() +
  labs(title = "PAC vs. Number of Clusters", x = "k", y = "PAC") +
  custom_theme2
p_gini <- ggplot(plot_df_indices, aes(x = k, y = Gini, color = Source)) +
  geom_line() + geom_point() +
  labs(title = "Gini Index vs. Number of Clusters", x = "k", y = "Gini Index") +
  custom_theme2

ggplot2::ggsave(filename = file.path(results_dir, "Suppl. Fig., Global_consensus.pdf"), plot = p_global, width = 4, height = 5, dpi = 300)
ggplot2::ggsave(filename = file.path(results_dir, "Fig 5d_LEFT panel_pac_index.pdf"), plot = p_pac, width = 4, height = 5, dpi = 300)
ggplot2::ggsave(filename = file.path(results_dir, "Fig 5d_RIGHT panel_gini_index.pdf"), plot = p_gini, width = 4, height = 5, dpi = 300)

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()

# e) Fig 5e: Contingency analyses ----
data <- annotation_col
data <- subset(data, M3C_cluster != 7 & M3C_cluster != 4)
data$M3C_cluster <- factor(data$M3C_cluster, levels = c(2, 1, 8, 3, 5, 9, 11, 6, 10))
data$Neftel_Metaprogram2 <- factor(data$Neftel_Metaprogram2, levels = c("NPC", "OPC", "MES", "AC"))

cont.table <- table(data$M3C_cluster, data$Neftel_Metaprogram2)
chi_res <- chisq.test(cont.table)
residuals <- chi_res$stdres
residuals_thresholded <- residuals
residuals_thresholded[abs(residuals_thresholded) <= 1.96] <- NA

pdf(file = file.path(results_dir, "Fig 5e_Contingency_clusters_vs_metaprograms.pdf"),
    width = 6 , height = 6)
corrplot::corrplot(residuals_thresholded, is.cor = FALSE,
                   col = colorRampPalette(rev(RColorBrewer::brewer.pal(8, "RdBu")))(10),
                   na.label = " ",
                   tl.col = "black",
                   na.label.col = "white")
dev.off()

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()

cont.table.p <- table(data$M3C_cluster, data$progression_status)
data.s <- subset(data, sex != "None")
cont.table.s <- table(data.s$M3C_cluster, data.s$sex)

chi_res.p <- chisq.test(cont.table.p)
chi_res.s <- chisq.test(cont.table.s)
residuals.p <- chi_res.p$stdres
residuals.s <- chi_res.s$stdres

residuals_thresholded.p <- residuals.p
residuals_thresholded.p[abs(residuals_thresholded.p) <= 1.96] <- NA
residuals_thresholded.s <- residuals.s
residuals_thresholded.s[abs(residuals_thresholded.s) <= 1.96] <- NA

pdf(file = file.path(results_dir, "Suppl. Fig., Contingency_clusters_vs_progression_status.pdf"),
    width = 6, height = 6)
corrplot::corrplot(residuals_thresholded.p, is.cor = FALSE,
                   col = colorRampPalette(rev(RColorBrewer::brewer.pal(8, "RdBu")))(10),
                   na.label = " ",
                   tl.col = "black",
                   na.label.col = "white")
dev.off()

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()


pdf(file = file.path(results_dir, "Suppl. Fig., Contingency_clusters_vs_sex.pdf"),
    width = 6, height = 6)
corrplot::corrplot(residuals_thresholded.s, is.cor = FALSE,
                   col = colorRampPalette(rev(RColorBrewer::brewer.pal(8, "RdBu")))(10),
                   na.label = " ",
                   tl.col = "black",
                   na.label.col = "white")
dev.off()

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()


# f) Fig 5f: Marker protein heatmap ----
# --- 1. Get data
prot_mat <- selected_data_list[["log_1.5"]] # from STEP5a
meta_df   <- annotation_col # from STEP6

# --- 2. Filter metadata to exclude cluster 4 and relabel clusters
meta_df_filtered <- meta_df %>%
  filter(M3C_cluster != 7 & M3C_cluster != 4) %>%
  mutate(M3C_cluster = paste0("C", M3C_cluster))

# --- 3. Filter prot_mat and match/reorder metadata in one step
keep <- intersect(colnames(prot_mat), meta_df_filtered$Cluster)
prot_mat_filtered <- prot_mat[, keep]
meta_df_filtered  <- meta_df_filtered[match(keep, meta_df_filtered$Cluster), ]

# --- 4. Create column annotations for heatmap
annot_col.hm <- data.frame(
  ClusterLabel = meta_df_filtered$M3C_cluster,
  sample       = meta_df_filtered$Cluster,
  row.names    = meta_df_filtered$Cluster
)

# Sanity check
stopifnot(all(rownames(annot_col.hm) == colnames(prot_mat_filtered)))

# --- 5. Order samples by cluster  

#Specify order
cluster_order <- paste0("C", c(2, 1, 8, 3, 5, 9, 11, 6, 10))

annot_col.hm$ClusterLabel <- factor(annot_col.hm$ClusterLabel, levels = cluster_order)
ordered_samples <- rownames(annot_col.hm)[order(annot_col.hm$ClusterLabel)]
prot_mat_filtered <- prot_mat_filtered[, ordered_samples]
annot_col.hm        <- annot_col.hm[ordered_samples, , drop = FALSE]

# --- 6. Create design matrix using existing labels
cluster_labels <- annot_col.hm$ClusterLabel
design         <- model.matrix(~ 0 + cluster_labels)
colnames(design) <- levels(cluster_labels)

# --- 7. Create contrasts
clusters <- colnames(design)
contrast_list <- lapply(clusters, function(c) {
  others <- setdiff(clusters, c)
  paste0(c, " - (", paste(others, collapse = " + "), ")/", length(others))
})
names(contrast_list) <- paste0(clusters, "_vs_rest")

contrast_matrix <- makeContrasts(contrasts = unlist(contrast_list), levels = design)

# --- 8. Fit model and extract results
fit   <- lmFit(prot_mat_filtered, design)
fit2  <- contrasts.fit(fit, contrast_matrix)
fit2  <- eBayes(fit2)

all_results_df <- lapply(colnames(contrast_matrix), function(name) {
  topTable(fit2, coef = name, number = Inf, adjust.method = "fdr") %>%
    mutate(contrast = name, gene = rownames(.))
}) %>% bind_rows()

# --- 9. Select top 15 markers by t statistic per cluster
top_markers_per_cluster <- all_results_df %>%
  group_by(contrast) %>%
  slice_max(order_by = t, n = 15, with_ties = FALSE) %>%
  ungroup()

top_genes <- unique(top_markers_per_cluster$gene)

# --- 10. Subset and order proteins

top_markers_per_cluster <- top_markers_per_cluster %>%
  mutate(
    cluster = sub("^(C\\d+).*", "\\1", contrast),
    cluster = factor(cluster, levels = cluster_order)
  ) %>%
  arrange(cluster)

gene_order <- top_markers_per_cluster %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

gene_order <- gene_order[gene_order %in% rownames(prot_mat_filtered)]
prot_mat_top <- prot_mat_filtered[gene_order, , drop = FALSE]

# --- 11. Visualize heatmap
breaks <- seq(-3, 3, length.out = 50)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaks) - 1)

ann_colors <- list(
  ClusterLabel = c(
    "C1"  = "#2ca02c", "C2"  = "#FFDD44", "C3"  = "#FF7F0E",
    "C5"  = "#1f77b4", "C6"  = "#D62728",
    "C8"  = "#8c564b", "C9"  = "#17becf", "C10" = "#e377c2", "C11" = "black"
  )
)

top15_protein_activity_markers <- file.path(results_dir, "Fig 5f_Top_markers_heatmap_COMPASS.png")

# open a high-res raster device
png(
  filename = top15_protein_activity_markers,
  width    = 2400,    # pixels
  height   = 1500,    # pixels
  res      = 300      # dpi
)

pheatmap(
  prot_mat_top,
  annotation_col    = annot_col.hm["ClusterLabel"],
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  scale             = "row",
  color             = colors,
  breaks            = breaks,
  fontsize_row      = 3
)

dev.off()


Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()


# g) Fig 5g: EGFR/CDK4/hypoxia bubble plot ----
prot_mat2 <- gsva_results_list[["log_1.5"]]
meta_df2 <- annotation_col
ordered_data2 <- m3c_results_list.1[["realdataresults"]][[11]][["ordered_data"]]
meta_df2$Cluster <- colnames(ordered_data2)
meta_df2$M3C_cluster <- factor(
  combined_m3c_results[["compass_log_1.5"]][["realdataresults"]][[11]][["ordered_annotation"]]$consensuscluster
)
keep2 <- intersect(
  colnames(prot_mat2),
  meta_df2$Cluster[meta_df2$M3C_cluster != 7 & meta_df2$M3C_cluster != 4]
)
prot_mat_filtered2 <- prot_mat2[, keep2]
meta_df_filtered2 <- meta_df2 %>% dplyr::filter(Cluster %in% keep2)
N2 <- nrow(prot_mat_filtered2)
prot_mat_filtered2 <- scale(prot_mat_filtered2, scale = apply(prot_mat_filtered2, 2, sd) * sqrt((N2 - 1)/N2))

cluster_info <- meta_df_filtered2 %>%
  dplyr::select(sample = Cluster, cluster = M3C_cluster)

genes_of_interest <- list(EGFR = 1, CDK4 = 1, CA12 = 1)
selected_rows <- purrr::map_chr(
  names(genes_of_interest),
  ~ {
    matches <- grep(paste0("^", .x, "(_|$)"), rownames(prot_mat_filtered2), value = TRUE)
    idx <- genes_of_interest[[.x]]
    if (length(matches) >= idx) matches[idx] else NULL
  }
)
selected_rows <- purrr::compact(selected_rows)

expr_df <- prot_mat_filtered2[selected_rows, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(-feature, names_to = "sample", values_to = "expression") %>%
  dplyr::inner_join(cluster_info, by = "sample") %>%
  dplyr::mutate(gene = sub("(_.*)?$", "", feature))

plot_data_final <- expr_df %>%
  dplyr::group_by(cluster, gene) %>%
  dplyr::summarise(
    avg_expression = mean(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(z_score = as.numeric(scale(avg_expression))) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    expr_df %>%
      dplyr::count(cluster) %>%
      dplyr::rename(sample_count = n),
    by = "cluster"
  ) %>%
  dplyr::mutate(
    gene = factor(gene, levels = c("EGFR", "CDK4", "CA12")),
    cluster = factor(cluster, levels = rev(c(2, 1, 8, 3, 5, 9, 6, 10, 11)))
  )

ggplot(plot_data_final, aes(x = gene, y = cluster, size = sample_count, fill = z_score)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(5, 10), name = "Sample Count") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Z-Score") +
  labs(title = "Normalized Gene Expression Across Clusters", x = "Gene", y = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot2::ggsave(
  filename = here("reproducibility", "figure5", "results", "Fig 5g_Bubble Plot_CDK4_EGFR_hypoxia_bubble_COMPASS.pdf"),
  width = 6, height = 4, dpi = 300, device = "pdf"
)

Sys.sleep(0.1)
while (grDevices::dev.cur() > 1) grDevices::dev.off()


# Summary outputs ----
cell_counts <- as.data.frame(table(seurat_metadata$orig.ident))
message(paste0("samples, n = ", length(unique(seurat_metadata$orig.ident))))
paired_samples <- length(unique(na.omit(seurat_metadata$pair)))
unpaired_samples <- length(unique(seurat_metadata$orig.ident[is.na(seurat_metadata$pair)])) - 1
message(paste0("patients, n = ", paired_samples + unpaired_samples))
message(paste0("cells, n = ", length(seurat_metadata$orig.ident)))
message(paste0("pseudobulked clusters, n = ", n_clu))
message(paste0("pseudobulked clusters after removing consensus cluster 7 and 4, n = ", nrow(meta_df_filtered)))
