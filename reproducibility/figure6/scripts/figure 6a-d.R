#!/usr/bin/env Rscript

# ====== LIBRARIES ======
library(GSEABase)  # Gene set container
library(GSVA)      # GSVA/gsvaParam/gsva
library(limma)     # linear modeling for differential expression
library(ggplot2)   # plotting
library(here)      # project-rooted paths

# ====== PATHS & INPUTS ======
# ExpressionSet for Figure 6 is expected under reproducibility/figure6/data
eset_path <- here::here("reproducibility", "figure6", "data", "ExpressionSet_GSE145128.rds")
eset <- readRDS(eset_path)

# ====== QC: Expression distribution ======
# Boxplot of raw expression values (no outliers shown)
box <- boxplot(eset@assayData[["exprs"]], outline = FALSE, main = "Raw expression distribution")

# ====== COMPASS GENE SETS FUNCTION ======
# Extract gene signatures for a given context (copied/adapted from earlier)
compass_gsc <- function(context,
                        subset_dir   = here::here("reproducibility", "figure3", "data", "subsets"),
                        n            = 200,
                        min_conf     = 1,
                        targets      = NULL,
                        driver_filter = FALSE,
                        output       = c("list", "df", "gsc")) {
  output <- match.arg(output)
  require(cmapR)
  require(GSEABase)
  
  # Load context-specific subsetted GCT object
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) {
    stop("No subset file found for context: ", context)
  }
  gct <- readRDS(subset_file)
  
  # Filter by total CPS confidence
  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  
  # Optionally restrict to cancer-driver-associated signatures
  if (driver_filter) {
    cds <- gct@cdesc$cancer_driver_summary
    keep_idx <- keep_idx[ cds[keep_idx] != "None" ]
  }
  
  if (length(keep_idx) == 0) {
    if (output == "gsc") return(GeneSetCollection(list()))
    if (output == "df")  return(data.frame())
    return(list())
  }
  
  # Extract matrix and set gene symbols as rownames
  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol
  
  # Build signature names including confidence suffix
  ids      <- gct@cdesc$id[keep_idx]
  conf_tot <- gct@cdesc$cps_conf_total[keep_idx]
  out_names <- paste0(ids, "_c", conf_tot)
  
  # Optionally subset to specific targets
  if (!is.null(targets)) {
    keep_tgt <- gct@cdesc$cmap_name[keep_idx] %in% targets
    if (!any(keep_tgt)) {
      warning("No signatures found for targets: ", paste(targets, collapse = ", "))
      if (output == "gsc") return(GeneSetCollection(list()))
      if (output == "df")  return(data.frame())
      return(list())
    }
    mat       <- mat[, keep_tgt, drop = FALSE]
    out_names <- out_names[keep_tgt]
  }
  
  # Get bottom n genes per signature (low-valued genes)
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    ord  <- suppressWarnings(order(vals, decreasing = FALSE, na.last = "keep"))
    head(rownames(mat)[ord], n)
  })
  names(res_list) <- out_names
  
  # Format output
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

# ====== BUILD COMPASS GENE SETS (DRIVERS ONLY) ======
# Use glioma context, restrict to driver-associated signatures
compass_sets_drivers <- compass_gsc(
  context = "glioma",
  n = 250,
  output = "gsc",
  min_conf = 1,
  driver_filter = TRUE
)

# ====== GSVA ======
# Build GSVA parameter object and run GSVA on the ExpressionSet
gsvapar <- gsvaParam(eset, compass_sets_drivers)
gene_sets <- geneSets(gsvapar)          # process the gene sets to deduplicate etc.
gsvapar <- gsvaParam(eset, gene_sets)    # reinitialize with cleaned gene sets
gsva_result <- gsva(gsvapar)            # run GSVA
gsva_data <- as.data.frame(exprs(gsva_result))  # enrichment scores matrix

# ====== CLINICAL FILTERING ======
# Keep only cR and N paired samples: drop experimental relapsed samples (eR) and specific patients
gsva_result_clinical <- gsva_result[, gsva_result$relapse_TYPE != "eR"]
gsva_result_clinical <- gsva_result_clinical[, !(gsva_result_clinical$patient_ID %in% c("B046", "B078"))]

# Check distribution after clinical filtering
box2 <- boxplot(gsva_result_clinical@assayData[["exprs"]], main = "GSVA scores after clinical filtering")

# ====== NORMALIZATION ======
# Extract matrix to normalize by population SD
tonorm_eset <- gsva_result_clinical@assayData[["exprs"]]
N <- nrow(tonorm_eset)  # number of genes/features
# Scale each column (sample) using pooled/population standard deviation
norm_eset <- scale(tonorm_eset, scale = apply(tonorm_eset, 2, sd) * sqrt(N - 1 / N))
box3 <- boxplot(norm_eset, main = "Normalized GSVA scores")

# Update the ExpressionSet with normalized values
gsva_result_clinical.n <- gsva_result_clinical
exprs(gsva_result_clinical.n) <- norm_eset

# ====== DIFFERENTIAL PROTEIN ACTIVITY ANALYSIS ======
# Build design matrix for relapse type comparison
relapse_factor <- factor(gsva_result_clinical@phenoData@data[["relapse_TYPE"]])
# Relevel so that reference is the second level (preserves original intent)
relapse_factor_releveled <- relevel(relapse_factor, ref = levels(relapse_factor)[2])
mod <- model.matrix(~ relapse_factor_releveled)

# Fit linear model and compute statistics
fit <- lmFit(gsva_result_clinical.n, mod)
# Apply treat with logFC threshold of 0.05
fit <- eBayes(fit)
res <- decideTests(fit)

# Get full top table for coefficient 2 (contrast of interest)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "t")

tt <- topTable(fit, coef = 2, number = Inf) |>
  arrange(dplyr::desc(t))

# Save result as supplementary table to results folder
out_fname <- here::here("reproducibility", "figure6", "results",
                        "Suppl_Table_250717_GSE145128_naive_vs_cR_Compass_U251MGc1-3_b250_DRIVERS_ONLY_gsva_scaled.csv")
dir.create(dirname(out_fname), recursive = TRUE, showWarnings = FALSE)
write.csv(tt, out_fname, row.names = TRUE)

# Figure 6b STRING source data are stored under reproducibility/figure6/data.
# "Source data_Fig6b_string_annotations.txt" contains the 35-protein STRING input and annotations.
# "Source data_Fig6b_string_output.pdf" records the STRING v12.0 network and settings: full evidence
# network, experiments/databases, text mining excluded, minimum interaction score 0.400,
# no added interactors, and disconnected nodes hidden.
# "Source data_Fig6b_string_interactions.csv" contains the quantitative STRING interaction scores.

# ====== VOLCANO PLOT ======
# Prepare volcano plot annotations
tt$logP <- -log10(tt$P.Value)
tt$color_group <- "NS"
tt$color_group[tt$logFC > 0.02 & tt$logP > 1.3]  <- "Up"
tt$color_group[tt$logFC < -0.02 & tt$logP > 1.3] <- "Down"
tt$color_group <- factor(tt$color_group, levels = c("Down", "NS", "Up"))

# Define publication-friendly color scheme
volcano_colors <- c(
  Down = "#4575B4",
  NS   = adjustcolor("gray70", alpha.f = 0.1),
  Up   = "#D73027"
)

# Plot volcano
volcano_plot <- ggplot(tt, aes(x = logFC, y = logP, color = color_group)) +
  geom_point(alpha = 0.4, size = 2) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-0.02, 0.02), color = "gray50", linetype = "dashed") +
  geom_hline(yintercept = 1.3, color = "gray50", linetype = "dashed") +
  labs(
    x = "GSVA enrichment score difference",
    y = expression(-log[10]~~Raw~P-value),
    color = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

# Optionally save volcano plot
ggplot2::ggsave(
  filename = here::here("reproducibility", "figure6", "results", "Fig6a_Volcano_naive_vs_cR.pdf"),
  plot = volcano_plot,
  width = 4, height = 4, dpi = 300
)

# ====== BETA-CATENIN TRANSCRIPTOME vs PROTEIN ACTIVITY ======

# --- Prepare clinical subset of the expression Eset (already have this earlier) ---
# Keep only cR and N paired samples (remove experimental relapsed samples)
eset_clinical <- eset[, eset$relapse_TYPE != "eR"]
eset_clinical <- eset_clinical[, eset_clinical$patient_ID != "B046" & eset_clinical$patient_ID != "B078"]

# --- 1. Transcript-level CTNNB1 expression ----
# Extract CTNNB1 expression (single gene)
ctnnb1_expr <- exprs(eset_clinical)["CTNNB1", , drop = TRUE]

# Build phenotype / sample metadata from the expression eset
pheno_expr <- pData(eset_clinical)[, c("relapse_TYPE", "patient_ID", "geo_accession"), drop = FALSE]
pheno_expr_df <- as.data.frame(pheno_expr) %>%
  tibble::rownames_to_column("sample_id") %>%        # capture original sample identifier if needed
  dplyr::select(sample_id, relapse_TYPE, patient_ID, geo_accession)

# Combine into transcriptome data.frame
ctnnb1_expr_df <- data.frame(
  CTNNB1_expression = ctnnb1_expr,
  relapse_TYPE      = pheno_expr_df$relapse_TYPE,
  patientID         = pheno_expr_df$patient_ID,
  geo_ACCESSION     = pheno_expr_df$geo_accession,
  row.names         = NULL
)
# Optionally keep sample_id if you want explicit mapping:
ctnnb1_expr_df$sample_id <- pheno_expr_df$sample_id

# --- 2. Protein activity for the specific CTNNB1 gene set ----
target_set <- "CTNNB1_XPR041_U251MG.311_96H:O21_c2"

# Extract its activity vector
ctnnb1_act_vec <- gsva_result_clinical.n@assayData[["exprs"]][target_set, , drop = TRUE]

# Build phenotype data.frame from the GSVA clinical object
pheno_activity <- pData(gsva_result_clinical.n)[, c("relapse_TYPE", "patient_ID", "geo_accession"), drop = FALSE]
pheno_activity_df <- as.data.frame(pheno_activity) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, relapse_TYPE, patient_ID, geo_accession)

# Combine into protein activity data.frame
ctnnb1_act_df <- data.frame(
  CTNNB1_activity = ctnnb1_act_vec,
  relapse_TYPE    = pheno_activity_df$relapse_TYPE,
  patientID       = pheno_activity_df$patient_ID,
  geo_ACCESSION   = pheno_activity_df$geo_accession,
  row.names       = NULL
)
ctnnb1_act_df$sample_id <- pheno_activity_df$sample_id

# Define output directory (adjust if needed)
out_dir <- here::here("reproducibility", "figure6", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save CTNNB1 transcriptome expression dataframe
write.csv(
  ctnnb1_expr_df,
  file = file.path(out_dir, "Source_data_CTNNB1_transcript_expression.csv"),
  row.names = FALSE
)

# Save CTNNB1 protein activity dataframe
write.csv(
  ctnnb1_act_df,
  file = file.path(out_dir, "Source_data_CTNNB1_protein_activity.csv"),
  row.names = FALSE
)

library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Prepare transcriptome data ----
ctnnb1_df_plot <- ctnnb1_expr_df %>%
  dplyr::filter(relapse_TYPE %in% c("n", "cR")) %>%
  dplyr::rename(value = CTNNB1_expression) %>%
  dplyr::mutate(
    modality = "transcript",
    relapse_TYPE = factor(relapse_TYPE, levels = c("n", "cR"))
  ) %>%
  dplyr::select(patientID, relapse_TYPE, value, modality)

# ---- Prepare protein activity data ----
ctnnb1_act_df_plot <- ctnnb1_act_df %>%
  dplyr::filter(relapse_TYPE %in% c("n", "cR")) %>%
  dplyr::rename(value = CTNNB1_activity) %>% # adjust name if needed
  dplyr::mutate(
    modality = "protein_activity",
    relapse_TYPE = factor(relapse_TYPE, levels = c("n", "cR"))
  ) %>%
  dplyr::select(patientID, relapse_TYPE, value, modality)

# ---- Compute panel-specific ranges ----
expr_mat_transcriptome <- exprs(eset_clinical)
transcript_range <- range(expr_mat_transcriptome, na.rm = TRUE)

activity_mat <- gsva_result_clinical.n@assayData[["exprs"]]  # rows=gene sets, cols=samples
activity_range <- range(activity_mat, na.rm = TRUE)

transcript_limits <- c(transcript_range[1], transcript_range[2])
activity_limits   <- c(activity_range[1],   activity_range[2])


# ---- Plot transcript panel ----
p_transcript <- ggplot(ctnnb1_df_plot, aes(x = relapse_TYPE, y = value, group = patientID)) +
  geom_line(alpha = 0.6, color = "black") +
  geom_point(shape = 21, fill = NA, size = 3, stroke = 1, color = "black") +
  labs(x = "Relapse Type", y = "CTNNB1 transcript", title = "Transcript") +
  scale_y_continuous(limits = transcript_limits) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

# ---- Plot protein activity panel ----
p_activity <- ggplot(ctnnb1_act_df_plot, aes(x = relapse_TYPE, y = value, group = patientID)) +
  geom_line(alpha = 0.6, color = "black") +
  geom_point(shape = 21, fill = NA, size = 3, stroke = 1, color = "black") +
  labs(x = "Relapse Type", y = "CTNNB1 activity", title = "Protein activity") +
  scale_y_continuous(limits = activity_limits) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

# ---- Combine and save ----
combined <- p_transcript + p_activity + patchwork::plot_layout(nrow = 1)

ggplot2::ggsave(
  filename = here::here("reproducibility", "figure6", "results", "Fig6c_CTNNB1_paired_dotplot.pdf"),
  plot = combined,
  width = 3.5,
  height = 4,
  dpi = 300
)

# ---- Paired ttest ----

library(dplyr)

# Helper to run paired t-test given a dataframe with columns: patientID, relapse_TYPE, value
paired_ttest <- function(df, value_name) {
  # Pivot wider so each patient has both values side by side
  wide <- df %>%
    dplyr::filter(relapse_TYPE %in% c("n", "cR")) %>%
    tidyr::pivot_wider(
      id_cols = patientID,
      names_from = relapse_TYPE,
      values_from = value
    ) %>%
    # keep only patients with both measurements
    dplyr::filter(!is.na(n) & !is.na(cR))
  
  if (nrow(wide) < 2) {
    warning(glue::glue("Not enough paired samples for {value_name}"))
    return(NULL)
  }
  
  test <- t.test(wide$cR, wide$n, paired = TRUE)
  
  # Build summary
  summary_df <- tibble::tibble(
    modality        = value_name,
    n_pairs         = nrow(wide),
    t_statistic     = test$statistic,
    df              = test$parameter,
    p_value         = test$p.value,
    mean_diff       = test$estimate,          # cR - n
    conf_low        = test$conf.int[1],
    conf_high       = test$conf.int[2]
  )
  list(wide = wide, test = test, summary = summary_df)
}

# Prepare the two data frames (assumes these exist):
# ctnnb1_expr_df: transcript with columns CTNNB1_expression, relapse_TYPE, patientID, geo_ACCESSION
# ctnnb1_act_df: activity with columns CTNNB1_activity, relapse_TYPE, patientID, geo_ACCESSION

# Normalize column names for convenience
transcript_plot_df <- ctnnb1_expr_df %>%
  dplyr::rename(value = CTNNB1_expression,
                relapse_TYPE = relapse_TYPE,
                patientID = patientID) %>%
  dplyr::mutate(relapse_TYPE = factor(relapse_TYPE, levels = c("n", "cR")))

activity_plot_df <- ctnnb1_act_df %>%
  dplyr::rename(value = CTNNB1_activity,
                relapse_TYPE = relapse_TYPE,
                patientID = patientID) %>%
  dplyr::mutate(relapse_TYPE = factor(relapse_TYPE, levels = c("n", "cR")))

# Run paired t-tests
transcript_res <- paired_ttest(transcript_plot_df, "transcript")
activity_res   <- paired_ttest(activity_plot_df,   "protein_activity")

# Combine summaries
combined_summary <- dplyr::bind_rows(
  transcript_res$summary,
  activity_res$summary
)

print(combined_summary)

# Optionally save summary
write.csv(
  combined_summary,
  file = here::here("reproducibility", "figure6", "results", "Source_data_CTNNB1_paired_ttest_summary.csv"),
  row.names = FALSE
)
