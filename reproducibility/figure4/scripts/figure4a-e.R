#!/usr/bin/env Rscript
# Script to reproduce Figure 4 analyses (Figure 4a–4d)

# Load required packages
library(echarts4r)           # Interactive charts (e.g., donut plots)
library(dplyr)               # Data manipulation (pipes, select, filter)
library(htmlwidgets)         # Save interactive widgets as HTML files
library(webshot2)            # Capture screenshots of HTML widgets
library(here)                # Simplify file paths (project-root)
library(readr)               # Fast CSV reading/writing
library(data.table)          # Fast data loading (e.g., fread)
library(DESeq2)              # Differential expression analysis
library(cmapR)               # LINCS data parsing (parse GCTX, LINC signatures)
library(GSEABase)            # Gene set data structures (GeneSetCollection)
library(fgsea)               # Fast gene set enrichment analysis (fgseaMultilevel)
library(stringr)             # String operations (regex, etc.)
library(ggplot2)             # Plotting system
library(patchwork)           # Combine multiple ggplots
library(precrec)             # Precision–Recall and ROC curve evaluation
library(tibble)              # Modern data frames (tribble, rownames_to_column)
library(purrr)               # Functional programming (map, map_dfr)

# Identify this script for the 'here' package
here::i_am("reproducibility/figure4/scripts/figure4a-e.R")

# Create output directories if they do not exist
data_dir    <- here("reproducibility", "figure4", "data")
results_dir <- here("reproducibility", "figure4", "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Define a utility function to generate COMPASS gene set collections
compass_gsc <- function(context,
                        subset_dir   = here("reproducibility", "figure3", "data", "subsets"),
                        n            = 200,
                        min_conf     = 1,
                        targets      = NULL,
                        driver_filter= FALSE,
                        output       = c("list", "df", "gsc")) {
  output <- match.arg(output)
  require(cmapR)
  require(GSEABase)
  
  # Load the precomputed subset GCT for the given context
  subset_file <- file.path(subset_dir, paste0(context, "_subset.rds"))
  if (!file.exists(subset_file)) {
    stop("Missing subset for context: ", context)
  }
  gct <- readRDS(subset_file)
  
  # Filter signatures by confidence and optionally by driver status
  idx <- which(gct@cdesc$cps_conf_total >= min_conf)
  if (driver_filter) {
    idx <- idx[gct@cdesc$cancer_driver_summary[idx] != "None"]
  }
  if (length(idx) == 0) {
    if (output == "gsc") return(GeneSetCollection(list()))
    if (output == "df")  return(data.frame())
    return(list())
  }
  
  # Build matrix of gene-level stats for selected signatures
  mat <- gct@mat[ , idx, drop = FALSE]
  rownames(mat) <- gct@rdesc$symbol
  ids  <- gct@cdesc$id[idx]
  conf <- gct@cdesc$cps_conf_total[idx]
  names_out <- paste0(ids, "_c", conf)
  
  # Optional target filter: keep only signatures targeting specified genes
  if (!is.null(targets)) {
    keep_tgt <- gct@cdesc$cmap_name[idx] %in% targets
    if (!any(keep_tgt)) {
      warning("No signatures for targets: ", paste(targets, collapse = ", "))
      if (output == "gsc") return(GeneSetCollection(list()))
      if (output == "df")  return(data.frame())
      return(list())
    }
    mat <- mat[ , keep_tgt, drop = FALSE]
    names_out <- names_out[keep_tgt]
  }
  
  # Extract bottom n genes (most down-regulated) for each signature
  res_list <- lapply(seq_len(ncol(mat)), function(j) {
    v <- mat[ , j]
    o <- order(v, decreasing = FALSE, na.last = TRUE)
    head(rownames(mat)[o], n)
  })
  names(res_list) <- names_out
  
  # Return result in requested format
  if (output == "gsc") {
    gs_list <- lapply(names_out, function(nm) GeneSet(geneIds = res_list[[nm]], setName = nm))
    return(GeneSetCollection(gs_list))
  } else if (output == "df") {
    maxlen <- max(lengths(res_list))
    df <- do.call(cbind, lapply(res_list, function(x) { length(x) <- maxlen; x }))
    return(as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE))
  } else {
    return(res_list)
  }
}


# Figure 4a: Donut plots for GBM63 dataset ----------------------------------------

# Step 1: Setup paths to input data
df1_path <- file.path(data_dir, "annotation_protein localization_GBM63.csv")
df2_path <- file.path(data_dir, "annotation_protein function_GBM63.csv")

# Step 2: Read files and prepare data
df1 <- read.csv(df1_path)
df2 <- read.csv(df2_path)
# Select relevant columns and rename for plotting
plot_data1 <- df1 %>%
  dplyr::select(Subcellular.localization, Glioblastoma) %>%
  dplyr::rename(Compass = Glioblastoma)
plot_data2 <- df2 %>%
  dplyr::select(Function, Glioblastoma63) %>%
  dplyr::rename(Compass = Glioblastoma63)

# Define custom color palettes
custom_colors  <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                    "#e6ab02", "#a6761d", "#666666", "#d53e4f", "#3288bd")  # for localization
custom_colors2 <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", 
                    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")  # for function

# Step 3: Create donut plots (function and localization)
# Top: Function donut (Figure 4a top)
p1 <- plot_data2 %>%
  e_charts(Function, renderer = "svg") %>%
  e_pie(Compass, radius = c("30%", "50%")) %>%   # Donut: inner radius 30%, outer radius 50%
  e_labels(show = TRUE, position = "outside", formatter = "{b}", fontSize = 20) %>%   # Label with function names
  e_legend(show = FALSE) %>%                     # Hide legend
  e_color(custom_colors2) %>%                    # Apply custom colors
  e_tooltip(trigger = "item") %>%
  e_show_loading()

# Bottom: Localization donut (Figure 4a bottom)
p2 <- plot_data1 %>%
  e_charts(Subcellular.localization, renderer = "svg") %>%
  e_pie(Compass, radius = c("30%", "50%")) %>%   # Donut: inner radius 30%, outer radius 50%
  e_labels(show = TRUE, position = "outside", formatter = "{b}", fontSize = 20) %>%   # Label with localization names
  e_legend(show = FALSE) %>%                     # Hide legend
  e_color(custom_colors) %>%                     # Apply custom colors
  e_tooltip(trigger = "item") %>%
  e_show_loading()

# Step 4: Save donut plots as standalone HTML files
html_file1 <- file.path(results_dir, "Fig4a_donut_function.html")
html_file2 <- file.path(results_dir, "Fig4a_donut_localization.html")
saveWidget(p1, file = html_file1, selfcontained = TRUE)
saveWidget(p2, file = html_file2, selfcontained = TRUE)

# Figure 4b: Random walk enrichment plots ---------------------------------------

# Define target genes and corresponding GEO series with sample codes
params <- tibble::tribble(
  ~GENE, ~GSE,         ~gsms,
  "ASNS", "GSE171163", "XXXXXXXXXXXXXXXX000111XXXXXX",
  "RHOA", "GSE111571", "XXX000XXXXXXXXXXXX11",
  "CD109","GSE169418", "XXXXXXXXXXXX0101XXXX",
  "MYC",  "GSE86518",  "0X1XXXX0XXXX10X1X0X10X1XXXXXX"
)
# Path to COMPASS subset data (for glioma context)
subset_dir <- here("reproducibility", "figure3", "data", "subsets")

# Analysis loop for each target gene
results <- list()  # store NES, p-value, and plot for each gene

for (i in seq_len(nrow(params))) {
  # Extract parameters for this gene
  GENE <- params$GENE[i]
  GSE  <- params$GSE[i]
  gsms <- params$gsms[i]
  
  # Load raw counts downloaded from Figshare
  count_file <- file.path(data_dir, paste0(GSE, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
  if (!file.exists(count_file)) {
    stop("Missing raw count file: ", count_file)
  }
  counts <- data.table::fread(count_file)
  mat <- as.matrix(counts[ , -1, with = FALSE])
  rownames(mat) <- counts[[1]]
  
  # Load gene annotation (GRCh38) for mapping GeneIDs to symbols
  annot_file <- file.path(data_dir, "Human.GRCh38.p13.annot.tsv.gz")
  if (!file.exists(annot_file)) {
    stop("Missing gene annotation file: ", annot_file)
  }
  annot_df <- data.table::fread(annot_file, data.table = FALSE, stringsAsFactors = FALSE)
  annot_df$GeneID <- as.character(annot_df$GeneID)
  rownames(annot_df) <- annot_df$GeneID
  
  # Filter samples by 'X' in gsms string (X = exclude sample)
  sml <- strsplit(gsms, "")[[1]]
  keep_idx <- which(sml != "X")
  mat <- mat[ , keep_idx]
  sml <- sml[keep_idx]
  # Assign group labels: 'ctrl' vs 'treat'
  groups <- factor(ifelse(sml == "1", "treat", "ctrl"), levels = c("ctrl", "treat"))
  coldata <- data.frame(Group = groups, row.names = colnames(mat))
  
  # Pre-filter low-count genes (require >=1 count in at least min group size samples)
  min_samps <- min(table(groups))
  keep_genes <- rowSums(mat >= 1) >= min_samps
  mat <- mat[keep_genes, ]
  
  # Differential expression analysis using DESeq2
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~Group)
  dds <- DESeq(dds, test = "Wald", sfType = "poscount")
  res <- results(dds, contrast = c("Group", "treat", "ctrl"), alpha = 0.05)
  
  # Prepare named vector of test statistics for fgsea
  stat_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("GeneID") %>%
    filter(!is.na(stat)) %>%
    left_join(annot_df, by = "GeneID")
  stat_df_unique <- stat_df %>%
    group_by(Symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  stat_vec <- setNames(stat_df_unique$stat, stat_df_unique$Symbol)
  stat_vec <- sort(stat_vec, decreasing = TRUE)
  
  # Load context-specific COMPASS gene sets (top 250 genes per set for glioma)
  gs_list <- compass_gsc(context = "glioma", subset_dir = subset_dir, n = 250, output = "list")
  # Filter to gene sets matching the target gene (names start with e.g. "ASNS_")
  gs_list <- gs_list[ grepl(paste0("^", GENE, "_"), names(gs_list)) ]
  
  # Run FGSEA on the statistics vector
  set.seed(123)
  fg <- fgseaMultilevel(pathways = gs_list, stats = stat_vec, minSize = 1, maxSize = 500)
  
  # 2.10) get all matching gene‐sets for this GENE
  gs_names <- grep(paste0("^", GENE, "_"), fg$pathway, value = TRUE)
  
  # 2.11) extract their confidence N from the "_cN" suffix
  confs <- as.integer(sub(".*_c(\\d+)$", "\\1", gs_names))
  
  # 2.12) find the max confidence and sample one if there’s a tie
  max_conf     <- max(confs, na.rm = TRUE)
  best_indices <- which(confs == max_conf)
  sel_name     <- if (length(best_indices) == 1) {
    gs_names[best_indices]
  } else {
    sample(gs_names[best_indices], 1)
  }
  
  # 2.13) now just do one enrichment plot for sel_name
  nes  <- fg$NES[ fg$pathway == sel_name ]
  pval <- fg$pval[ fg$pathway == sel_name ]
  
  p <- plotEnrichment(
    pathway = gs_list[[sel_name]],
    stats   = stat_vec
  ) +
    labs(
      title    = paste0(GENE, " | ", sel_name, " (NES=", round(nes,2), ")"),
      subtitle = paste0("p=", signif(pval,3))
    ) +
    theme_minimal() +
    theme(
      axis.text.y       = element_text(size = 12),
      axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
      axis.line         = element_line(color = "black", size = 0.5),
      axis.ticks        = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid        = element_blank(),
      legend.position   = "bottom"
    )
  
  # 2.14) Store only that one plot, NES & pval
  results[[GENE]] <- list(
    plot = p,
    NES  = nes,
    pval = pval
  )
}

# extract the single plot per gene
all_plots <- purrr::map(results, "plot")

# arrange them in a 2×2 grid
combined_plots <- patchwork::wrap_plots(all_plots, ncol = 2)

# save
ggplot2::ggsave(
  filename = file.path(results_dir, "Fig4b_random_walks_allSets.pdf"),
  plot     = combined_plots,
  width    = 10, 
  height   = 8,
  bg       = "white"
)

# Figure 4c and d: Precision-Recall (PR) curves for GBM63 ----------------------

# Load benchmarking z-scores for COMPASS (Glioblastoma 63 dataset)
input_file <- file.path(data_dir, "250212 BENCHMARKING_GBM63_DESeq2MEAN_COMPASS&SigCOM_FILTERED_FOR_ROC_zscores.csv")
data3 <- read_csv(input_file, col_types = cols())

# Flip signs for z-score columns 9–16 (for PR curves, higher = better score)
data3[, 9:16] <- -data3[, 9:16]

# PR curve: Nuclear access vs No access
# Subset data by nuclear localization flag
data_nuc   <- data3 %>% filter(Protein_Access_to_Nucleus == "Yes")
data_nonuc <- data3 %>% filter(Protein_Access_to_Nucleus == "No")
# Compute PRC for each subset using Z_250 scores
pr_nuc <- evalmod(
  scores = data_nuc$Z_250, labels = data_nuc$Label,
  modnames = "Nuclear", mode = "rocprc"
)
pr_nonuc <- evalmod(
  scores = data_nonuc$Z_250, labels = data_nonuc$Label,
  modnames = "NoNuclear", mode = "rocprc"
)
# Extract and combine PR curve data for plotting
df_nuc   <- fortify(pr_nuc)   %>% mutate(Source = "Nuclear access")
df_nonuc <- fortify(pr_nonuc) %>% mutate(Source = "No nuclear access")
combined_nuc_df <- bind_rows(df_nuc, df_nonuc)
# Calculate AUC values and sample counts
auc_nuc   <- round(attr(pr_nuc,   "auc")[2,4], 3)
auc_nonuc <- round(attr(pr_nonuc, "auc")[2,4], 3)
count_nuc   <- attr(pr_nuc,   "data_info")[4]
count_nonuc <- attr(pr_nonuc, "data_info")[4]
# Plot PR curves for nuclear vs no-nuclear access groups
p_nuc <- ggplot(combined_nuc_df %>% filter(curvetype == "PRC"),
                aes(x = x, y = y, color = Source)) +
  geom_line(size = 1) +
  labs(title = paste0("Nuclear access: ", auc_nuc, " vs No access: ", auc_nonuc),
       x = "Recall", y = "Precision",
       color = paste0("n=", count_nuc, ",", count_nonuc)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title     = element_text(size = 10),
        legend.key.size  = unit(0.5, "cm"),
        legend.text      = element_text(size = 12),
        panel.border     = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title       = element_text(size = 14),
        axis.text        = element_text(size = 14),
        axis.ticks       = element_line(color = "black", size = 0.5),
        axis.ticks.length= unit(0.25, "cm")) +
  scale_color_manual(values = c("Nuclear access" = "#0072B2",
                                "No nuclear access" = "#D55E00"))

# PR curves: by protein primary function group
# Re-categorize certain functions for grouping
data_sign  <- data3 %>% filter(Protein_Primary_Function == "Signaling")
data_reg   <- data3 %>% filter(Protein_Primary_Function == "Regulatory")
data_enzf  <- data3 %>% 
  filter(Protein_Primary_Function %in% c("Metabolic Enzyme", "Protein Folding")) %>%
  mutate(Protein_Primary_Function = "Enzyme&Folding")
data_extr  <- data3 %>% 
  filter(Protein_Primary_Function %in% c("Extracellular Effectors", "Extracellular Signaling", "Adhesion")) %>%
  mutate(Protein_Primary_Function = "Extracellular")
data_trans <- data3 %>% filter(Protein_Primary_Function == "Transport")
# Compute PRC for each function category
pr1 <- evalmod(
  scores = data_sign$Z_250, labels = data_sign$Label,
  modnames = "Signaling", mode = "rocprc"
)
pr2 <- evalmod(
  scores = data_reg$Z_250, labels = data_reg$Label,
  modnames = "Regulatory", mode = "rocprc"
)
pr3 <- evalmod(
  scores = data_enzf$Z_250, labels = data_enzf$Label,
  modnames = "Enzyme&Folding", mode = "rocprc"
)
pr4 <- evalmod(
  scores = data_extr$Z_250, labels = data_extr$Label,
  modnames = "Extracellular", mode = "rocprc"
)
pr5 <- evalmod(
  scores = data_trans$Z_250, labels = data_trans$Label,
  modnames = "Transport", mode = "rocprc"
)
# Combine PR curve data from all groups
df1 <- fortify(pr1) %>% mutate(Source = "Signaling")
df2 <- fortify(pr2) %>% mutate(Source = "Regulatory")
df3 <- fortify(pr3) %>% mutate(Source = "Enzyme&Folding")
df4 <- fortify(pr4) %>% mutate(Source = "Extracellular")
df5 <- fortify(pr5) %>% mutate(Source = "Transport")
combined_func_df <- bind_rows(df1, df2, df3, df4, df5)
# Calculate AUC for each functional group
auc_vals <- c(
  Signaling     = round(attr(pr1, "auc")[2,4], 3),
  Regulatory    = round(attr(pr2, "auc")[2,4], 3),
  EnzymeFolding = round(attr(pr3, "auc")[2,4], 3),
  Extracellular = round(attr(pr4, "auc")[2,4], 3),
  Transport     = round(attr(pr5, "auc")[2,4], 3)
)
# Plot PR curves for each function group
p_func <- ggplot(combined_func_df %>% filter(curvetype == "PRC"),
                 aes(x = x, y = y, color = Source)) +
  geom_line() +
  labs(title = "PRC by Protein Function",
       subtitle = paste0("AUC: ", paste(auc_vals, collapse = ", ")),
       x = "Recall", y = "Precision") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title     = element_text(size = 10),
        legend.key.size  = unit(0.5, "cm"),
        legend.text      = element_text(size = 12),
        panel.border     = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title       = element_text(size = 14),
        axis.text        = element_text(size = 14),
        axis.ticks       = element_line(color = "black", size = 0.5),
        axis.ticks.length= unit(0.25, "cm")) +
  scale_color_manual(values = c("Signaling"      = "#0072B2",
                                "Regulatory"     = "#D55E00",
                                "Enzyme&Folding" = "#009E73",
                                "Extracellular"  = "#D62728",
                                "Transport"      = "#CC79A7"))

# Save PR curve plots to files
ggsave(filename = file.path(results_dir, "Fig4c_PRC_function.pdf"),
       plot = p_func, width = 6, height = 4, bg = "white")
ggsave(filename = file.path(results_dir, "Fig4d_PRC_nuclear_access.pdf"),
       plot = p_nuc, width = 6, height = 4, bg = "white")

# Figure 4e: Drug-target interaction NES dot plot ------------------------------

# Load and clean drug signature data
gct_file   <- file.path(data_dir, "drugs_U251MG_withLINCsignatures.gct")
meta_file  <- file.path(data_dir, "metadata drugs_U251MG_withLINCsignatures.txt")
annot_file <- file.path(data_dir, "LINC gene annotations.csv")
drug_file  <- file.path(data_dir, "COMPASS_drug_target interactions.csv")
# Parse the GCT file containing drug perturbation signatures
gct <- parse_gctx(gct_file)
# Annotate rows with gene symbols and inference data
gene_map <- read_csv(annot_file, col_types = cols())
entrez2sym <- setNames(gene_map$'gene name', gene_map$entrez)
infer_map  <- setNames(gene_map$inference, gene_map$'gene name')
gct@rdesc$symbol    <- unname(entrez2sym[gct@rid])
gct@rdesc$inference <- unname(infer_map[gct@rdesc$symbol])
gct@cdesc$sig_id    <- gct@cid
# Merge additional metadata (dose, etc.) into column descriptors
meta <- read.delim(meta_file, stringsAsFactors = FALSE)
gct@cdesc <- left_join(as.data.frame(gct@cdesc), meta, by = "sig_id")
# Load drug–target interactions and incorporate into metadata
drug_df <- read_csv(drug_file, col_types = cols(
  Drug = col_character(), Target = col_character(), Mode_of_action = col_character()
))
drug_meta <- drug_df %>%
  group_by(Drug) %>%
  summarise(Targets = list(unique(Target)), MoA = list(unique(Mode_of_action)), .groups = "drop")
gct@cdesc <- gct@cdesc %>%
  left_join(drug_meta, by = c("cmap_name" = "Drug")) %>%
  mutate(Targets = ifelse(is.na(Targets), list(character()), Targets),
         MoA     = ifelse(is.na(MoA),     list(character()), MoA))

# Optionally save the cleaned GCT object for reuse
saveRDS(gct, file = file.path(data_dir, "drugs_U251MG_withLINCsignatures.cleaned.rds"))

# Perform gene set enrichment analysis (FGSEA) on each drug signature
mat <- as.data.frame(gct@mat)
rownames(mat) <- make.unique(gct@rdesc$symbol)
# Compile all unique targets from the drug metadata
all_targets <- unique(unlist(drug_meta$Targets))
# Load COMPASS gene sets for glioma context (top 250 genes), filtered to relevant targets
comp_sets <- compass_gsc("glioma", n = 250, targets = all_targets, output = "list")
# Run FGSEA for each drug signature (each column in the matrix)
set.seed(123)
fgsea_list <- lapply(colnames(mat), function(sig) {
  stats <- mat[ , sig]
  names(stats) <- rownames(mat)
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  res <- fgseaMultilevel(pathways = comp_sets, stats = stats,
                         minSize = 10, maxSize = 500, nPermSimple = 5000)
  res$sig_id <- sig
  res
})
fgsea_res <- bind_rows(fgsea_list)

# Annotate FGSEA results with drug and target info
sig_info <- gct@cdesc %>%
  dplyr::select(sig_id, cmap_name, pert_idose, Targets, MoA) %>%
  dplyr::mutate(Targets = purrr::map_chr(Targets, ~ paste(., collapse = ",")),
         MoA     = purrr::map_chr(MoA,     ~ paste(., collapse = ",")))

gf <- fgsea_res %>% left_join(sig_info, by = "sig_id")

# Select specific drugs and targets of interest for plotting
sel_drugs   <- c("palbociclib", "thioridazine", "afatinib", "simvastatin",
                 "equilin", "AZD-8055", "olaparib", "BI-2536", "daunorubicin")
sel_targets <- c("CDK4", "CDK6", "DRD1", "DRD2", "EGFR", "ERBB2", "ERBB4",
                 "MTOR", "HMGCR", "HSD17B1", "PARP1", "PARP2", "PLK1", "TOP2A")

# Filter results for selected drugs (at specified doses) and targets
gf_filtered <- gf %>%
  mutate(pert_dose = str_replace_all(pert_idose, "µ", "u") %>% str_squish()) %>%
  filter(cmap_name %in% sel_drugs,
         pert_dose == c(palbociclib = "10 uM", thioridazine = "10 uM",
                        afatinib   = "10 uM", simvastatin  = "10 uM",
                        equilin    = "10 uM", `AZD-8055`   = "0.12 uM",
                        olaparib   = "10 uM", `BI-2536`    = "10 uM",
                        daunorubicin = "10 uM")[cmap_name]) %>%
  mutate(Target = sub("^(.*?)_.*$", "\\1", pathway)) %>%
  filter(Target %in% sel_targets)
# For each drug–target pair, choose the gene set with highest confidence (conf)
plot_df <- gf_filtered %>%
  mutate(conf = as.integer(str_extract(pathway, "(?<=_c)\\d+$"))) %>%
  group_by(cmap_name, Target) %>%
  filter((Target == "HMGCR" & conf == 2) | (Target != "HMGCR" & conf == max(conf))) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  dplyr::select(-conf)

write.csv(plot_df[, -8], paste0(results_dir, "/Source Data_Fig4e_NES_dotplot.csv"))
write.csv(gf[, -8], paste0(results_dir, "/Supplementary Table 4_all inhibitors tested.csv"))


# Compute plot aesthetics (point size, color, transparency based on NES and significance)
plot_data <- plot_df %>%
  mutate(NES_size = abs(NES),
         color_group = case_when(
           padj < 0.05 & NES < 0 ~ "#1f78b4",   # significant & negative NES (blue)
           padj < 0.05 & NES > 0 ~ "#e31a1c",   # significant & positive NES (red)
           TRUE                  ~ "grey"      # not significant
         ),
         log_padj = -log10(padj),
         transparency = ifelse(padj < 0.055,
                               scales::rescale(log_padj, to = c(0.2, 1)),
                               0.4)) %>%
  mutate(cmap_name = factor(cmap_name, levels = rev(sel_drugs)))

# Create dot plot of NES values for selected drug–target pairs
p_dot <- ggplot(plot_data, aes(x = cmap_name, y = Target,
                               size = NES_size, color = color_group, alpha = transparency)) +
  geom_point() +
  coord_flip() +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line   = element_line(color = "black"),
        axis.ticks  = element_line(color = "black")) +
  labs(x = "Drug", y = "Protein",
       title = "Fig 4e | NES Dot Plot",
       color = "Significance", size = "|NES|",
       alpha = "Transparency\n(-log10 p_adj)")

# Save dot plot to file
ggsave(file.path(results_dir, "Fig4e_NES_dotplot.pdf"),
       plot = p_dot, width = 6, height = 4, bg = "white")
