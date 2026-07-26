#paste here esch validated step

#!/usr/bin/env Rscript
# Figure 2 Pipeline: RNA-Seq Reproducibility & GSEA across Cancer Contexts
# Author: Igor Cima
# Date:   2025-07-08

# ====== LIBRARIES ======
library(here)                     
library(data.table)               
library(dplyr)                    
library(tibble)                   
library(purrr)                    
library(DESeq2)                   
library(fgsea)                    
library(ggplot2)                  
library(scales)                  
library(Seurat)
library(AUCell)
library(GSEABase)
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(forcats)

# ====== SETUP ======
fig2_dir    <- here("reproducibility","figure2","data")
results_dir <- here("reproducibility","figure2","results")
dir.create(fig2_dir,    recursive=TRUE, showWarnings=FALSE)
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

# ====== STEP 1: Load and preprocess gene sets and annotations  ======

# AKT Gene set collection
pathways <- readRDS(file.path(fig2_dir,"Figure2a AKT gene set collection.rds"))

# Annotation for AKT gene set collection
annotation_gs <- read.csv(
  file.path(fig2_dir, "250617 Figure 2a_all gene sets_annotated_FINAL_github.csv"),
  stringsAsFactors = FALSE
)

# drop excluded sets
to_exclude <- annotation_gs$setName[annotation_gs$exclude == "y"]
pathways   <- pathways[! names(pathways) %in% to_exclude]

# Entrez->Symbol annotation
annot_file <- file.path(fig2_dir, "Human.GRCh38.p13.annot.tsv.gz")
if (!file.exists(annot_file)) {
  stop("Missing gene annotation file: ", annot_file)
}
annot <- fread(annot_file, quote="", data.table=FALSE)
annot$GeneID <- as.character(annot$GeneID)
rownames(annot) <- annot$GeneID

# ====== STEP 2: Helper Functions ======

# Load raw counts downloaded from Figshare
fetch_counts <- function(gse) {
  count_file <- file.path(fig2_dir, paste0(gse, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
  if (!file.exists(count_file)) {
    stop("Missing raw count file: ", count_file)
  }
  dt <- data.table::fread(count_file)
  mat <- as.matrix(dt[, -1, with=FALSE])
  rownames(mat) <- dt[[1]]
  mat
}

# Run DESeq2
run_deseq2 <- function(counts, gsms, annot, alpha=0.05) {
  sml_raw <- strsplit(gsms, "")[[1]]
  keep    <- sml_raw != "X"
  mat     <- counts[, keep]
  sml     <- factor(sml_raw[keep], levels=c("0","1"), labels=c("control","treatment"))
  coldata <- data.frame(Group=sml, row.names=colnames(mat))
  keep_genes <- rowSums(mat >= 1) >= min(table(sml))
  mat        <- mat[keep_genes, ]
  dds <- DESeqDataSetFromMatrix(countData=mat, colData=coldata, design=~Group)
  dds <- DESeq(dds, test="Wald", sfType="poscount")
  res <- results(dds, contrast=c("Group","treatment","control"), alpha=alpha)
  df  <- as.data.frame(res) %>%
    rownames_to_column("GeneID") %>%
    left_join(annot, by="GeneID") %>%
    dplyr::select(GeneID, Symbol, stat, log2FoldChange, pvalue, padj, baseMean)
  stats <- df$stat; names(stats) <- df$Symbol; sort(stats, decreasing=TRUE)
  list(results=df, stats=stats)
}

# Classify GSEA hits
classify_hits <- function(df) {
  df %>% mutate(
    classification = case_when(
      padj < 0.1 & NES >  1.5 & expected_NES == "pos" ~ "True Positive",
      padj < 0.1 & NES < -1.5 & expected_NES == "neg" ~ "True Positive",
      padj < 0.1 & NES >  1.5 & expected_NES == "neg" ~ "False Positive",
      padj < 0.1 & NES < -1.5 & expected_NES == "pos" ~ "False Positive",
      padj >= 0.1 | abs(NES) < 1.5                  ~ "False Negative",
      TRUE                                          ~ "Unclassified"
    )
  )
}

# Preprocess Seurat object: normalize, cell cycle regression, PCA+UMAP
env_preprocess_sc <- function(seurat_obj, s_genes, g2m_genes) {
  obj <- NormalizeData(seurat_obj)
  obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, vars.to.regress = c("S.Score","G2M.Score"))
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj)
  obj <- RunUMAP(obj, dims = 1:30)
  obj
}

# Build null AUCell model or load cache
env_build_null_model <- function(expr_mat,
                                 genes,
                                 n_sets     = 1000,
                                 set_size   = 200,
                                 cache_file = "AUCell_null_model.rds",
                                 seed = 1) {
  if (file.exists(cache_file)) {
    message("Loading cached null model from ", cache_file)
    return(readRDS(cache_file))
  }
  # build n_sets random gene lists (length set_size) excluding your AKT sets
  set.seed(seed)
  message("Computing and caching ", cache_file)
  null_sets <- replicate(n_sets,
                         sample(setdiff(genes, unique(unlist(pathways))),
                                set_size),
                         simplify = FALSE)
  names(null_sets) <- paste0("Random_", seq_along(null_sets))
  
  # compute rankings once
  rankings <- AUCell_buildRankings(expr_mat, splitByBlocks = TRUE)
  
  # for each random list, compute AUC and extract as a per‐cell vector
  null_aucs <- lapply(null_sets, function(gs) {
    auc_obj <- AUCell_calcAUC(gs, rankings, aucMaxRank = nrow(rankings) * 0.05)
    vec     <- as.numeric(getAUC(auc_obj))            # length = nCells
    names(vec) <- colnames(getAUC(auc_obj))           # cell barcodes
    vec
  })
  
  # bind into a matrix: rows = cells, cols = random sets
  null_matrix <- do.call(cbind, null_aucs)
  colnames(null_matrix) <- names(null_aucs)
  # save for caching
  saveRDS(null_matrix, cache_file)
  null_matrix
}

# Run AUCell and normalize by null model
env_run_aucell <- function(expr_mat, gene_sets, null_matrix = NULL, null_csv = NULL) {
  if (!is.null(null_matrix) && is.matrix(null_matrix)) {
    # use provided null_matrix
  } else if (!is.null(null_csv) && file.exists(null_csv)) {
    message("Loading null matrix from CSV: ", null_csv)
    null_matrix <- as.matrix(read.csv(null_csv, row.names = 1, check.names = FALSE))
  } else {
    stop("Must supply either a null_matrix object or a valid null_csv path.")
  }
  rankings <- AUCell_buildRankings(expr_mat, splitByBlocks = TRUE, plotStats = FALSE)
  cells_AUC <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = nrow(rankings) * 0.05)
  auc_raw <- getAUC(cells_AUC)
  raw_mat <- t(auc_raw)
  null_means <- rowMeans(null_matrix)
  null_sds   <- apply(null_matrix, 1, sd)
  norm1 <- sweep(raw_mat, 1, null_means)
  norm2 <- sweep(norm1, 1, null_sds, "/")
  list(raw = raw_mat, normalized = norm2)
}

# Plot UMAP by treatment (vehicle vs zstk474)
env_plot_umap <- function(seurat_obj, title_size = 8) {
  obj <- subset(seurat_obj, subset = dose != "1")
  DimPlot(
    obj,
    reduction = "umap",
    group.by  = "treatment",
    cols      = c("#e31a1c80", "#1f78b480")
  ) + theme(plot.title = element_text(size = title_size))
}

# Density plot
env_plot_density <- function(norm_mat, sig_name, meta, doses = c("0","10")) {
  df <- tibble(
    Expression = c(
      norm_mat[ meta$dose == doses[1], sig_name ],
      norm_mat[ meta$dose == doses[2], sig_name ]
    ),
    Dose = factor(
      rep(doses, times = c(
        sum(meta$dose == doses[1]),
        sum(meta$dose == doses[2])
      )),
      levels = doses,
      labels = paste0("Dose ", doses)
    )
  )
  ggplot(df, aes(x = Expression, fill = Dose)) +
    geom_density(alpha = 0.5, position = "identity") +
    scale_fill_manual(values = c("Dose 0" = "#e31a1c", "Dose 10" = "#1f78b4")) +
    labs(
      title = sig_name,
      x = "AUCell scores (norm.)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 5),
      panel.grid    = element_blank(),
      axis.line     = element_line(color = "black", size = 0.5),
      axis.ticks    = element_line(color = "black", size = 0.5),
      axis.text.x   = element_text(angle = 45, hjust = 1)
    )
}


# ====== STEP 3: Compute and Save ======

bubble_data <- list()
for (ctx in c("NSCLC","GBM","Melanoma")) {
  gse <- switch(ctx,
                NSCLC    = "GSE149626",
                GBM      = "GSE133116",
                Melanoma = "GSE111005")
  gsm <- switch(ctx,
                NSCLC    = "XX0011",
                GBM      = "00001111XXX",
                Melanoma = "1111XXXX0000")
  
  counts <- fetch_counts(gse)
  res    <- run_deseq2(counts, gsm, annot)
  
  set.seed(123)
  fg <- fgseaMultilevel(pathways, stats=res$stats, minSize=1, maxSize=500)
  fg_df <- as_tibble(fg) %>%
    inner_join(
      annotation_gs %>% dplyr::select(setName, expected_NES, context_gs),
      by = c("pathway" = "setName")
    ) %>%
    mutate(context_query = ctx)
  
  df0 <- classify_hits(fg_df) %>% filter(!is.na(NES), !is.na(padj))
  bubble_data[[ctx]] <- df0

##>>>> Fig 2a, Bubbleplots contexts <<<<  ----
    
  cont_cols <- c(
    "Human glioma"                      = "#619CFF",
    "Human melanoma"                     = "#F8766D",
    "Human lung cancer"                 = "#7CAE00",
    "Human prostate cancer"             = "#00B4F0",
    "Human breast cancer"               = "#00BFC4",
    "Consensus"                         = "#DE8C00",
    "Mouse prostate"                   = "#00C08B",
    "Mouse T cells"                    = "#00BA38",
    "Curated"                          = "#FF64B0",
    "Curated, normal T cells"          = "#C77CFF",
    "Curated, human prostate cancer"   = "#F564E3",
    "Human hepatocellular carcinoma"   = "#B79F00"
  )
  
    p <- ggplot(df0, aes(x=abs(NES), y=-log10(padj), color=context_gs, size=size)) +
    geom_point(data=df0 %>% filter(classification=="True Positive"), aes(fill=context_gs), shape=21) +
    geom_point(data=df0 %>% filter(classification=="False Positive"), shape=21, fill=NA) +
    geom_point(data=df0 %>% filter(classification=="False Negative"), aes(fill=context_gs), shape=21, alpha=0.2) +
    geom_vline(xintercept=1.5, linetype="dashed") +
    geom_hline(yintercept=-log10(0.1), linetype="dashed") +
    scale_x_continuous(limits=c(0,2.5), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,5),   expand=c(0,0)) +
    labs(title=paste0(ctx, " (", gse, ")"), x="|NES|", y="-log10(padj)") +
    scale_color_manual(values=cont_cols) + scale_fill_manual(values=cont_cols) +
    theme_minimal() + theme(
      legend.position="none", panel.grid.minor=element_blank(),
      axis.title=element_text(size=14), axis.text=element_text(size=14),
      axis.line=element_line(color="black",size=0.5), axis.ticks=element_line(color="black",size=0.5),
      axis.ticks.length=unit(0.25,"cm"), axis.text.x=element_text(angle=45,hjust=1)
    )+ scale_size_continuous(range = c(2, 12))
  
  ggsave(file.path(results_dir,paste0("Fig 2a_",ctx,"_bubble.pdf")), p, width=5, height=5, bg="white")
}

# Combine all contexts, drop leadingEdge, and save
all_bubble_df <- bind_rows(bubble_data, .id="context")
if("leadingEdge" %in% names(all_bubble_df)) {
  all_bubble_df <- dplyr::select(all_bubble_df, -leadingEdge)
}

library(cowplot)
legend <- cowplot::get_legend(
  p + 
    guides(color = guide_legend(override.aes = list(shape = 21, size = 5)),
           fill  = guide_legend(override.aes = list(shape = 21, size = 5))) +
    theme(legend.position = "right")
)
# then display or compose into a figure
legend <- cowplot::plot_grid(legend)

ggsave(
  filename = file.path(results_dir,paste0("Fig 2a_legend_bubble.pdf")),
  plot = legend,
  device = "pdf"
)



write.csv(
  all_bubble_df,
  file.path(results_dir, "Source_data_Fig 2a_bubble.csv"),
  row.names=FALSE
)

# End of Figure 2a pipeline

# ====== STEP 4: GSVA Regression for Figure 2b ======

# Load annotated ExpressionSet for GSE145128
#      (rows = unique genes, cols = samples; phenoData already in there)

siedSet <- readRDS(file.path(fig2_dir,
                             "Fig 2b_GSE145128 ExpressionSet_Annotated_UniqueGeneNames.rds"))

# Define the four AKT1 signatures + “AKT1” consensus
desired_sets <- c(
  "AKT1_XPR015_U251MG.311_96H:E14_c3",  # GBM U251MG:E14
  "AKT1_XPR008_A549.311_96H:G11_c3",  # NSCLC A549:G11
  "AKT1_XPR008_A375.311_96H:G11_c3",  # Melanoma A375:G11
  "AKT1"                             # consensus
)

#  Subset the Figure2a AKT gene‐set collection
#       (loaded already into `pathways` in STEP 1)

# AKT Gene set collection
pathways <- readRDS(file.path(fig2_dir,"Figure2a AKT gene set collection.rds"))


gs_list <- pathways[desired_sets]              # pick only those four

# 4b.4: Select the eight samples for regression and their “reference” values
selected_cols <- c(
  "GSM4307009","GSM4307012","GSM4307018","GSM4307014",
  "GSM4307010","GSM4307013","GSM4307019","GSM4307015"
)

reference <- c(11, 300, 60, 0, 450, 900, 1050, 220)

# 4b.5: Extract the expression matrix and run GSVA
expr_sel <- exprs(siedSet)[, selected_cols]

library(GSVA)
env_run_gsva <- function(expr_mat, gsets) {
  # wrap GSVA call in a helper to ensure reproducibility
  gsvapar <- gsvaParam(expr_mat, gsets, maxDiff = TRUE)
  gsva(gsvapar)
}
gsva_scores <- env_run_gsva(expr_sel, gs_list)

# 4b.6: Compute regression R² and one‐sided P for each gene‐set row
library(tibble)
library(purrr)
env_compute_regression <- function(scores, ref, cols) {
  sub <- scores[, cols, drop=FALSE]
  map_dfr(seq_len(nrow(sub)), function(i) {
    y <- as.numeric(sub[i, ])
    x <- ref
    m <- lm(y ~ x)
    tibble(
      setName = rownames(sub)[i],
      R2      = summary(m)$r.squared * sign(coef(m)[2]),
      P       = pt(coef(m)[2] / summary(m)$coefficients[2,2],
                   df.residual(m), lower.tail = FALSE)
    )
  })
}
reg_results <- env_compute_regression(gsva_scores, reference, selected_cols)

# 4b.7: FDR‐correct and save table
reg_results$P_adj <- p.adjust(reg_results$P, method = "fdr")
write.csv(
  reg_results,
  file.path(results_dir, "Source_data_Fig 2b_regression_results.csv"),
  row.names = FALSE
)

###>>>> Fig 2b, Scatterplots regression <<<<  ----

env_plot_regression <- function(gsva_row, ref, title) {
  df <- tibble(reference = ref, response = as.numeric(gsva_row))
  model <- lm(response ~ reference, df)
  r2    <- round(summary(model)$r.squared, 2)
  ggplot(df, aes(x = reference, y = response)) +
    geom_smooth(method="lm", color="black", fill="gray", alpha=0.6, size=0.5) +
    geom_point(color="darkred", size=2) +
    ggtitle(paste0(title, " — R²=", r2)) +
    labs(x="Reference Data", y="GSVA score") +
    theme_minimal() +
    theme(
      plot.title     = element_text(size=10),
      legend.position= "none",
      panel.grid     = element_blank(),
      axis.title     = element_text(size=14),
      axis.text      = element_text(size=14),
      axis.line      = element_line(color="black", size=0.5),
      axis.ticks     = element_line(color="black", size=0.5),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x    = element_text(angle=45, hjust=1)
    )
}

plots <- map(desired_sets, function(setName) {
  if (! setName %in% rownames(gsva_scores)) {
    warning("Missing GSVA row: ", setName)
    return(NULL)
  }
  env_plot_regression(gsva_scores[setName,], reference, setName)
}) %>% compact()

final_plot <- plot_grid(plotlist = plots, ncol = 2)
ggsave(
  file.path(results_dir, "Fig 2b_GSVA_regression.pdf"),
  final_plot, width=6, height=6, bg="white"
)

# ====== STEP 5: Single-cell AKT1‐AUCell (Fig 2d–f) ======

# Define the four AKT1 signatures + “AKT1” consensus
desired_sets <- c(
  "AKT1_XPR015_U251MG.311_96H:E14_c3",  # GBM U251MG:E14
  "AKT1_XPR008_A549.311_96H:G11_c3",  # NSCLC A549:G11
  "AKT1_XPR008_A375.311_96H:G11_c3",  # Melanoma A375:G11
  "AKT1"                             # consensus
)

#  Subset the Figure2a AKT gene‐set collection
#       (loaded already into `pathways` in STEP 1)

# AKT Gene set collection
pathways <- readRDS(file.path(fig2_dir,"Figure2a AKT gene set collection.rds"))


gs_list <- pathways[desired_sets]              # pick only those four


for (ctx in c("u87", "a172")) {
  
  # 5.1 Load processed Seurat object and pick U87 cells
  #     (change “u87” to “a172” if desired)
  combined_s <- readRDS(file.path(fig2_dir, "241105 combined_seurat_objects_subsetted and processed_vechicle&zstk474.rds"))
  s.obj      <- combined_s[[ctx]]
  
  # 5.2 Preprocess: normalize, regress out cell‐cycle, PCA, clusters & UMAP
  s.obj <- env_preprocess_sc(
    seurat_obj = s.obj,
    s_genes     = cc.genes$s.genes,
    g2m_genes   = cc.genes$g2m.genes
  )
  
  # 5.3 Build or load null AUCell model for random gene‐sets
  expr_mat    <- Seurat::GetAssayData(s.obj, assay="RNA", slot="data")
  null_file   <- file.path(fig2_dir, paste0(ctx, " AUCell_null_model.rds"))
  null_matrix <- env_build_null_model(
    expr_mat   = expr_mat,
    genes      = rownames(expr_mat),
    cache_file = null_file
  )
  
  # 5.4 Run AUCell for the four AKT1 gene‐sets from Fig 2a
  #     (desired_sets and pathways already defined in STEP 4b)
  auc_res <- env_run_aucell(
    expr_mat    = expr_mat,
    gene_sets   = gs_list,
    null_matrix = null_matrix
  )
  
#>>>> Fig 2c and e, UMAP cell lines <<<<  ----
  
  # 5.5 Plot UMAP colored by treatment (vehicle vs zstk474) → Fig 2d
  p_umap <- env_plot_umap(s.obj)
  ggsave(
    filename = file.path(results_dir, paste0("Fig 2c,e_UMAP_treatment_", ctx, ".pdf")),
    plot     = p_umap,
    width    = 6,
    height   = 5,
    bg       = "white"
  )

  #>>>> Fig 2c and e, Density plots glioma <<<<  ----
  
  # 5.6 Density plots for Dose 0 vs Dose 10 → Fig 2e
  density_plots <- purrr::map(
    desired_sets,
    function(sig) {
      if (!sig %in% colnames(auc_res$normalized)) {
        warning("Skipping density plot, signature not found: ", sig)
        return(NULL)
      }
      env_plot_density(
        norm_mat = auc_res$normalized,
        sig_name = sig,
        meta     = s.obj@meta.data
      )
    }
  ) %>% compact()
  
  # now combine the ones we actually got
  p_density <- cowplot::plot_grid(plotlist = density_plots, ncol = 2)
  
  ggsave(
    filename = file.path(results_dir, paste0("Fig 2c,e_Density_", ctx, ".pdf")),
    plot     = p_density,
    width    = 8,
    height   = 6,
    bg       = "white"
  )
  
  
  # 5.7 Compute Wilcoxon p‐values for Dose 10 vs Dose 0 → Fig 2f
  library(tibble)
  
  stats_df <- purrr::map_df(
    desired_sets,
    function(sig) {
      # skip any missing
      if (!sig %in% colnames(auc_res$normalized)) {
        warning("Signature not found, skipping: ", sig)
        return(tibble(Gene_Set = sig, p_wilcox = NA_real_))
      }
      vals <- auc_res$normalized[,sig]         # one row per signature
      md   <- s.obj@meta.data
      x0   <- vals[md$dose == 0]
      x10  <- vals[md$dose == 10]
      p    <- if (length(x0) > 1 && length(x10) > 1) {
        wilcox.test(x10, x0)$p.value
      } else {
        NA_real_
      }
      tibble(Gene_Set = sig, p_wilcox = p)
    }
  ) %>%
    dplyr::mutate(p_adj = p.adjust(p_wilcox, method = "fdr"))
  
  # save as source data
  write.csv(
    stats_df,
    file    = file.path(results_dir, paste0("Source_Fig 2d_Stats_", ctx, ".csv")),
    row.names = FALSE
  )
  
  #>>>> Fig 2d and f, Barplots contexts <<<<  ----
  #>
  # 5.8 Barplot of –log10(FDR) → Fig 2f
  p_bar <- ggplot(stats_df, aes(
    x = forcats::fct_reorder(Gene_Set, p_adj),
    y = -log10(p_adj)
  )) +
    geom_col(fill = "black") +
    labs(x = NULL, y = "-log10(FDR)") +
    theme_minimal() +
    theme(
      panel.grid    = element_blank(),
      axis.line     = element_line(color = "black", size = 0.5),
      axis.ticks    = element_line(color = "black", size = 0.5),
      axis.text.x   = element_text(angle = 45, hjust = 1),
      axis.title    = element_text(size = 14)
    )
  
  ggsave(
    filename = file.path(results_dir, paste0("Fig 2d,f_Barplot_", ctx, ".pdf")),
    plot     = p_bar,
    width    = 6,
    height   = 5,
    bg       = "white"
  )
  
}
