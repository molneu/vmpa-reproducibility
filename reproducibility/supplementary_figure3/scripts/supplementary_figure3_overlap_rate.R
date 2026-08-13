#!/usr/bin/env Rscript

# ----------------------------
# Overlap rate boxplots (compact, parallel)
# ----------------------------

# Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require(GSEABase)) BiocManager::install("GSEABase")
if (!require(GSVA))     BiocManager::install("GSVA")
if (!require(dplyr))    install.packages("dplyr")
if (!require(ggplot2))  install.packages("ggplot2")
if (!require(here))     install.packages("here")
if (!require(future))   install.packages("future")
if (!require(furrr))    install.packages("furrr")

library(GSEABase)
library(GSVA)
library(dplyr)
library(ggplot2)
library(here)
library(future)
library(furrr)

here::i_am("reproducibility/supplementary_figure3/scripts/supplementary_figure3_overlap_rate.R")

# compass_gsc() function
source(here("reproducibility", "figure3", "scripts", "preprocessing", "5_compass_gsc_function.R"))

# Directories
data_dir    <- here("reproducibility", "figure3", "data", "subsets")
results_dir <- here("reproducibility", "supplementary_figure3", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
contexts <- c("glioma","melanoma","nsclc","ovarian", "crc","breast","prostate","pdac")
sizes    <- c(25, 50, 100, 150, 200, 250, 300)

# Helper: Tukey-style boxplot stats
compute_box_stats <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  iqr <- q[3] - q[1]
  lower_whisker <- max(min(x, na.rm = TRUE), q[1] - 1.5 * iqr)
  upper_whisker <- min(max(x, na.rm = TRUE), q[3] + 1.5 * iqr)
  tibble(
    ymin   = lower_whisker,
    lower  = q[1],
    middle = q[2],
    upper  = q[3],
    ymax   = upper_whisker
  )
}

# Parallel plan
plan(cluster, workers = parallel::detectCores() - 1)
options(future.globals.maxSize = 5 * 1024^3)

# Expand grid of all tasks
params <- expand.grid(Context = contexts, Size = sizes, stringsAsFactors = FALSE)

# Parallel computation
boxplot_stats <- future_pmap(params, function(Context, Size) {
  message("Processing ", Context, " size ", Size)
  
  # Build GeneSetCollection with your function
  gsc <- compass_gsc(context = Context,
                     subset_dir = data_dir,
                     n = Size,
                     min_conf = 1,
                     driver_filter = FALSE,
                     output = "gsc")
  
  if (length(gsc) == 0) return(NULL)
  
  # Compute overlaps
  all_genes <- unique(unlist(lapply(gsc, geneIds)))
  overlap <- computeGeneSetsOverlap(gsc, all_genes)
  diag(overlap) <- NA
  
  # Exclude intra-prefix overlaps
  set_prefix <- sub(":.*", "", rownames(overlap))
  keep <- outer(set_prefix, set_prefix, FUN = "!=")
  overlap <- overlap * keep
  
  overlaps <- overlap[upper.tri(overlap)]
  overlaps <- overlaps[!is.na(overlaps)]
  
  if (length(overlaps) == 0) return(NULL)
  
  stats <- compute_box_stats(overlaps) %>%
    mutate(Context = Context, Size = Size)
  stats
})

# Combine
boxplot_df <- bind_rows(boxplot_stats)

# Plot
p <- ggplot(boxplot_df, aes(x = factor(Size), fill = Context)) +
  geom_boxplot(
    aes(ymin = ymin, lower = lower, middle = middle,
        upper = upper, ymax = ymax),
    stat = "identity"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Gene set size", y = "Overlap rate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p
# Save
ggsave(file.path(results_dir, "Supplementary_Figure_3_overlap_gene_sets_boxplots.pdf"),
       p, width = 10, height = 6, bg = "white")

# Save source data for Supplementary Figure 3
write.csv(
  boxplot_df,
  file = file.path(results_dir, "Source_Data_Supplementary_Figure_3_overlap_gene_sets_boxplots.csv"),
  row.names = FALSE
)
