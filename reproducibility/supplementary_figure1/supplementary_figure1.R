#!/usr/bin/env Rscript

library(cmapR)
library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

fig1_dir <- here("reproducibility", "figure1", "data")
results_dir <- here("reproducibility", "supplementary_figure1", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

file_names <- c(
  "level5_beta_trt_xpr_n142901x12328.gctx",
  "siginfo_beta.txt",
  "LINC gene annotations.csv"
)

missing_files <- file_names[!file.exists(file.path(fig1_dir, file_names))]
if (length(missing_files) > 0) {
  stop(
    "Missing Supplementary Figure 1 input files in ", fig1_dir, ": ",
    paste(missing_files, collapse = ", "),
    ". Download the Figshare data and place the Figure 1 inputs in this data folder."
  )
}

trt_xpr <- parse_gctx(file.path(fig1_dir, "level5_beta_trt_xpr_n142901x12328.gctx"))
meta_full <- read.delim(file.path(fig1_dir, "siginfo_beta.txt"))
meta_xpr <- meta_full %>% filter(pert_type == "trt_xpr")
annotated <- annotate_gct(trt_xpr, meta_xpr, dim = "column", keyfield = "sig_id")

gene_map <- read.csv(file.path(fig1_dir, "LINC gene annotations.csv"), stringsAsFactors = FALSE)
entrez2sym <- setNames(gene_map$gene.name, gene_map$entrez)
annotated@rdesc$symbol <- unname(entrez2sym[annotated@rid])

keep <- annotated@cdesc$tas > 0.2 &
  !is.na(annotated@cdesc$cmap_name) &
  annotated@cdesc$cmap_name != ""
a02 <- subset_gct(annotated, cid = which(keep))
a02@cdesc$id <- paste0(a02@cdesc$cmap_name, "_", a02@cdesc$sig_id)

df_mat <- as.data.frame(a02@mat)
rownames(df_mat) <- make.unique(a02@rdesc$symbol)
colnames(df_mat) <- a02@cdesc$id

akt1_ids <- a02@cdesc %>%
  as_tibble() %>%
  filter(cmap_name == "AKT1") %>%
  pull(id)
akt1_sigs <- df_mat[, colnames(df_mat) %in% akt1_ids, drop = FALSE]

lowest_is_akt1 <- apply(akt1_sigs, 2, function(column) {
  which.min(column) == which(rownames(akt1_sigs) == "AKT1")
})
akt1_sigs_top <- akt1_sigs[, lowest_is_akt1, drop = FALSE]
message("Selected signatures with AKT1 lowest, n=", ncol(akt1_sigs_top))

df_bot5 <- akt1_sigs_top %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "signature", values_to = "expression") %>%
  group_by(signature) %>%
  slice_min(expression, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_AKT1 = ifelse(gene == "AKT1", "AKT1", "Other"))

plot <- ggplot(df_bot5, aes(x = reorder(gene, expression), y = expression, fill = is_AKT1)) +
  geom_col() +
  facet_wrap(~signature, scales = "free_x") +
  scale_fill_manual(values = c(AKT1 = "red", Other = "grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(title = "Bottom 5 Genes per AKT1 Signature", x = "Gene", y = "Expression")

ggsave(
  file.path(results_dir, "Selected_AKT1_signatures.png"),
  plot,
  width = 8,
  height = 6,
  bg = "white"
)
