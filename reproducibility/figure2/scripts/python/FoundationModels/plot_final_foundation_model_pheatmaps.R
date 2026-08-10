suppressPackageStartupMessages({
  library(grid)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: Rscript plot_final_foundation_model_pheatmaps.R <gsva.csv> <bulkformer.csv> <aucell.csv> <uce4.csv> <uce33.csv> <geneformer_cell.csv> <geneformer_gene.csv> <output_dir>")
}

gsva <- read.csv(args[[1]], check.names = FALSE)
bulkformer <- read.csv(args[[2]], check.names = FALSE)
aucell <- read.csv(args[[3]], check.names = FALSE)
uce4 <- read.csv(args[[4]], check.names = FALSE)
uce33 <- read.csv(args[[5]], check.names = FALSE)
geneformer_cell <- read.csv(args[[6]], check.names = FALSE)
geneformer_gene <- read.csv(args[[7]], check.names = FALSE)
output_dir <- args[[8]]

signature_order <- c(
  "AKT1",
  "AKT1_XPR015_U251MG.311_96H:E14_c3",
  "AKT1_XPR008_A549.311_96H:G11_c3",
  "AKT1_XPR008_A375.311_96H:G11_c3"
)
context_labels <- c("Consensus", "Glioma", "NSCLC", "Melanoma")
signature_labels <- c("AKT1", "U251MG:E14", "A549:G11", "A375:G11")
colors <- colorRampPalette(c("#f7f7f7", "#d9e2ef", "#7f9bbd", "#1f4e79"))(100)

make_matrix <- function(dat, row_levels, value_column) {
  mat <- matrix(
    NA_real_,
    nrow = length(row_levels),
    ncol = length(signature_order),
    dimnames = list(row_levels, signature_order)
  )
  for (i in seq_len(nrow(dat))) {
    mat[as.character(dat$plot_method[i]), as.character(dat$signature[i])] <- dat[[value_column]][i]
  }
  mat
}

make_heatmap <- function(mat, cap, column_labels, title = NA_character_, legend = TRUE) {
  mat_plot <- pmin(mat, cap)
  labels <- matrix(sprintf("%.2f", mat), nrow = nrow(mat), dimnames = dimnames(mat))
  labels[is.na(mat)] <- ""
  pheatmap(
    mat_plot,
    color = colors,
    breaks = seq(0, cap, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = labels,
    number_color = "black",
    fontsize_number = 7,
    fontsize = 8,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    labels_col = column_labels,
    cellwidth = 38,
    cellheight = 18,
    border_color = "white",
    main = title,
    legend = legend,
    silent = TRUE
  )
}

save_single <- function(ph, path, width, height) {
  pdf(path, width = width, height = height, useDingbats = FALSE)
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
}

save_pair <- function(left, right, path, width, height) {
  pdf(path, width = width, height = height, useDingbats = FALSE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1, 1.08), "null"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(left$gtable)
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.draw(right$gtable)
  popViewport(2)
  dev.off()
}

gsva_plot <- data.frame(
  signature = gsva$signature,
  plot_method = "GSVA",
  value = -log10(gsva$p_value)
)
bulk_plot <- subset(
  bulkformer,
  model == "BulkFormer" & method %in% c("cosine", "gene_membership_spearman")
)
bulk_plot <- data.frame(
  signature = bulk_plot$signature,
  plot_method = ifelse(bulk_plot$method == "cosine", "Cosine", "Spearman"),
  value = bulk_plot$neglog10_p_two_sided
)
bulk_data <- rbind(gsva_plot, bulk_plot)
bulk_matrix <- make_matrix(bulk_data, c("GSVA", "Cosine", "Spearman"), "value")
panel_a <- make_heatmap(bulk_matrix, 2, context_labels)

aucell_plot <- data.frame(
  dataset = aucell$dataset,
  signature = aucell$signature,
  plot_method = "AUCell",
  value = -log10(aucell$mw_p / 2)
)
uce_plot <- rbind(
  data.frame(
    dataset = uce4$dataset,
    signature = uce4$signature,
    plot_method = "UCE-4",
    value = uce4$expected_neglog10_p
  ),
  data.frame(
    dataset = uce33$dataset,
    signature = uce33$signature,
    plot_method = "UCE-33",
    value = uce33$expected_neglog10_p
  )
)
uce_data <- rbind(aucell_plot, uce_plot)
uce_rows <- c("AUCell", "UCE-4", "UCE-33")
uce_a172_matrix <- make_matrix(subset(uce_data, dataset == "A172"), uce_rows, "value")
uce_u87_matrix <- make_matrix(subset(uce_data, dataset == "U87"), uce_rows, "value")
panel_b_a172 <- make_heatmap(uce_a172_matrix, 120, context_labels, "A172", FALSE)
panel_b_u87 <- make_heatmap(uce_u87_matrix, 120, context_labels, "U87", TRUE)

geneformer_cell_plot <- subset(geneformer_cell, model == "Geneformer cell embedding")
geneformer_cell_plot <- data.frame(
  dataset = geneformer_cell_plot$dataset,
  signature = geneformer_cell_plot$signature,
  plot_method = "Cell emb. cosine",
  value = geneformer_cell_plot$expected_neglog10_p
)
geneformer_gene_plot <- subset(
  geneformer_gene,
  representation == "pretrained_gene_embeddings" &
    method %in% c("cosine", "spearman")
)
geneformer_gene_plot <- data.frame(
  dataset = geneformer_gene_plot$dataset,
  signature = geneformer_gene_plot$signature,
  plot_method = ifelse(
    geneformer_gene_plot$method == "cosine",
    "Gene emb. cosine",
    "Gene emb. Spearman"
  ),
  value = geneformer_gene_plot$expected_neglog10_p
)
geneformer_data <- rbind(aucell_plot, geneformer_cell_plot, geneformer_gene_plot)
geneformer_rows <- c("AUCell", "Cell emb. cosine", "Gene emb. cosine", "Gene emb. Spearman")
geneformer_a172_matrix <- make_matrix(
  subset(geneformer_data, dataset == "A172"),
  geneformer_rows,
  "value"
)
geneformer_u87_matrix <- make_matrix(
  subset(geneformer_data, dataset == "U87"),
  geneformer_rows,
  "value"
)
panel_c_a172 <- make_heatmap(geneformer_a172_matrix, 120, signature_labels, "A172", FALSE)
panel_c_u87 <- make_heatmap(geneformer_u87_matrix, 120, signature_labels, "U87", TRUE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_single(panel_a, file.path(output_dir, "Figure_S2a_BulkFormer.pdf"), 4.2, 2.5)
save_pair(panel_b_a172, panel_b_u87, file.path(output_dir, "Figure_S2b_UCE.pdf"), 8.4, 2.7)
save_pair(panel_c_a172, panel_c_u87, file.path(output_dir, "Figure_S2c_Geneformer.pdf"), 8.6, 3.0)

pdf(
  file.path(output_dir, "Figure_S2_foundation_models.pdf"),
  width = 7.2,
  height = 7.3,
  useDingbats = FALSE
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  3,
  2,
  heights = unit(c(0.78, 1, 1.12), "null"),
  widths = unit(c(1, 1.08), "null")
)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(panel_a$gtable)
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(panel_b_a172$gtable)
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(panel_b_u87$gtable)
popViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
grid.draw(panel_c_a172$gtable)
popViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
grid.draw(panel_c_u87$gtable)
popViewport(2)
grid.text("a", x = unit(1.5, "mm"), y = unit(183, "mm"), just = c("left", "top"), gp = gpar(fontsize = 9, fontface = "bold"))
grid.text("b", x = unit(1.5, "mm"), y = unit(137, "mm"), just = c("left", "top"), gp = gpar(fontsize = 9, fontface = "bold"))
grid.text("c", x = unit(1.5, "mm"), y = unit(77, "mm"), just = c("left", "top"), gp = gpar(fontsize = 9, fontface = "bold"))
dev.off()
