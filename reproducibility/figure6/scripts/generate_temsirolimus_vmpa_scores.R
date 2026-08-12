#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(vmpaR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript generate_temsirolimus_vmpa_scores.R <annotated-raw-counts.csv> <output-directory>",
    call. = FALSE
  )
}

raw_file <- normalizePath(args[1], mustWork = TRUE)
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(raw_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("ensembl_gene_id", "gene_symbol") %in% names(raw))) {
  stop(
    "Raw-count input must contain ensembl_gene_id and gene_symbol columns.",
    call. = FALSE
  )
}
if (anyDuplicated(raw$ensembl_gene_id)) {
  stop("Raw-count input contains duplicated Ensembl identifiers.", call. = FALSE)
}

sample_ids <- setdiff(names(raw), c("ensembl_gene_id", "gene_symbol"))
counts <- as.matrix(raw[, sample_ids, drop = FALSE])
storage.mode(counts) <- "numeric"
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
  stop("Raw counts must be finite, non-negative and non-missing.", call. = FALSE)
}

symbols <- raw$gene_symbol
keep <- !is.na(symbols) & nzchar(symbols)
counts <- counts[keep, , drop = FALSE]
symbols <- symbols[keep]

counts_sum <- rowsum(counts, group = symbols, reorder = TRUE)
counts_sum <- counts_sum[rowSums(counts_sum) > 0, , drop = FALSE]

dge <- edgeR::DGEList(counts = counts_sum)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
expression_tmm_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

write.csv(
  data.frame(hgnc_symbol = rownames(counts_sum), counts_sum, check.names = FALSE),
  file.path(output_dir, "temsirolimus_gene_symbol_counts_collapsed_sum.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(hgnc_symbol = rownames(expression_tmm_logcpm), expression_tmm_logcpm, check.names = FALSE),
  file.path(output_dir, "temsirolimus_gene_symbol_TMM_logCPM.csv"),
  row.names = FALSE
)

for (scaling in c("none", "sample_z", "sample_pop_sd", "signature_z")) {
  result <- vmpaR::vmpa(
    input = expression_tmm_logcpm,
    context = "glioma",
    algorithm = "gsva",
    gsva_score_scaling = scaling,
    unique = FALSE,
    n = 250L,
    min_conf = 1L,
    driver_filter = FALSE,
    verbose = TRUE,
    gsva_min_size = 1L
  )

  write.csv(
    as.data.frame(result, check.names = FALSE),
    file.path(
      output_dir,
      paste0("vmpa_glioma_n250_", scaling, "_unique_FALSE_raw_result.csv")
    ),
    row.names = FALSE
  )
}
