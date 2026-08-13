#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(precrec)
  library(ggplot2)
  library(here)
})

here::i_am("reproducibility/figure3/scripts/figure3d.R")

data_dir <- here("reproducibility", "figure3", "data")
results_dir <- here("reproducibility", "figure3", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

observed_file <- file.path(
  data_dir,
  "figure3d_observed_common_intersection.csv"
)
viper_null_file <- file.path(
  data_dir,
  "figure3d_viper_100_null_summary.csv"
)
fgsea_null_file <- file.path(
  data_dir,
  "250715_fGSEA_100NULLmodels_SUMMARY.csv"
)
size_file <- file.path(
  data_dir,
  "250710 Benchmarking_fgsea results_ALL_CONTEXTS.csv"
)

input_files <- c(observed_file, viper_null_file, fgsea_null_file, size_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing Figure 3d input files:\n", paste(missing_files, collapse = "\n"))
}

observed <- read_csv(observed_file, show_col_types = FALSE)
viper_null <- read_csv(viper_null_file, show_col_types = FALSE)
fgsea_null <- read_csv(fgsea_null_file, show_col_types = FALSE)
size_observed <- read_csv(size_file, show_col_types = FALSE)

size_benchmark <- size_observed %>%
  filter(str_detect(collection, "^b[0-9]+$")) %>%
  mutate(
    null_collection = paste0(Context, "_", collection),
    GeneSetName = new_id
  ) %>%
  left_join(
    fgsea_null %>%
      rename(null_collection = Collection),
    by = c(
      "null_collection",
      "GSE" = "GEO",
      "GENE" = "TargetGene",
      "GeneSetName"
    )
  ) %>%
  mutate(
    z_score = if_else(
      !is.na(grand_sd_NES) & grand_sd_NES > 0,
      (NES - grand_mean_NES) / grand_sd_NES,
      NA_real_
    ),
    confidence = as.numeric(str_match(pathway.x, "_c([0-9]+)$")[, 2])
  ) %>%
  group_by(GENE, Context.x) %>%
  filter(is.na(confidence) | confidence == max(confidence, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    score = case_when(
      expected_nes_sign == "pos" ~ z_score,
      expected_nes_sign == "neg" ~ -z_score,
      TRUE ~ NA_real_
    ),
    label = as.integer(random_flag == 1)
  ) %>%
  filter(!is.na(score))

gene_set_sizes <- c(25, 50, 100, 150, 200, 250, 300)
size_methods <- paste0("b", gene_set_sizes)
size_scores <- split(size_benchmark$score, size_benchmark$collection)[size_methods]
size_labels <- split(size_benchmark$label, size_benchmark$collection)[size_methods]

if (any(lengths(size_scores) != 218L)) {
  stop("Each gene-set-size benchmark must contain exactly 218 observations.")
}

size_pr_object <- evalmod(
  scores = size_scores,
  labels = size_labels,
  modnames = as.character(gene_set_sizes),
  mode = "rocprc"
)

size_auc_table <- auc(size_pr_object) %>%
  filter(curvetypes == "PRC") %>%
  transmute(gene_set_size = as.integer(modnames), AUCPR = aucs) %>%
  arrange(gene_set_size)

size_plot <- ggplot(
  fortify(size_pr_object) %>% filter(curvetype == "PRC"),
  aes(x = x, y = y, color = modname)
) +
  geom_line(linewidth = 1) +
  labs(x = "Recall", y = "Precision", color = "Gene-set size") +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "bottom"
  )

write_csv(
  size_benchmark,
  file.path(results_dir, "Source data_Fig 3d_gene_set_size_zscores.csv")
)
write_csv(
  size_auc_table,
  file.path(results_dir, "Source data_Fig 3d_gene_set_size_AUCPR.csv")
)
ggsave(
  file.path(results_dir, "Fig 3d_VMPA_PR_gene_set_size.pdf"),
  size_plot,
  width = 6.5,
  height = 4.5,
  bg = "white"
)

normalize_context <- function(x) {
  x <- trimws(x)
  recode(
    x,
    "CRC" = "crc",
    "Glioblastoma" = "glioma",
    "PDAC" = "pdac",
    "Breast Cancer" = "breast",
    "Ovarian Cancer" = "ovarian",
    "Melanoma" = "melanoma",
    "Prostate Cancer" = "prostate",
    "NSCLC" = "nsclc",
    .default = tolower(x)
  )
}

expected_methods <- c(
  "Collectri_ULM",
  "SigCom250",
  "VIPER_default",
  "VIPER_noFilter"
)

observed <- observed %>%
  mutate(
    Context = normalize_context(Context),
    GENE = toupper(GENE),
    method = case_when(
      str_detect(collection, "_b250$") ~ "VMPA_250",
      TRUE ~ collection
    )
  )

method_counts <- observed %>%
  count(method, random_flag, name = "n")

if (!setequal(unique(observed$method), c("VMPA_250", expected_methods))) {
  stop("Unexpected methods in the observed common-intersection table.")
}
if (any(method_counts$n != 64L)) {
  stop("Each method and class must contain exactly 64 observations.")
}

viper_rows <- observed %>%
  filter(method %in% c("VIPER_default", "VIPER_noFilter")) %>%
  left_join(
    viper_null %>%
      transmute(
        GSE,
        GENE = toupper(GENE),
        Context = normalize_context(Context),
        method = collection,
        grand_mean_NES,
        grand_sd_NES
      ),
    by = c("GSE", "GENE", "Context", "method")
  )

fgsea_null <- fgsea_null %>%
  transmute(
    GSE = GEO,
    GENE = toupper(TargetGene),
    Context = normalize_context(Context),
    collection = Collection,
    GeneSetName,
    grand_mean_NES,
    grand_sd_NES
  )

vmpa_rows <- observed %>%
  filter(method == "VMPA_250") %>%
  mutate(GeneSetName = str_remove(pathway, "_c[0-9]+$")) %>%
  left_join(
    fgsea_null,
    by = c("GSE", "GENE", "Context", "collection", "GeneSetName")
  )

comparator_rows <- observed %>%
  filter(method %in% c("SigCom250", "Collectri_ULM")) %>%
  left_join(
    fgsea_null %>% select(-GeneSetName),
    by = c("GSE", "GENE", "Context", "collection")
  )

benchmark <- bind_rows(vmpa_rows, comparator_rows, viper_rows) %>%
  mutate(
    z_score = if_else(
      !is.na(grand_sd_NES) & grand_sd_NES > 0,
      (NES - grand_mean_NES) / grand_sd_NES,
      NA_real_
    ),
    score = case_when(
      expected_nes_sign == "pos" ~ z_score,
      expected_nes_sign == "neg" ~ -z_score,
      TRUE ~ NA_real_
    ),
    label = as.integer(random_flag == 1)
  )

if (anyNA(benchmark$score)) {
  stop("Missing null-model statistics or perturbation directions after merging.")
}

benchmark_counts <- benchmark %>%
  count(method, label, name = "n")
if (any(benchmark_counts$n != 64L)) {
  stop("The merged benchmark does not contain 64 real and 64 null scores per method.")
}

score_list <- split(benchmark$score, benchmark$method)
label_list <- split(benchmark$label, benchmark$method)

pr_object <- evalmod(
  scores = score_list,
  labels = label_list,
  modnames = names(score_list),
  mode = "rocprc"
)

auc_table <- auc(pr_object) %>%
  filter(curvetypes == "PRC") %>%
  transmute(method = modnames, AUCPR = aucs) %>%
  arrange(desc(AUCPR))

method_colors <- c(
  "VMPA_250" = "#d73027",
  "VIPER_default" = "#1f78b4",
  "VIPER_noFilter" = "#9ad0ec",
  "SigCom250" = "#2e7d32",
  "Collectri_ULM" = "#9ccc65"
)

method_labels <- c(
  "VMPA_250" = "VMPA",
  "VIPER_default" = "VIPER (default)",
  "VIPER_noFilter" = "VIPER (no filter)",
  "SigCom250" = "SigCom",
  "Collectri_ULM" = "CollecTRI-ULM"
)

plot_data <- fortify(pr_object) %>%
  filter(curvetype == "PRC")

benchmark_plot <- ggplot(
  plot_data,
  aes(x = x, y = y, color = modname)
) +
  geom_line(linewidth = 1.2) +
  labs(x = "Recall", y = "Precision", color = "Method") +
  scale_color_manual(values = method_colors, labels = method_labels) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "bottom"
  )

write_csv(
  benchmark,
  file.path(results_dir, "Source data_Fig 3d_zscores.csv")
)
write_csv(
  auc_table,
  file.path(results_dir, "Source data_Fig 3d_AUCPR.csv")
)
ggsave(
  file.path(results_dir, "Fig 3d_VMPA_PR_benchmark.pdf"),
  benchmark_plot,
  width = 6.5,
  height = 4.5,
  bg = "white"
)

print(auc_table)
