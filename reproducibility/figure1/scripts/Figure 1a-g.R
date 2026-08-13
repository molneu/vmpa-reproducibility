#!/usr/bin/env Rscript
# Figure 1 Pipeline: CMap LV5 ➔ AKT-specific & cancer-driver signatures
# Author: Igor Cima
# Date:   2025-07-05


# ====== LIBRARIES ======
library(cmapR)
library(here)
library(httr)
library(zen4R)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(fgsea)
library(data.table)
library(org.Hs.eg.db)
library(msigdbr)
library(stringr)


# ====== SETUP ======
fig1_dir    <- here("reproducibility","figure1","data")
results_dir <- here("reproducibility","figure1","results")
dir.create(fig1_dir,    recursive=TRUE, showWarnings=FALSE)
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

# ====== STEP 1: Check data folder ======
file_names <- c(
  "level5_beta_trt_xpr_n142901x12328.gctx",
  "siginfo_beta.txt",
  "250120 AKT1-related gene sets in MsigDB.csv",
  "250603 AKT_related gene sets_MSigDB.csv",
  "250604 multiple_TARGETS-related gene sets in MsigDB.csv",
  "LINC gene annotations.csv"
)

missing_files <- file_names[!file.exists(file.path(fig1_dir, file_names))]
if (length(missing_files) > 0) {
  stop(
    "Missing Figure 1 input files in ", fig1_dir, ": ",
    paste(missing_files, collapse = ", "),
    ". Download the Figshare data and place the Figure 1 inputs in this data folder."
  )
}

# ====== STEP 2: Load & Annotate Column Metadata ======
gctx_file    <- file.path(fig1_dir, "level5_beta_trt_xpr_n142901x12328.gctx")
siginfo_file <- file.path(fig1_dir, "siginfo_beta.txt")

trt_xpr   <- parse_gctx(gctx_file)
meta_full <- read.delim(siginfo_file)
meta_xpr  <- meta_full %>% filter(pert_type=="trt_xpr")
annotated <- annotate_gct(trt_xpr, meta_xpr, dim="column", keyfield="sig_id")

rm(trt_xpr, meta_full)

# ====== STEP 3: Annotate Rows ======
gene_map   <- read.csv(file.path(fig1_dir,"LINC gene annotations.csv"), stringsAsFactors=FALSE)
entrez2sym <- setNames(gene_map$gene.name, gene_map$entrez)
infer_map  <- setNames(gene_map$inference, gene_map$gene.name)

rids <- annotated@rid
annotated@rdesc$symbol    <- unname(entrez2sym[rids])
annotated@rdesc$inference <- unname(infer_map[annotated@rdesc$symbol])

# ====== STEP 3b: Build two filtered objects ======

# -- for broader pipeline (TAS > 0.15), s. STEP16, fig 1e,f --
keep     <- annotated@cdesc$tas > 0.15
a015     <- subset_gct(annotated, cid=which(keep))

rm(annotated)
gc()

keep     <- !is.na(a015@cdesc$cmap_name) & a015@cdesc$cmap_name!=""
a015     <- subset_gct(a015, cid=which(keep))

rm_idx   <- which(a015@cdesc$cmap_name=="UnTrt" | grepl("^BRDN0", a015@cdesc$cmap_name))
if (length(rm_idx)) {
  keep2   <- setdiff(seq_along(a015@cdesc$cmap_name), rm_idx)
  a015    <- subset_gct(a015, cid=keep2)
}

# -- for AKT1 pipeline (TAS > 0.2) --
keep   <- a015@cdesc$tas > 0.2
a02    <- subset_gct(a015, cid=which(keep))

keep   <- !is.na(a02@cdesc$cmap_name) & a02@cdesc$cmap_name!=""
a02    <- subset_gct(a02, cid=which(keep))

# change column id annotation to include gene name
a02@cdesc$id <- paste0(a02@cdesc$cmap_name, "_", a02@cdesc$sig_id)

# ====== STEP 4: Subset AKT1 Signatures ======
df_mat    <- as.data.frame(a02@mat)
syms02    <- a02@rdesc$symbol
rownames(df_mat) <- make.unique(syms02)
colnames(df_mat) <- a02@cdesc$id

akt1_ids  <- a02@cdesc %>% as_tibble() %>% filter(cmap_name=="AKT1") %>% pull(id)
akt1_sigs <- df_mat[, colnames(df_mat) %in% akt1_ids]

# ====== STEP 5: Select Top AKT1-Downregulated =====
lowest_is_akt1   <- apply(akt1_sigs, 2, function(col)
  which.min(col) == which(rownames(akt1_sigs)=="AKT1")
)
akt1_sigs_TOP    <- akt1_sigs[, lowest_is_akt1]
message("Selected signatures with AKT1 lowest, n=", ncol(akt1_sigs_TOP))

# ====== STEP 7: Prepare Gene Sets =====
msig      <- msigdbr(species="Homo sapiens")
akt1_sets <- subset(msig, grepl("AKT1", gs_name))
akt_sets <- subset(msig, grepl("AKT", gs_name))
akt1_list <- split(akt1_sets$gene_symbol, akt1_sets$gs_name)
akt_list  <- split(akt_sets$gene_symbol,  akt_sets$gs_name)

# ====== STEP 8: Run fGSEA (AKT1 sets) =====
set.seed(1)
stats_list <- lapply(akt1_sigs_TOP, function(col) setNames(col, rownames(akt1_sigs_TOP)))
res_list   <- lapply(names(stats_list), function(s) {
  res<-fgseaMultilevel(
    pathways   = akt1_list,
    stats      = stats_list[[s]],
    minSize    = 1,
    maxSize    = 500,
    nPermSimple= 5000
  )
  res$signature <- s
  res
})
fgsea_all  <- bind_rows(res_list)

# ====== STEP 9: Add Expected NES & Classify =====
gs1_df    <- read.csv(file.path(fig1_dir,"250120 AKT1-related gene sets in MsigDB.csv"), row.names=1)
names(gs1_df)[ names(gs1_df) == "expected.NES.sign.in.AKT.KO.cells" ] <- "expected_NES_sign"

fgsea_all <- fgsea_all %>%
  left_join(
    gs1_df %>% dplyr::select(gs_name,
                             expected_NES_sign),
    by = c("pathway" = "gs_name")
  ) 

fgsea_cls <- fgsea_all %>%
  filter(!is.na(expected_NES_sign)) %>%
  mutate(result = case_when(
    padj > 0.1                                                ~ "FN",
    padj <= 0.1 & NES >  0 & expected_NES_sign=="neg"        ~ "FP",
    padj <= 0.1 & NES <  0 & expected_NES_sign=="pos"        ~ "FP",
    padj <= 0.1 & NES >  0 & expected_NES_sign=="pos"        ~ "TP",
    padj <= 0.1 & NES <  0 & expected_NES_sign=="neg"        ~ "TP",
    TRUE                                                      ~ NA_character_
  ))

fgsea_wide <- fgsea_cls %>%
  dplyr::select(pathway, signature, result) %>%
  pivot_wider(names_from=signature, values_from=result)

fgsea_sum <- fgsea_wide %>%
  rowwise() %>%
  mutate(
    FN = sum(c_across(where(is.character))=="FN", na.rm=TRUE),
    FP = sum(c_across(where(is.character))=="FP", na.rm=TRUE),
    TP = sum(c_across(where(is.character))=="TP", na.rm=TRUE),
    Precision = ifelse((TP+FP)>0, TP/(TP+FP), 0),
    Recall    = ifelse((TP+FN)>0, TP/(TP+FN), 0),
    F1        = ifelse((Precision+Recall)>0,
                       2*(Precision*Recall)/(Precision+Recall), 0)
  ) %>%
  ungroup()

# ====== STEP 10: Plot Fig 1a−c =====

###>>>> Fig 1a, Dotplot representative GSEA <<<<  ----

fgsea_U251   <- fgsea_all %>%
  filter(signature=="AKT1_XPR015_U251MG.311_96H:E14", !is.na(expected_NES_sign))

fgsea_U251_o <- fgsea_U251[order(fgsea_U251$NES, decreasing=TRUE),]
fgsea_U251_o$significant <- ifelse(fgsea_U251_o$padj < 0.05, "Significant","Not Significant")

dot_p <- ggplot(fgsea_U251_o, aes(
  x = reorder(pathway, NES), y = NES,
  size= size, color=significant
)) +
  geom_point() +
  coord_flip() +
  labs(x="Gene Set", y="Normalized Enrichment Score (NES)") +
  ggtitle("XPR015_U251MG.311_96H:E14") +
  theme_minimal() +
  theme(
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size=14),
    axis.text            = element_text(size=9),
    axis.line            = element_line(color="black", size=0.5),
    axis.ticks           = element_line(color="black", size=0.5),
    axis.ticks.length    = unit(0.15,"cm"),
    axis.text.x          = element_text(angle=0,hjust=0.5),
    legend.text          = element_text(size=14),
    legend.title         = element_text(size=14)
  ) +
  scale_color_manual(values=c("Not Significant"="gray","Significant"="red"))

ggsave(
  file.path(results_dir,"Fig 1a_Representative_fGSEA_result.pdf"),
  dot_p, width=8, height=3, bg="white"
)

# View(fgsea_sum)

###>>>> Fig 1b, Dotplot Prec-Rec AKT1 <<<<  ----
p_pr <- ggplot(fgsea_sum, aes(x=Recall,y=Precision)) +
  geom_jitter(size=8, width=0.05, height=0.05, shape=1, color="black") +
  scale_x_continuous(limits=c(-0.1,1.1),expand=c(0,0)) +
  scale_y_continuous(limits=c(-0.1,1.1),expand=c(0,0)) +
  labs(x="Recall", y="Precision") +
  theme_minimal() +
  theme(
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size=14),
    axis.text          = element_text(size=14),
    axis.line          = element_line(color="black", size=0.5),
    axis.ticks         = element_line(color="black", size=0.5),
    axis.ticks.length  = unit(0.25,"cm"),
    axis.text.x        = element_text(angle=45,hjust=1),
    legend.text        = element_text(size=14),
    legend.title       = element_text(size=14)
  )

ggsave(
  file.path(results_dir,"Fig 1b_PR plot_AKT1.pdf"),
  p_pr, width=5, height=5, bg="white"
)

###>>>> Fig 1c, Pie chart  FP/FN/TP AKT1 <<<<  ----

# ensure Category factor order if you care
tb <- fgsea_sum %>%
  summarise(TP = sum(TP), FP = sum(FP), FN = sum(FN)) %>%
  pivot_longer(everything(), names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("FN", "FP", "TP")))

# compute percentages for labels (optional)
tb <- tb %>%
  mutate(
    pct = Count / sum(Count),
    lbl = paste0(Category, ": ", scales::percent(pct, accuracy = 0.1))
  )

# pie chart
p_pie <- ggplot(tb, aes(x = "", y = Count, fill = Category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(TP = "#619CFF", FP = "#F8766D", FN = "grey")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank()
  ) +
  labs(title = "Confusion Components") +
  geom_text(
    aes(label = lbl),
    position = position_stack(vjust = 0.5),
    size = 3
  )


ggsave(
  file.path(results_dir,"Fig 1c_Pie_AKT1.pdf"),
  p_pie, width=4, height=5, bg="white"
)

# ====== STEPS 11–15: repeat analogous blocks for broader AKT gene sets =====

# perform same computations with akt_gene_sets_list instead of akt1_gene_sets_list
gs2_df <- read.csv(file.path(fig1_dir,"250603 AKT_related gene sets_MSigDB.csv"),row.names=1)

# Ensure reproducibility
set.seed(1)

# Convert filtered AKT signatures to named statistic vectors
dge_signatures <- lapply(akt1_sigs_TOP, function(col) {
  stats <- setNames(col, rownames(akt1_sigs_TOP))
  stats[!is.na(stats)]
})

# Run fgseaMultilevel for each signature
set.seed(1)
fgsea_results_list <- lapply(names(dge_signatures), function(sig_name) {
  res <- fgseaMultilevel(
    pathways = akt_list,
    stats    = dge_signatures[[sig_name]],
    minSize  = 1,
    maxSize  = 500,
    nPermSimple = 5000
  )
  res$signature <- sig_name
  res
})

fgsea_all_akt <- bind_rows(fgsea_results_list)
fgsea_all_akt <- fgsea_all_akt[, -8]


# merge expected signs
fgsea_all_akt <- fgsea_all_akt %>%
  left_join(gs2_df %>% dplyr::select(gs_name,`expected.NES.sign.in.AKT.KO.cells`), by=c("pathway"="gs_name")) %>%
  dplyr::rename(expected_NES_sign=`expected.NES.sign.in.AKT.KO.cells`)

# classify
fgsea_cls_akt <- fgsea_all_akt %>% filter(!is.na(expected_NES_sign)) %>% mutate(result=case_when(
  padj>0.1~"FN",
  padj<=0.1&NES>0&expected_NES_sign=="neg"~"FP",
  padj<=0.1&NES<0&expected_NES_sign=="pos"~"FP",
  padj<=0.1&NES>0&expected_NES_sign=="pos"~"TP",
  padj<=0.1&NES<0&expected_NES_sign=="neg"~"TP",
  TRUE~NA_character_
))

# wide
fgsea_wide <- fgsea_cls_akt %>% dplyr::select(pathway,signature,result) %>% pivot_wider(names_from=signature,values_from=result)

# summary
fgsea_sum_akt <- fgsea_wide %>% rowwise() %>% mutate(
  FN=sum(c_across(where(is.character))=="FN",na.rm=TRUE),
  FP=sum(c_across(where(is.character))=="FP",na.rm=TRUE),
  TP=sum(c_across(where(is.character))=="TP",na.rm=TRUE),
  Precision=ifelse((TP+FP)>0,TP/(TP+FP),0),
  Recall   =ifelse((TP+FN)>0,TP/(TP+FN),0),
  F1       =ifelse((Precision+Recall)>0,2*(Precision*Recall)/(Precision+Recall),0)
) %>% ungroup()

# classify results based on padj 0.1 threshold
fgsea_classified_akt <- fgsea_all_akt %>%
  filter(!is.na(expected_NES_sign)) %>%
  mutate(
    result = case_when(
      padj > 0.1 ~ "FN",
      padj <= 0.1 & NES > 0 & expected_NES_sign == "neg" ~ "FP",
      padj <= 0.1 & NES < 0 & expected_NES_sign == "pos" ~ "FP",
      padj <= 0.1 & NES > 0 & expected_NES_sign == "pos" ~ "TP",
      padj <= 0.1 & NES < 0 & expected_NES_sign == "neg" ~ "TP",
      TRUE ~ NA_character_
    )
  )


fgsea_wide <- fgsea_classified_akt %>%
  dplyr::select(pathway, signature, result) %>%
  pivot_wider(names_from = signature, values_from = result)

# compute summary metrics
fgsea_summary_akt <- fgsea_wide %>%
  rowwise() %>%
  mutate(
    FN = sum(c_across(where(is.character)) == "FN", na.rm = TRUE),
    FP = sum(c_across(where(is.character)) == "FP", na.rm = TRUE),
    TP = sum(c_across(where(is.character)) == "TP", na.rm = TRUE),
    Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
    Recall    = ifelse((TP + FN) > 0, TP / (TP + FN), 0),
    F1 = ifelse((Precision + Recall) > 0,
                2 * (Precision * Recall) / (Precision + Recall),
                0)
  ) %>%
  ungroup()

# save summary
write.csv(
  fgsea_summary_akt,
  file.path(results_dir, "Source data Fig 1d,e:250603_AKT_pathway_fGSEA_classification_and_PR.csv"),
  row.names = FALSE
)

###>>>> Fig 1d, Dotplot Prec-Rec Plot AKT <<<<  ----

p_pr_akt <- ggplot(fgsea_summary_akt, aes(x = Recall, y = Precision)) +
  geom_jitter(size = 8, width = 0.05, height = 0.05, shape = 1, color = "black") +
  scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0,0)) +
  labs(x = "Recall", y = "Precision") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 14),
    axis.line  = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(results_dir, "Fig 1d_PR plot_AKT.pdf"),
  plot     = p_pr_akt,
  width    = 5,
  height   = 5,
  device   = "pdf",
  bg       = "white"
)

###>>>> Fig 1e, Pie chart FP/FN/TP AKT <<<<  ----

# order factors
tb_a <- fgsea_summary_akt %>%
  summarise(TP = sum(TP), FP = sum(FP), FN = sum(FN)) %>%
  pivot_longer(everything(), names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("FN", "FP", "TP")))

# compute percentages for labels 
tb_a <- tb_a %>%
  mutate(
    pct = Count / sum(Count),
    lbl = paste0(Category, ": ", scales::percent(pct, accuracy = 0.1))
  )

# pie chart
p_pie <- ggplot(tb_a, aes(x = "", y = Count, fill = Category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(TP = "#619CFF", FP = "#F8766D", FN = "grey")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank()
  ) +
  labs(title = "Confusion Components") +
  geom_text(
    aes(label = lbl),
    position = position_stack(vjust = 0.5),
    size = 3
  )


ggsave(
  file.path(results_dir,"Fig 1e_Pie_AKT.pdf"),
  p_pie, width=4, height=5, bg="white"
)

# ====== STEP 16: Additional cancer-driver signatures (Fig 1f & 1g) =====

#free up memory
rm(a02)
gc()

gs3_df <- read.csv(file.path(fig1_dir,"250604 multiple_TARGETS-related gene sets in MsigDB.csv"),row.names=1)


#selected common cancer drivers
target_genes <- c("AKT1","ERBB2","KRAS","TP53","CCND3","CDK4","EGFR","NF1","RHOA","STAT3")

#keep only targets
keep_targets    <- a015@cdesc$cmap_name %in% target_genes
a015 <- subset_gct(a015, cid = which(keep_targets))

#keep only exemplar signatures
keep_exemplars     <- a015@cdesc$is_exemplar_sig == 1
a015     <- subset_gct(a015, cid = which(keep_exemplars))

# -- change column id annotation to include gene name
a015@cdesc$id <- paste0(a015@cdesc$cmap_name, "_", a015@cdesc$sig_id)


#build expression matrix from a015 (TAS>0.15)
expr_all <- as.data.frame(a015@mat)
syms015 <- a015@rdesc$symbol
rownames(expr_all) <- make.unique(syms015)
col_ids <- a015@cdesc$id
colnames(expr_all) <- col_ids

# select signatures where their own gene z-score < –1.5
sig2gene <- setNames(a015@cdesc$cmap_name, col_ids)

keep2 <- vapply(col_ids, function(id) {
  tgt <- sig2gene[[id]]
  sc  <- expr_all[tgt, id]
  !is.na(sc) && sc < -1.5
}, logical(1))

expr_filt <- expr_all[, keep2, drop = FALSE]

# tabulate number of signatures for each target:
gene_prefix <- sub("_.*", "", colnames(expr_filt))
sig_counts <- as.data.frame(table(gene_prefix), stringsAsFactors = FALSE)
colnames(sig_counts) <- c("gene", "n_signatures")
print(sig_counts)

# KRAS ans NF1 have only one signature each, exclude, because difficult to estimate precision/recall on n=1)
# in Methods: report on chosing common cancer drivers with >1 signature available

expr_filt <- expr_filt %>%
  dplyr::select(-starts_with("KRAS_"), -starts_with("NF1_"))

# 16.4. save for downstream FGSEA
write.csv(
  expr_filt,
  file.path(results_dir,"Source data Fig1 f,g_cancer_drivers_expr_zlt-1.5_tas0.15.csv"),
  row.names=TRUE
)

# fgsea
set.seed(1)

# prepare signatures
dge_signatures <- lapply(expr_filt, function(col) {
  stats <- setNames(col, rownames(expr_filt))
  stats[!is.na(stats)]  # remove NAs if any
})

# prepare gene sets
# For paper,  running only H, CGP (C2), CP(C2) and C6
h_gene_sets_H = msigdbr(species = "Homo sapiens", collection = "H") # hallmark gene sets
h_gene_sets_CP = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP") # canonical pathways
h_gene_sets_CGP = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP") # chemical and genetic perturbations
h_gene_sets_C6  = msigdbr(species = "Homo sapiens", collection = "C6") 

msigdbr_H = split(x = h_gene_sets_H$gene_symbol, f = h_gene_sets_H$gs_name)
msigdbr_CP = split(x = h_gene_sets_CP$gene_symbol, f = h_gene_sets_CP$gs_name)
msigdbr_CGP = split(x = h_gene_sets_CGP$gene_symbol, f = h_gene_sets_CGP$gs_name)
msigdbr_C6 = split(x = h_gene_sets_C6$gene_symbol, f = h_gene_sets_C6$gs_name)

msigdbr_H_sel <- msigdbr_H[grepl("CCND3|AKT1|P53|CDK4|EGFR|ERBB2|RHOA|STAT3", names(msigdbr_H), ignore.case = TRUE)]
msigdbr_CGP_sel <- msigdbr_CGP[grepl("CCND3|AKT1|P53|CDK4|EGFR|ERBB2|RHOA|STAT3", names(msigdbr_CGP), ignore.case = TRUE)]
msigdbr_CP_sel <- msigdbr_CP[grepl("CCND3|AKT1|P53|CDK4|EGFR|ERBB2|RHOA|STAT3", names(msigdbr_CP), ignore.case = TRUE)]
msigdbr_C6_sel <- msigdbr_C6[grepl("CCND3|AKT1|P53|CDK4|EGFR|ERBB2|RHOA|STAT3", names(msigdbr_C6), ignore.case = TRUE)]

combined_gene_sets <- c(msigdbr_H_sel, msigdbr_CGP_sel, msigdbr_CP_sel, msigdbr_C6_sel)
length(combined_gene_sets)

# Run fgseaMultilevel for each signature and store results with an ID
set.seed(1)
fgsea_results_list <- lapply(names(dge_signatures), function(sig_name) {
  stats <- dge_signatures[[sig_name]]
  res <- fgseaMultilevel(
    pathways = combined_gene_sets,
    stats = stats,
    minSize = 1,
    maxSize = 500,
    nPermSimple = 5000
  )
  res$signature <- sig_name  # add column for signature ID
  res
})

fgsea_all_results <- bind_rows(fgsea_results_list)

fgsea_all_results <- fgsea_all_results[, -8]

# add a column with target gene

fgsea_all_results <- fgsea_all_results %>% mutate(sig_target = str_extract(signature, "^[^_]+"))

#add column with input gene ID related to gene sets

target_genes <- c("CCND3", "AKT1", "TP53", "CDK4", "EGFR", "ERBB2", "RHOA", "STAT3")

# Detect presence of gene name as substring (case-insensitive)
fgsea_all_results$sig_input <- sapply(fgsea_all_results$pathway, function(pw) {
  matches <- sapply(target_genes, function(gene) {
    grepl(gene, pw, ignore.case = TRUE)
  })
  matched_genes <- target_genes[matches]
  if (length(matched_genes) == 0) {
    return(NA)
  } else {
    return(paste(matched_genes, collapse = ";"))
  }
})

fgsea_all_results$sig_input <- gsub("\\bP53\\b", "TP53", fgsea_all_results$sig_input) #change P53 to TP53 in the sig_input column

#Filter rows from fgsea_all_results where sig_input == sig_target and combine all those rows into one new dataframe.
matched_rows_df <- fgsea_all_results %>%
  filter(sig_input == sig_target)

#Add expected NES sign column based on manually curated gene set analysis
matched_rows_df <- matched_rows_df %>%
  left_join(gs3_df, by = "pathway")

fgsea_classified <- matched_rows_df %>%
  filter(!is.na(expected_NES_sign)) %>%
  mutate(
    result = case_when(
      padj > 0.1 ~ "FN",
      padj <= 0.1 & NES > 0 & expected_NES_sign == "neg" ~ "FP",
      padj <= 0.1 & NES < 0 & expected_NES_sign == "pos" ~ "FP",
      padj <= 0.1 & NES > 0 & expected_NES_sign == "pos" ~ "TP",
      padj <= 0.1 & NES < 0 & expected_NES_sign == "neg" ~ "TP",
      TRUE ~ NA_character_
    )
  )

print("Number of CRISPR perturbation experiments used for multiple cancer driver figure")
length(unique(fgsea_classified$signature))

#summarize results

fgsea_pathway_summary <- fgsea_classified %>%
  filter(!is.na(sig_target)) %>%
  group_by(pathway) %>%
  summarise(
    sig_target = dplyr::first(sig_target),  # assumes one sig_target per pathway
    TP = sum(result == "TP", na.rm = TRUE),
    FP = sum(result == "FP", na.rm = TRUE),
    FN = sum(result == "FN", na.rm = TRUE),
    
    Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
    Recall    = ifelse((TP + FN) > 0, TP / (TP + FN), 0),
    
    F1 = ifelse((Precision + Recall) > 0,
                2 * (Precision * Recall) / (Precision + Recall),
                0),
    .groups = "drop"
  )

#per gene

fgsea_gene_summary <- fgsea_classified %>%
  group_by(sig_target) %>%
  summarise(
    TP = sum(result == "TP", na.rm = TRUE),
    FP = sum(result == "FP", na.rm = TRUE),
    FN = sum(result == "FN", na.rm = TRUE),
    
    Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
    Recall    = ifelse((TP + FN) > 0, TP / (TP + FN), 0),
    
    F1 = ifelse((Precision + Recall) > 0,
                2 * (Precision * Recall) / (Precision + Recall),
                0),
    .groups = "drop"
  )


###>>>> Fig 1f, Dotplot Prec-Rec Plot Multiple cancer genes <<<<  ----

# Select the relevant columns (Precision, Recall, and Gene)
data_cpg <- fgsea_pathway_summary[,c(2,6,7)]  # Column 2 is the "Gene" column
table(data_cpg$sig_target)

# Create the plot with colored dots and black outlines
p_pr_multi <- ggplot(data_cpg, aes(x = Recall, y = Precision)) +
  # Add colored dots with black outlines
  geom_jitter(data = head(data_cpg, -1), aes(fill = sig_target), size = 8, width = 0.03, height = 0.03, shape = 21, color = "black") +
  scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  labs(x = "Recall", y = "Precision", fill = "Gene") +  # Use `fill` for the legend
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 14),  # Adjust legend text size
    legend.title = element_text(size = 14)  # Adjust legend title size
  )

#save
ggsave(
  filename = file.path(results_dir, "Fig 1f_PR Plot_Multiple Cancer genes.pdf"),
  plot     = p_pr_multi,
  width    = 7,
  height   = 6,
  device   = "pdf",
  bg       = "white"
)


# Total number of gene sets

n_gene_sets <- nrow(data_cpg)

# Print the result
print(paste("Number of gene sets analyzed:", n_gene_sets))


# Filter rows where Recall = 0 and Precision = 0
zero_recall_precision <- data_cpg %>%
  filter(Recall == 0 & Precision == 0)

# Count the number of rows
num_zero_recall_precision <- nrow(zero_recall_precision)

# Print the result
print(paste("Number of data points with 0 Recall and 0 Precision:", num_zero_recall_precision))

# Count data points with Recall < 0.25
low_recall <- data_cpg %>%
  filter(Recall < 0.25)

num_low_recall <- nrow(low_recall)

# Count data points with both Precision < 0.25 and Recall < 0.25
low_precision_recall <- data_cpg %>%
  filter(Precision < 0.25 & Recall < 0.25)

num_low_precision_recall <- nrow(low_precision_recall)

# Print the results
print(paste("Number of data points with Recall < 0.25:", num_low_recall))
print(paste("Number of data points with Precision < 0.25 and Recall < 0.25:", num_low_precision_recall))


###>>>> Fig 1g, Stacked Plot FP/FN/TP Multiple cancer genes <<<<  ----


# Reshape the data for ggplot2
plot_data <- pivot_longer(fgsea_gene_summary, 
                          cols = c(TP, FP, FN), 
                          names_to = "Category", 
                          values_to = "Count")

# Create the stacked barplot
p_bar_multi <- ggplot(plot_data, aes(x = sig_target, y = Count, fill = Category)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("TP" = "#619CFF", "FP" = "#F8766D", "FN" = "grey")) +
  labs(
    y = "Percentage",
    x = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    #axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0))

#save
ggsave(
  filename = file.path(results_dir, "Fig 1g_StackedBar_Multiple Cancer genes.pdf"),
  plot     = p_bar_multi,
  width    = 6,
  height   = 6,
  device   = "pdf",
  bg       = "white"
)
