# ====== STEP 8: Random Forest validation ======

library(here)
library(caret)
library(ggplot2)

results_dir <- here("reproducibility", "figure5", "results")
supp_results_dir <- here("reproducibility", "supplementary_figure7", "results")
dir.create(supp_results_dir, recursive = TRUE, showWarnings = FALSE)

gsva_results_list <- readRDS(file.path(results_dir, "STEP 4 output_gsva_results_list.rds"))

# Use the GSVA scores from compass (log_1.5 subset)
compass250_gsva_scores <- gsva_results_list[["log_1.5"]]

# --- Step A: Bianca scaling ---
N <- nrow(compass250_gsva_scores)
compass250_gsva_scores_bianca <- scale(
  compass250_gsva_scores,
  scale = apply(compass250_gsva_scores, 2, sd) * sqrt((N - 1) / N)
)

# --- Step B: Remove highly correlated rows ---
cor_matrix <- cor(t(compass250_gsva_scores_bianca))
high_cor <- findCorrelation(cor_matrix, cutoff = 0.8)
compass250_gsva_scores_bianca.corr <- compass250_gsva_scores_bianca[-high_cor, ]

# --- Step C: Select top 10% variance features ---
num_genes <- nrow(compass250_gsva_scores_bianca.corr)
top_10_percent <- ceiling(num_genes * 0.10)
top_variance_genes <- head(order(apply(compass250_gsva_scores_bianca.corr, 1, var), decreasing = TRUE), top_10_percent)

compass250_gsva_scores_scaled_TOP <- compass250_gsva_scores_bianca.corr[top_variance_genes, ]

##### IV) -> Random Forest Validation of Clusters (100x) ----

library(caret)
library(randomForest)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)

# ---- Function to run RF validation ----
rf_validation <- function(data_matrix, cluster_labels, n_iter = 10, ntree = 500, seed = 42) {
  set.seed(seed)
  
  training_df <- data.frame(t(data_matrix), cluster = as.factor(cluster_labels))
  results_list <- list()
  
  for (i in 1:n_iter) {
    # Train/test split
    train_index <- createDataPartition(training_df$cluster, p = 0.8, list = FALSE)
    train_set <- training_df[train_index, ]
    test_set  <- training_df[-train_index, ]
    
    # Train RF
    rf_model <- randomForest(cluster ~ ., data = train_set, ntree = ntree, importance = TRUE)
    
    # Predict
    test_predictions <- predict(rf_model, test_set)
    
    # Confusion matrix
    cm <- confusionMatrix(test_predictions, test_set$cluster)
    
    # F1 scores
    precision <- cm$byClass[, "Precision"]
    recall    <- cm$byClass[, "Recall"]
    f1_score  <- 2 * (precision * recall) / (precision + recall)
    
    # Collect metrics
    cluster_metrics <- as.data.frame(cm$byClass) %>%
      mutate(
        Cluster   = rownames(.),
        Iteration = i,
        `F1_Score` = f1_score
      ) %>%
      dplyr::select(Iteration, Cluster, Sensitivity, Specificity, `Balanced Accuracy`, `F1_Score`)
    
    results_list[[i]] <- cluster_metrics
  }
  
  results_df <- bind_rows(results_list)
  results_df$Cluster <- factor(results_df$Cluster, levels = unique(results_df$Cluster))
  
  # Summarise medians
  medians_df <- results_df %>%
    group_by(Cluster) %>%
    summarise(
      Sensitivity       = round(median(Sensitivity, na.rm = TRUE), 2),
      Specificity       = round(median(Specificity, na.rm = TRUE), 2),
      `Balanced Accuracy` = round(median(`Balanced Accuracy`, na.rm = TRUE), 2),
      F1_Score          = round(median(`F1_Score`, na.rm = TRUE), 2),
      .groups = "drop"
    )
  
  list(results_df = results_df, medians_df = medians_df)
}

# ---- Function to plot metrics ----
plot_metric <- function(results_df, medians_df, metric_name) {
  ggplot(results_df, aes(x = Cluster, y = .data[[metric_name]])) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste(metric_name, "Distribution Across Clusters"), 
         y = metric_name, x = "Cluster") +
    ylim(0, 1) +
    geom_text(data = medians_df, 
              aes(x = Cluster, y = 1, label = round(.data[[metric_name]], 2)), 
              vjust = -0.5, color = "blue")
}

# ---- Run validation for COMPASS ----

combined_m3c_results <- readRDS(file.path(results_dir, "STEP 5 output_M3C object_combined_m3c_results.rds"))
m3c_results_list.1 <- combined_m3c_results[[1]]


# Extract cluster assignments
ordered_annot <- m3c_results_list.1[["realdataresults"]][[11]][["ordered_annotation"]]
cluster_labels <- ordered_annot$consensuscluster
table(cluster_labels)

# Exclude clusters 4 and 7
clusters_to_remove <- c("4", "7")
keep_idx <- !(cluster_labels %in% clusters_to_remove)

filtered_data   <- compass250_gsva_scores_scaled_TOP[, keep_idx]
filtered_labels <- cluster_labels[keep_idx]
filtered_labels <- droplevels(filtered_labels)

# Run RF validation on filtered data
rf_gsva <- rf_validation(filtered_data, filtered_labels, n_iter = 100)

# Display median performance table ---
ggtexttable(rf_gsva$medians_df, rows = NULL,
            theme = ttheme("minimal", base_size = 7,
                           colnames.style = list(face = "bold")))


# Run RF validation with permuted labels ---
set.seed(42)
permuted_labels <- sample(filtered_labels)  # shuffle labels
rf_perm <- rf_validation(filtered_data, permuted_labels, n_iter = 100)
rf_perm$results_df[is.na(rf_perm$results_df)] <- 0

# Display median performance table
ggtexttable(rf_perm$medians_df, rows = NULL,
            theme = ttheme("minimal", base_size = 7,
                           colnames.style = list(face = "bold")))

# combine, visualize and test 

# Add type labels
rf_gsva$results_df$Type <- "Real"
rf_perm$results_df$Type <- "Permuted"
rf_perm$results_df[is.na(rf_perm$results_df)] <- 0

# Combine
all_results <- bind_rows(rf_gsva$results_df, rf_perm$results_df)

#plot

rfplot <- ggplot(all_results, aes(x = Cluster, y = F1_Score, fill = Type)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.8, width = 0.6,
               position = position_identity()) +
  scale_fill_manual(values = c("Real" = "#e31a1c", "Permuted" = "gray")) +  
  theme_minimal(base_size = 12) +
  labs(
    title = "Random Forest validation: Real vs Permuted labels",
    x = "Cluster",
    y = "F1 Score",
    fill = "Condition"
  ) +
  ylim(0, 1) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(color = "black"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggplot2::ggsave(
  filename = file.path(supp_results_dir, "Supplementary_Figure_7_random_forest_validation.pdf"),
  plot = rfplot,
  width = 5, height = 4, dpi = 300, device = "pdf"
)

# ====== STEP 8b: Statistical testing ======

# --- Global comparison (median F1 per iteration) ---
real_iter <- rf_gsva$results_df %>%
  group_by(Iteration) %>%
  summarise(medF1 = median(F1_Score, na.rm = TRUE), .groups = "drop")

perm_iter <- rf_perm$results_df %>%
  group_by(Iteration) %>%
  summarise(medF1 = median(F1_Score, na.rm = TRUE), .groups = "drop")

wilcox_global <- wilcox.test(real_iter$medF1, perm_iter$medF1)
cat("Global Wilcoxon test p =", wilcox_global$p.value, "\n")

# --- Cluster-wise tests ---
clusters <- unique(rf_gsva$results_df$Cluster)
pvals <- sapply(clusters, function(cl) {
  real_vals <- rf_gsva$results_df %>% filter(Cluster == cl) %>% pull(F1_Score)
  perm_vals <- rf_perm$results_df %>% filter(Cluster == cl) %>% pull(F1_Score)
  tryCatch(
    wilcox.test(real_vals, perm_vals)$p.value,
    error = function(e) NA
  )
})
pvals_fdr <- p.adjust(pvals, method = "fdr")

cluster_tests <- data.frame(
  Cluster = clusters,
  p_value = pvals,
  FDR = pvals_fdr
)

print(cluster_tests)
print(min(rf_gsva$medians_df$F1_Score))
print(max(rf_gsva$medians_df$F1_Score))
print(min(rf_perm$medians_df$F1_Score))
print(max(rf_perm$medians_df$F1_Score))

median(rf_gsva$medians_df$F1_Score)
median(rf_perm$medians_df$F1_Score)
