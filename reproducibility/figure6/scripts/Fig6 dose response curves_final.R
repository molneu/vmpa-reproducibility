#!/usr/bin/env Rscript

# Packages
library(drc)
library(dplyr)
library(gridExtra)
library(ggplotify)
library(here)

# Define folders
data_dir    <- here("reproducibility", "figure6", "data")
results_dir <- here("reproducibility", "figure6", "results")

# source data
files <- c(
  "250315_MTT_data_BN91.csv",
  "250618_MTT_data_BN118.csv",
  "250714_MTT_data_BN91.csv",
  "250714_MTT_data_BN118.csv",
  "250721_MTT_data_BN118.csv"
)

# Helper: draw your original plot
draw_main_plot <- function(df, fit2, main_title) {
  plot(fit2,
       type      = "average",
       pch       = c(16, 16),
       col       = c("#4575B4", "#D73027"),
       lty       = c(1, 1),
       lwd       = 1,
       cex       = 1.2,
       log       = "x",
       xlab      = "Dose (µM)",
       ylab      = "Viability (%)",
       main      = main_title,
       xlim      = c(0, max(df$Dose)),
       ylim      = c(0, max(fit2$coefficients + 15)),
       legendPos = c(0.4 * max(df$Dose, na.rm = TRUE),
                     0.20 * max(df$Response, na.rm = TRUE))
  )
  
  # Mean ± SEM
  sum_df <- aggregate(Response ~ Cell_line + Dose, data = df,
                      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                          sd = sd(x, na.rm = TRUE),
                                          n = length(x)))
  sum_df <- do.call(data.frame, sum_df)
  names(sum_df)[names(sum_df) == "Response.mean"] <- "mean"
  names(sum_df)[names(sum_df) == "Response.sd"]   <- "sd"
  names(sum_df)[names(sum_df) == "Response.n"]    <- "n"
  sum_df$sem <- with(sum_df, sd / sqrt(n))
  
  cell_levels <- unique(sum_df$Cell_line)
  col_map <- setNames(c("#4575B4", "#D73027")[seq_along(cell_levels)], cell_levels)
  
  for (cl in cell_levels) {
    d <- subset(sum_df, Cell_line == cl)
    if (nrow(d) == 0) next
    # Plot points at all doses, including 0
    points(d$Dose, d$mean, pch = 16, cex = 1.1, col = col_map[[cl]])
    # Plot SEM at all doses, including 0
    arrows(x0 = d$Dose, y0 = d$mean - d$sem,
           x1 = d$Dose, y1 = d$mean + d$sem,
           angle = 90, code = 3, length = 0.03,
           lwd = 1, col = col_map[[cl]])
  }
}

# Loop over files
for (f in files) {
  message("Processing file: ", f)
  
  df <- read.csv(file.path(data_dir, f), stringsAsFactors = FALSE)
  
  fit2 <- drm(Response ~ Dose,
              curveid = Cell_line,
              data    = df,
              fct     = LL.4())
  
  sig_df <- df %>%
    dplyr::group_by(Dose) %>%
    dplyr::summarize(p = t.test(Response ~ Cell_line)$p.value, .groups = "drop") %>%
    dplyr::mutate(
      p_adj = p.adjust(p, method = "holm"),
      sig = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
  
  p_main <- ggplotify::as.ggplot(function() draw_main_plot(df, fit2, sub("\\.csv$", "", f)))
  tbl <- gridExtra::tableGrob(sig_df, rows = NULL)
  
  # Conditional naming
  if (grepl("250714_MTT_data_BN91", f)) {
    pdf_name <- file.path(results_dir, "Fig6e_left_panel.pdf")
  } else if (grepl("250618_MTT_data_BN118", f)) {
    pdf_name <- file.path(results_dir, "Fig6e_right_panel.pdf")
  } else {
    base <- sub("\\.csv$", "", f)
    pdf_name <- file.path(results_dir, paste0("Suppl_Fig_MTT_", base, ".pdf"))
  }
  
  pdf(pdf_name, width = 6, height = 8)
  gridExtra::grid.arrange(p_main, tbl, heights = c(2/3, 1/3))
  dev.off()
}
