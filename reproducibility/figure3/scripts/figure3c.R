library(echarts4r)
library(dplyr)
library(htmlwidgets)
library(webshot2)
library(here)

here::i_am("reproducibility/figure3/scripts/figure3c.R")

# ====== STEP 1: Setup directories and paths ======

# Define output directories for figure 3 reproducibility
fig3_data_dir    <- here("reproducibility", "figure3", "data")
fig3_results_dir <- here("reproducibility", "figure3", "results")

# Create directories if they do not exist (recursive, suppress warnings)
dir.create(fig3_data_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(fig3_results_dir, recursive = TRUE, showWarnings = FALSE)

# Paths
df_path <- file.path(fig3_data_dir, "annotation_protein function.csv")
df2_path <- file.path(fig3_data_dir, "annotation_protein localization.csv")

# ====== STEP 2: Read files and plot ======

df <- read.csv(df_path)
df2 <- read.csv(df2_path)

# Prepare data for plotting
plot_data <- df %>%
  dplyr::select(Function, Compass)

plot_data3 <- df2 %>%
  dplyr::select(Subcellular.localization, Compass)

#colors

#used for localization
custom_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", 
                   "#e6ab02", "#a6761d", "#666666", "#d53e4f", "#3288bd")

#used for function
custom_colors2 <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", 
                   "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")

# Create the donut plot with custom colors and no legend

#figure 3c top (function)
p1 <- plot_data %>%
  e_charts(Function, renderer = "svg") %>%
  e_pie(Compass, radius = c("30%", "50%")) %>%  # Smaller hole by decreasing inner radius
  e_labels(show = TRUE, position = 'outside', formatter = "{b}", fontSize = 20) %>%  # Show only the function names
  e_legend(show = FALSE) %>%  # Remove the legend
  e_color(custom_colors2) %>%  # Apply custom colors
  e_tooltip(trigger = 'item') %>%
  e_show_loading()

#figure 3c bottom (localization)
p3 <- plot_data3 %>%
  e_charts(Subcellular.localization, renderer = "svg") %>%
  e_pie(Compass, radius = c("30%", "50%")) %>%  # Smaller hole by decreasing inner radius
  e_labels(show = TRUE, position = 'outside', formatter = "{b}", fontSize = 20) %>%  # Show only the function names
  e_legend(show = FALSE) %>%  # Remove the legend
  e_color(custom_colors) %>%  # Apply custom colors
  e_tooltip(trigger = 'item') %>%
  e_show_loading()

# ====== STEP 3: Save ======

html_file1 <- file.path(fig3_results_dir, "Fig 3c_donut_function.html")
html_file3 <- file.path(fig3_results_dir, "Fig 3c_donut_localization.html")

htmlwidgets::saveWidget(
  widget      = p1,
  file        = html_file1,
  selfcontained = TRUE   # embed dependencies into one file
)

htmlwidgets::saveWidget(
  widget      = p3,
  file        = html_file3,
  selfcontained = TRUE   # embed dependencies into one file
)
