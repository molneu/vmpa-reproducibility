library(echarts4r)
library(dplyr)
library(htmlwidgets)
library(here)

here::i_am("reproducibility/supplementary_figure6/scripts/supplementary_figure6_function_annotation.R")

data_path <- here("reproducibility", "figure3", "data", "annotation_protein function.csv")
results_dir <- here("reproducibility", "supplementary_figure6", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

plot_data <- read.csv(data_path) %>%
  dplyr::select(Function, Collectri)

custom_colors <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
  "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"
)

plot <- plot_data %>%
  e_charts(Function, renderer = "svg") %>%
  e_pie(Collectri, radius = c("30%", "50%")) %>%
  e_labels(show = TRUE, position = "outside", formatter = "{b}", fontSize = 20) %>%
  e_legend(show = FALSE) %>%
  e_color(custom_colors) %>%
  e_tooltip(trigger = "item") %>%
  e_show_loading()

htmlwidgets::saveWidget(
  widget = plot,
  file = file.path(results_dir, "Supplementary_Figure_6_function_annotation.html"),
  selfcontained = TRUE
)
