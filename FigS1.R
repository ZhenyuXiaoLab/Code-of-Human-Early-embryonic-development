# code for Fig.S1

library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare data
plot_data <- data.frame(
  slice_num = factor(seurat_obj$slice_num),
  nCount_Spatial = seurat_obj$nCount_Spatial,
  nFeature_Spatial = seurat_obj$nFeature_Spatial,
  percent.mito = seurat_obj$percent.mito
)

plot_data_long <- gather(plot_data, key = "Metric", value = "Value", -slice_num)

# Calculate mean and standard error for each group
plot_data_summary <- plot_data_long %>%
  group_by(slice_num, Metric) %>%
  summarise(
    mean_value = mean(Value),
    se = sd(Value)/sqrt(n())
  )

# Calculate maximum value for left y-axis
left_y_max <- max(plot_data_summary$mean_value[plot_data_summary$Metric != "percent.mito"]) * 1.1

# Create bar plot
p <- ggplot() +
  geom_col(data = plot_data_summary %>% filter(Metric != "percent.mito"), 
           aes(x = slice_num, y = mean_value, fill = Metric),
           position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_errorbar(data = plot_data_summary %>% filter(Metric != "percent.mito"),
                aes(x = slice_num, ymin = mean_value - se, ymax = mean_value + se, group = Metric),
                position = position_dodge(0.8),
                width = 0.2) +
  geom_line(data = plot_data_summary %>% filter(Metric == "percent.mito"),
            aes(x = slice_num, y = mean_value * (left_y_max / 5.5), group = 1, color = "percent.mito"),
            size = 1) +
  geom_point(data = plot_data_summary %>% filter(Metric == "percent.mito"),
             aes(x = slice_num, y = mean_value * (left_y_max / 5.5), color = "percent.mito"),
             size = 3) +
  scale_fill_manual(values = c("#aecfd4","#edb8b0"), name = "Metric") +
  scale_color_manual(values = c("percent.mito" = "#1b7837"), name = "Metric") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  labs(y = "count/feature_num") +
  scale_y_continuous(
    limits = c(0, left_y_max),
    sec.axis = sec_axis(~ . * (5.5 / left_y_max), name = "percent.mito", breaks = seq(0, 5.5, by = 1))
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line.y.right = element_line(color = "#1b7837"),
    axis.ticks.y.right = element_line(color = "#1b7837"),
    axis.text.y.right = element_text(color = "#1b7837"),
    axis.title.y.right = element_text(color = "#1b7837")
  )

print(p)
ggsave('./plot/qc_by_slices_with_mito.pdf', p, width = 14, height = 5)

library(ggpirate)

vln_plot_beauty <- function(meta_data, cols, color, group.by, mod = "A") {
  # Check input parameters
  if (missing(meta_data)) stop("meta_data is required")
  if (missing(cols)) stop("cols is required")
  if (missing(color)) stop("color is required")
  if (missing(group.by)) stop("group.by is required")
  
  # Extract required data
  data4plot <- meta_data[, c(group.by, cols)]
  data4plot$group <- meta_data[[group.by]]
  
  # Define plotting function
  plot_function <- function(col) {
    if (mod == "A") {
      p <- ggplot(data4plot, aes_string(x = "group", y = col, fill = "group")) +
        geom_violin(alpha = 0.4) + # Violin plot needs some transparency
        stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.1), width = 0.1) + # Add error bars
        geom_boxplot(alpha = 0.5, outlier.size = 0, size = 0.3, width = 0.3) + # Add boxplot
        scale_fill_manual(values = color) + # Fill colors
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # Remove background grid lines
              axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank()) + # Remove y-axis title
        labs(y = NULL) + # Remove default y-axis title
        annotate("text", x = Inf, y = Inf, label = col, hjust = 1.1, vjust = 2, size = 5, angle = 0) # Add y-axis title at top
    } else if (mod == "B") {
      p <- ggplot(data4plot, aes_string(x = "group", y = col, fill = "group")) +
        geom_pirate(aes_string(x = "group", y = col, fill = "group"), alpha = 0.4) +
        scale_fill_manual(values = color) + # Fill colors
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # Remove background grid lines
              axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank()) + # Remove y-axis title
        labs(y = NULL) + # Remove default y-axis title
        annotate("text", x = Inf, y = Inf, label = col, hjust = 1.1, vjust = 2, size = 5, angle = 0) # Add y-axis title at top
    }
    return(p)
  }
  
  # Create plot list
  plot_list <- lapply(cols, plot_function)
  
  # Combine plots
  combined_plot <- wrap_plots(plot_list, ncol = length(cols))
  
  return(combined_plot)
}

p = vln_plot_beauty(meta_data = seurat_obj, cols = 'percent.mito', group.by = 'celltype', mod = "A")

