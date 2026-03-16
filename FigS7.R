#code for Fig.S7
# Seurat	5.1.0
# SeuratData	0.2.2.9001
# ggplot2	3.4.4
#The code for slingshot is the same as that demonstrated in FIG1.R. 
#The heatmap is generated using the standard plotting method of pheatmap, based on scaled matrix data. 
#The spatial signaling pathway scoring plot employs the default parameters of the AddModuleScore function (https://satijalab.org/seurat/reference/addmodulescore).
library(Seurat)
library(ggplot2)
library(plotly)

# Load and examine data
seurat_obj <- readRDS('./sc_qc.rds')
dim(seurat_obj@assays$RNA@counts)
dim(seurat_obj@meta.data)

# Function to integrate expression data into meta.data
sl_genes_counts <- function(seurat_obj, features) {
  if (!is.character(features) || length(features) < 1) {
    stop("features must be a character vector of gene names")
  }
  
  counts_matrix <- Seurat::FetchData(seurat_obj, vars = features)
  seurat_obj@meta.data$end <- NA
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, counts_matrix)
  return(seurat_obj)
}

# Function to categorize expression data into three levels
color_split <- function(seurat_obj, features) {
  if (!is.character(features) || length(features) < 1) {
    stop("features must be a character vector of gene names")
  }
  
  for (gene in features) {
    if (!gene %in% names(seurat_obj@meta.data)) {
      stop(paste("Gene", gene, "not found in seurat_obj@meta.data"))
    }
    
    gene_min <- min(seurat_obj@meta.data[[gene]], na.rm = TRUE)
    gene_max <- max(seurat_obj@meta.data[[gene]], na.rm = TRUE)
    breakpoints <- seq(from = gene_min, to = gene_max, length.out = 4)
    
    seurat_obj@meta.data[[paste(gene, "level", sep = "_")]] <- cut(
      seurat_obj@meta.data[[gene]],
      breaks = breakpoints,
      labels = c("low", "medium", "high"),
      include.lowest = TRUE
    )
  }
  return(seurat_obj)
}

# Function for spatial feature plotting
spatial_feature_plot <- function(seurat_obj, features = c("CDH1"), mod = "expression", 
                                  save = "./test", show = FALSE, axis = TRUE, 
                                  min_size = 1.5, max_size = 3) {
  if (!file.exists(save)) {
    dir.create(save, recursive = TRUE)
  }
  
  seurat_obj <- sl_genes_counts(seurat_obj, features = features)
  seurat_obj <- color_split(seurat_obj, features = features)
  seurat_meta <- seurat_obj@meta.data
  
  library(viridis)
  print(paste("Drawing", features))
  
  if (mod == "expression") {
    i <- features
    col1 <- c("#bdbdbd", "#d9d9d9", "#fee5d9", "#fb6a4a", "#a50f15")
    seurat_meta$opacity <- ifelse(seurat_meta[[i]] == 0, 0.1, 1)
    
    if (is.numeric(seurat_meta[[i]])) {
      max_expr <- max(seurat_meta[[i]])
      seurat_meta$size <- min_size + (seurat_meta[[i]] / max_expr) * (max_size - min_size)
      seurat_meta$size[seurat_meta[[i]] == 0] <- min_size
    } else {
      seurat_meta$size <- min_size
    }
    
    p <- plot_ly(data = seurat_meta, x = ~x_coord, y = ~y_coord, 
                  z = ~slice_num, type = "scatter3d", mode = "markers", 
                  marker = list(size = ~size, opacity = ~opacity, line = list(width = 0)), 
                  color = ~get(i), colors = col1)
    
  } else if (mod == "levels") {
    i <- paste0(features, "_level")
    col1 <- unique(c("grey", "gold", "red"))
    seurat_meta$opacity <- 1
    seurat_meta$size <- max_size
    
    p <- plot_ly(data = seurat_meta, x = ~x_coord, y = ~y_coord, 
                  z = ~slice_num, type = "scatter3d", mode = "markers", 
                  marker = list(size = max_size, opacity = 1, line = list(width = 0)), 
                  color = ~get(i), colors = col1)
  }
  
  p <- p %>% layout(scene = list(
    aspectratio = list(x = 1, y = 1, z = 0.15), 
    bgcolor = "black", 
    xaxis = list(showgrid = axis, zeroline = axis, showticklabels = axis, 
                 title = if (axis) "x_coord" else ""), 
    yaxis = list(showgrid = axis, zeroline = axis, showticklabels = axis, 
                 title = if (axis) "y_coord" else ""), 
    zaxis = list(showgrid = axis, zeroline = axis, showticklabels = axis, 
                 title = if (axis) "slice_num" else "")
  ))
  
  if (!axis) {
    p <- p %>% layout(scene = list(
      xaxis = list(visible = FALSE), 
      yaxis = list(visible = FALSE), 
      zaxis = list(visible = FALSE)
    ))
  }
  
  if (show) {
    print(p)
  }
  
  htmlwidgets::saveWidget(p, file.path(save, paste0(features, "_", mod, ".html")))
}