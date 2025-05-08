#code for Fig.S7

### cytotrace2

library(Seurat)
#library(SeuratData)
library(ggplot2)
library(SeuratDisk)
packageVersion("Seurat")
packageVersion("Matrix")
seurat_g1 <- LoadH5Seurat('./data/g1.h5seurat')
library(CytoTRACE2) #loading
seurat_g1[['RNA']]<- seurat_g1[['Spatial']]
DefaultAssay(seurat_g1) <- 'RNA'
cytotrace2_result <- cytotrace2(seurat_g1,is_seurat = TRUE,species = 'human')
FeaturePlot(cytotrace2_result,features = 'CytoTRACE2_Score')
# library(SCP)
# FeatureDimPlot(
#   srt = cytotrace2_result, features = c("CytoTRACE2_Score"),
#   reduction = "UMAP", theme_use = "theme_blank"
# )
# CellDimPlot(
#   srt = cytotrace2_result, group.by = "celltype", stat.by = "CytoTRACE2_Potency",
#   reduction = "UMAP", theme_use = "theme_blank"
# )

annotation <- data.frame(phenotype = cytotrace2_result@meta.data$celltype) %>% set_rownames(., colnames(cytotrace2_result))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation, 
                  is_seurat = TRUE)

library(ggpirate)
library(patchwork)
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
              axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels vertically
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank()) + # Remove y-axis title
        labs(y = NULL) + # Remove default y-axis title
        annotate("text", x = Inf, y = Inf, label = col, hjust = 1.1, vjust = 2, size = 5, angle = 0) # Add y-axis title at the top
    } else if (mod == "B") {
      p <- ggplot(data4plot, aes_string(x = "group", y = col, fill = "group")) +
        geom_pirate(aes_string(x = "group", y = col, fill = "group"), alpha = 0.4) +
        scale_fill_manual(values = color) + # Fill colors
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # Remove background grid lines
              axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels vertically
              axis.title.x = element_blank(), # Remove x-axis title
              axis.title.y = element_blank()) + # Remove y-axis title
        labs(y = NULL) + # Remove default y-axis title
        annotate("text", x = Inf, y = Inf, label = col, hjust = 1.1, vjust = 2, size = 5, angle = 0) # Add y-axis title at the top
    }
    return(p)
  }
  
  # Create plot list
  plot_list <- lapply(cols, plot_function)
  
  # Combine plots
  combined_plot <- wrap_plots(plot_list, ncol = length(cols))
  
  return(combined_plot)
}
mycols <- readRDS('./data/fig1_cols_617.rds')
vplot <- vln_plot_beauty(meta_data=cytotrace2_result@meta.data, cols = c('CytoTRACE2_Score'), color = mycols, group.by = 'order_617',mod = 'A')


### score

seurat_g1 <- AddModuleScore(
    object = seurat_g1,
    features = genes,
    ctrl = 100, #默认值是100
    name = 'Stemness_Score'
)
vplot <- vln_plot_beauty(meta_data=seurat_am@meta.data, cols = c('Stemness_Score1'), color = mycols, group.by = 'celltype',mod = 'A')
ggsave('./plot/Stemness_score.pdf',vplot,width = 5.5, height = 5)
```

### Inferring the relationship between space and pseudotime reveals the order of differentiation and migration

``` R
library(Seurat)
library(Seurat)
library(ggplot2)
library(dplyr)
library(MASS)
library(hexbin)
library(ggplot2)
library(viridis)
library(monocle3)
seurat_obj <- readRDS('./data/seurat_g1_slingshot_629.rds')
unique(seurat_obj$celltype)
seurat_obj <- subset(seurat_obj,celltype %in%c('Inter.Epi','Anterior.Epi','PrCP','Posterior.Epi') )
seurat_obj$z_coord <- seurat_obj$z_coord*100
#seurat_obj <- subset(seurat_obj,slice_num%in%c(1:13))
seurat_end <- subset(seurat_obj,celltype%in%c('PrCP'))
get_spatial_center <- function(seurat_obj, spatial_x = NULL, spatial_y = NULL, spatial_z = NULL, fake_point = FALSE) {
  # Get coordinate data
  message('Getting the coordinate center for end cluster...')
  coords <- list(x = spatial_x, y = spatial_y, z = spatial_z)
  coord_data <- lapply(coords, function(coord) {
    if (!is.null(coord)) {
      FetchData(seurat_obj, vars = coord)
    } else {
      NULL
    }
  })
  
  # Remove NULL values
  coord_data <- coord_data[!sapply(coord_data, is.null)]
  
  # Merge coordinate data
  coord_df <- do.call(cbind, coord_data)
  
  # Get barcodes
  barcodes <- rownames(coord_df)
  
  # Calculate center point
  center <- colMeans(coord_df)
  
  # Create result dataframe with column names
  result_df <- as.data.frame(t(center))
  colnames(result_df) <- colnames(coord_df)
  
  if (fake_point) {
    # If fake_point is TRUE, directly return center point coordinates
    return(list(
      coordinates = result_df,
      barcode = "center_point",
      is_fake = TRUE
    ))
  } else {
    # Calculate Euclidean distance from each cell to center
    distances <- sqrt(rowSums((coord_df - center)^2))
    
    # Find closest cell
    closest_cell_index <- which.min(distances)
    closest_cell_barcode <- barcodes[closest_cell_index]
    closest_cell_coords <- coord_df[closest_cell_index, , drop = FALSE]
    
    # Return coordinates and barcode of closest cell
    return(list(
      coordinates = closest_cell_coords,
      barcode = closest_cell_barcode,
      is_fake = FALSE
    ))
  }
}
# Using virtual center point (fake_point = TRUE)
result_fake <- get_spatial_center(seurat_end, spatial_x = "x_coord", spatial_y = "y_coord", spatial_z = "z_coord",fake_point = TRUE)
# Print result
print(result_fake)

get_root_coord <- function(seurat_obj, spatial_x = NULL, spatial_y = NULL, spatial_z = NULL, end.x, end.y, end.z = NULL, distance = 0.95, mod = 'A') {
  if (mod == 'A') {
      message('mod A: Getting coordinates of the cell with minimum pseudotime...')
    # Find cell with minimum pseudotime
    min_pseudotime_index <- which.min(seurat_obj$pseudotime)
    min_pseudotime_barcode <- rownames(seurat_obj@meta.data)[min_pseudotime_index]
    
    # Get coordinates of that cell
    coords <- list(x = spatial_x, y = spatial_y, z = spatial_z)
    coord_data <- lapply(coords, function(coord) {
      if (!is.null(coord)) {
        FetchData(seurat_obj, vars = coord)
      } else {
        NULL
      }
    })
    coord_data <- coord_data[!sapply(coord_data, is.null)]
    coord_df <- do.call(cbind, coord_data)
    min_pseudotime_coords <- coord_df[min_pseudotime_barcode, , drop = FALSE]
    
    return(list(
      coordinates = min_pseudotime_coords,
      barcode = min_pseudotime_barcode,
      distance = NA  # Not calculating distance here, so return NA
    ))
  } else if (mod == 'B') {
    message('mod B: Calculating starting coordinates for root cell...')
    # Original calculation logic, compatible with 2D and 3D
    coords <- list(x = spatial_x, y = spatial_y, z = spatial_z)
    coord_data <- lapply(coords, function(coord) {
      if (!is.null(coord)) {
        FetchData(seurat_obj, vars = coord)
      } else {
        NULL
      }
    })
    
    coord_data <- coord_data[!sapply(coord_data, is.null)]
    coord_df <- do.call(cbind, coord_data)
    barcodes <- rownames(coord_df)
    
    # Determine end_coords based on input dimensions
    if (is.null(end.z)) {
      end_coords <- c(end.x, end.y)
    } else {
      end_coords <- c(end.x, end.y, end.z)
    }
    
    # Ensure coord_df and end_coords dimensions match
    if (ncol(coord_df) != length(end_coords)) {
      stop("Dimension mismatch between input coordinates and end coordinates")
    }
    
    distances <- sqrt(rowSums((coord_df - end_coords)^2))
    
    sorted_indices <- order(distances, decreasing = TRUE)
    percentile_index <- ceiling(length(distances) * (1-distance))
    
    chosen_index <- sorted_indices[percentile_index]
    chosen_barcode <- barcodes[chosen_index]
    chosen_coords <- coord_df[chosen_index, , drop = FALSE]
    
    return(list(
      coordinates = chosen_coords,
      barcode = chosen_barcode,
      distance = distances[chosen_index]
    ))
  } else {
    stop("Invalid mod parameter. Use 'A' or 'B'.")
  }
}

# Example usage
seurat_root <- subset(seurat_obj, celltype %in% c('Inter.Epi'))

# # Using get_root_coord function, mod='A'
# root_result_A <- get_root_coord(seurat_root, spatial_x = "x_coord", spatial_y = "y_coord", spatial_z = "z_coord", 
#                                 end.x = result_fake$coordinates$x_coord, 
#                                 end.y = result_fake$coordinates$y_coord, 
#                                 end.z = result_fake$coordinates$z_coord, 
#                                 distance = 0.95, mod = 'A')
# print(root_result_A)

# # Using get_root_coord function, mod='B', 3D case
root_result_B_3D <- get_root_coord(seurat_root, spatial_x = "x_coord", spatial_y = "y_coord", spatial_z = "z_coord", 
                                   end.x = result_fake$coordinates$x_coord, 
                                   end.y = result_fake$coordinates$y_coord, 
                                   end.z = result_fake$coordinates$z_coord, 
                                   distance = 0.95, mod = 'B')
print(root_result_B_3D)

# # Using get_root_coord function, mod='B', 2D case
root_result_B_2D <- get_root_coord(seurat_root, spatial_x = "x_coord", spatial_y = "y_coord", 
                                   end.x = result_fake$coordinates$x_coord, 
                                   end.y = result_fake$coordinates$y_coord, 
                                   distance = 0.7, mod = 'B')
print(root_result_B_2D)

root_result = root_result_B_2D

get_spatial_distance <- function(seurat_obj, spatial_x = NULL, spatial_y = NULL, spatial_z = NULL, root.x = NULL, root.y = NULL, root.z = NULL) {
    message('Calculating the distance from each point to the starting point and storing it in the Seurat object...')
  # Get coordinate data
  coords <- list(x = spatial_x, y = spatial_y, z = spatial_z)
  coord_data <- lapply(coords, function(coord) {
    if (!is.null(coord)) {
      FetchData(seurat_obj, vars = coord)
    } else {
      NULL
    }
  })
  
  # Remove NULL values
  coord_data <- coord_data[!sapply(coord_data, is.null)]
  
  # Merge coordinate data
  coord_df <- do.call(cbind, coord_data)
  
  # Create root coordinates vector
  root_coords <- c(root.x, root.y, root.z)
  root_coords <- root_coords[!sapply(root_coords, is.null)]
  
  # Check if number of coordinates matches
  if (length(root_coords) != ncol(coord_df)) {
    stop("The number of coordinates provided for spatial and root points must match.")
  }
  
  # Calculate Euclidean distance from each cell to root
  distances <- sqrt(rowSums((coord_df - root_coords)^2))
  
  # Store distances in Seurat object metadata
  seurat_obj@meta.data$spatial_distance <- distances
  
  # Return updated Seurat object
  return(seurat_obj)
}
# Using center point coordinates as starting point
seurat_obj <- get_spatial_distance(seurat_obj, spatial_x = "x_coord", spatial_y = "y_coord",
                                   root.x = root_result$coordinates["x_coord"], 
                                   root.y = root_result$coordinates["y_coord"]
                                   )#,root.z = root_result$coordinates["z_coord"], spatial_z = "z_coord"
DimPlot(seurat_obj,cells.highlight = root_result$barcode)

########Insert monocle3 function here, dependent on the barcode information obtained above
#root_result$barcode
call_monocle3 <- function(seurat_obj, group.by, root.cells, reduction = 'umap', assay = 'Spatial') {
  library(monocle3)
  library(tidyverse)
  library(patchwork)
  library(viridis)
  
  message("Starting Monocle3 analysis...")
  
  # Extract data from Seurat object
  message("Extracting data from Seurat object...")
  data <- GetAssayData(seurat_obj, assay = assay, slot = 'counts')
  cell_metadata <- seurat_obj@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  
  # Create cell_data_set object
  message("Creating cell_data_set object...")
  cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
  
  # Preprocess and reduce dimension
  message("Preprocessing and reducing dimensions...")
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, preprocess_method = "PCA")
  
  # Use UMAP from Seurat object
  message("Integrating UMAP from Seurat object...")
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(seurat_obj, reduction = reduction)
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  
  # Cluster cells and learn graph
  message("Clustering cells and learning graph...")
  cds <- cluster_cells(cds, reduction_method = 'UMAP', cluster_method = 'louvain', k = 350)
  cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
  
  # Plot cells
  message("Generating cell plot...")
  options(repr.plot.height = 7, repr.plot.width = 8)
  p1 <- plot_cells(cds, color_cells_by = group.by, label_groups_by_cluster = FALSE,group_label_size = 8,
                   graph_label_size = 4, cell_size = 1, trajectory_graph_segment_size = 1.2,               label_branch_points = FALSE, # Don't show branch point labels
                label_roots = FALSE, # Don't show root node labels
                label_leaves = FALSE ,trajectory_graph_color = 'black' ) #,trajectory_graph_color = '#bababa' 
  print(p1)
  #p1 = p1+scale_color_manual(values = mycols)
  # Order cells
  message("Ordering cells...")
  cds <- order_cells(cds, root_cells = root.cells)
  
  # Plot pseudotime
  message("Generating pseudotime plot...")
  mycol <- inferno(10)
  p2 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE, alpha = 0.7,
                   cell_stroke = 0.1, group_label_size = 8, graph_label_size = FALSE,
                   cell_size = 2.2, trajectory_graph_segment_size = 0.5) +
    scale_color_gradient2(low = "#8C3B77", mid = "#FF614D", high = "#FFEC5C", 
                          midpoint = median(cds@principal_graph_aux@listData$UMAP$pseudotime))
  print(p2)
  
  # Get pseudotime data
  message("Calculating pseudotime...")
  pseudotime_data <- pseudotime(cds, reduction_method = "UMAP")
  pseudotime_df <- data.frame(
    barcode = colnames(cds),
    pseudotime = pseudotime_data
  )
  
  # Add pseudotime to Seurat object
  message("Adding pseudotime to Seurat object...")
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = pseudotime_df$pseudotime,
    col.name = "pseudotime"
  )
  
  message("Monocle3 analysis complete!")
  return(seurat_obj)
}
seurat_obj <- call_monocle3(seurat_obj, 
                                            group.by = "celltype", 
                                            root.cells = root_result$barcode, 
                                            reduction = 'umap', 
                                            assay = 'Spatial')  # Can change assay if needed
########Remove outliers
remove_outliers <- function(seurat_obj, group.by) {
  data <- seurat_obj@meta.data
  message('About to remove outliers using the interquartile range method...')
  # Debug info: print original data row count
  message("Original number of cells:", nrow(data), "\n")
  
  data_filtered <- data %>%
    group_by(!!sym(group.by)) %>%
    mutate(
      is_outlier = scaled_spatial_distance < quantile(scaled_spatial_distance, 0.25) - 1.5 * IQR(scaled_spatial_distance) |
                   scaled_spatial_distance > quantile(scaled_spatial_distance, 0.75) + 1.5 * IQR(scaled_spatial_distance) |
                   scaled_pseudotime < quantile(scaled_pseudotime, 0.25) - 1.5 * IQR(scaled_pseudotime) |
                   scaled_pseudotime > quantile(scaled_pseudotime, 0.75) + 1.5 * IQR(scaled_pseudotime)
    ) %>%
    ungroup()
  
  # Get non-outlier barcodes
  kept_barcodes <- data_filtered %>%
    filter(!is_outlier) %>%
    pull(barcode)
  
  # Debug info: print filtered data row count
  message("Number of cells after outlier removal:", length(kept_barcodes), "\n")
  
  # Check if there are remaining cells
  if (length(kept_barcodes) == 0) {
    warning("No cells remain after removing outliers. Returning original object.")
    return(seurat_obj)
  }
  
  # Return filtered Seurat object
  seurat_obj_filtered <- subset(seurat_obj, cells = kept_barcodes)
  return(seurat_obj_filtered)
}
scale_to_01 <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

spatial_pseudo_scatter <- function(seurat_obj, group.by, root, end, cols = NULL, hex_plot = TRUE, hex_bins = 50, viridis_color = 'mako') {
  # Ensure necessary packages are loaded
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    install.packages("hexbin")
  }

  
  # Ensure necessary columns exist
  required_cols <- c("scaled_spatial_distance", "scaled_pseudotime", group.by)
  if (!all(required_cols %in% colnames(seurat_obj@meta.data))) {
    stop("One or more required columns are missing from the Seurat object metadata.")
  }
  seurat_obj$celltype <- as.character(seurat_obj@meta.data[[group.by]])
  
  # Subset Seurat object
  seurat_sub <- subset(seurat_obj, celltype %in% c(root, end))
  
  # Extract data
  data <- seurat_sub@meta.data[, c("scaled_spatial_distance", "scaled_pseudotime", group.by)]
  
  # Calculate density
  kde <- kde2d(data$scaled_spatial_distance, data$scaled_pseudotime, n = 100)
  density_df <- expand.grid(x = kde$x, y = kde$y)
  density_df$density <- as.vector(kde$z)
  
  # Assign density value to each point
  data$density <- NA
  for (i in 1:nrow(data)) {
    x <- data$scaled_spatial_distance[i]
    y <- data$scaled_pseudotime[i]
    closest_point <- which.min((density_df$x - x)^2 + (density_df$y - y)^2)
    data$density[i] <- density_df$density[closest_point]
  }
  
  # Add density values to Seurat object
  seurat_sub@meta.data$density <- NA
  seurat_sub@meta.data[rownames(data), "density"] <- data$density
  
  # Create base plot
  base_plot <- ggplot(data, aes(x = scaled_spatial_distance, y = scaled_pseudotime)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          legend.key = element_blank()) +
    labs(x = "Scaled Spatial Distance", y = "Scaled Pseudotime") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "dark red", alpha = 0.7)
  
  # Draw density scatter plot (p1)
  p1 <- base_plot +
    geom_point(aes(color = density)) +
    scale_color_viridis_c(option = viridis_color) +
    guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, label.position = "right", title = "Density", labels = c("Low", "High")))
  
  print(p1)
  
  # If hex_plot is TRUE, draw Hex plot
  if (hex_plot) {
    message('hex_plot = TRUE, about to draw hex density plot')
    p_hex <- base_plot +
      geom_hex(bins = hex_bins) +
      #scale_fill_viridis_c(option = viridis_color) +
      scale_fill_gradientn(colors = fangao_color) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, label.position = "right", title = "Density", labels = c("Low", "High")))
    
    print(p_hex)
  }
  
  # Draw scatter plot grouped by group.by (p2)
  if (!is.null(cols)) {
    unique_groups <- unique(data[[group.by]])
    if (length(cols) != length(unique_groups)) {
      stop("The number of colors provided does not match the number of groups.")
    }
    
    p2 <- base_plot +
      geom_point(aes(color = !!sym(group.by))) +
      scale_color_manual(values = cols[unique_groups]) +
      labs(color = group.by)
    
    print(p2)
  }
  # Return list of plots
  plot_list <- list(p1 = p1, p_hex = p_hex, p2 = p2)
  
  return(plot_list)
}

# Read color file
mycols <- readRDS('./data/fig1_cols_617.rds')

# First normalization to 0-1 range
seurat_obj@meta.data$scaled_spatial_distance <- scale_to_01(seurat_obj@meta.data$spatial_distance)
seurat_obj@meta.data$scaled_pseudotime <- scale_to_01(seurat_obj@meta.data$pseudotime)

# Remove outliers
seurat_obj <- remove_outliers(seurat_obj, group.by = "celltype")

# Second normalization to 0-1 range
seurat_obj@meta.data$scaled_spatial_distance <- scale_to_01(seurat_obj@meta.data$spatial_distance)
seurat_obj@meta.data$scaled_pseudotime <- scale_to_01(seurat_obj@meta.data$pseudotime)

# # Ensure mycols contains required colors
# if (!all(c('Inter.Epi', 'AM') %in% names(mycols))) {
#   stop("The provided color vector does not contain the required groups.")
# }

# Use spatial_pseudo_scatter function
result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = "celltype", root = "Inter.Epi", end = "AM", cols = mycols[c('Inter.Epi', 'AM')])
# # Use spatial_pseudo_scatter function
# result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = "celltype", root = "Inter.Epi", end = "Anterior.Epi", cols = mycols[c('Inter.Epi', 'Anterior.Epi')])
# result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = "celltype", root = "Inter.Epi", end = c("AM","Anterior.Epi"), cols = mycols[c("AM","Anterior.Epi","Inter.Epi")])

# # Use spatial_pseudo_scatter function
# result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = "celltype", root = "Inter.Epi", end = "PrCP", cols = mycols[c('Inter.Epi', "PrCP")])

# result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = "celltype", root = "Inter.Epi", end = "Anterior.Epi", cols = mycols[c('Inter.Epi', "Anterior.Epi")])

# result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = "celltype", root = "Inter.Epi", end = c("PrCP","Anterior.Epi"), cols = mycols[c('Inter.Epi', "Anterior.Epi","PrCP")])

Spatial_move <- function(seurat_obj, group.by, root_cluster, end_cluster, 
                         spatial_x = "x_coord", spatial_y = "y_coord", spatial_z = NULL, 
                         fake_point = TRUE, time_col = 'pseudotime', cols = NULL, 
                         remove_out_point = TRUE, use_root = 'min_time',assay = 'Spatial') {
  
  # Step 1: Get the end point coordinates
  seurat_end <- subset(seurat_obj, celltype %in% end_cluster)
  result_end <- get_spatial_center(seurat_end, spatial_x = spatial_x, spatial_y = spatial_y, 
                                   spatial_z = spatial_z, fake_point = fake_point)
  
  # Step 2: Get the root coordinates
  seurat_root <- subset(seurat_obj, celltype %in% root_cluster)
  
  if (use_root == 'min_time') {
    root_result <- get_root_coord(seurat_root, spatial_x = spatial_x, spatial_y = spatial_y, 
                                  spatial_z = spatial_z, mod = 'A')
  } else if (use_root %in% c('2D', '3D')) {
    root_result <- get_root_coord(seurat_root, spatial_x = spatial_x, spatial_y = spatial_y, 
                                  spatial_z = if(use_root == '3D') spatial_z else NULL, 
                                  end.x = result_end$coordinates[[spatial_x]], 
                                  end.y = result_end$coordinates[[spatial_y]], 
                                  end.z = if(use_root == '3D') result_end$coordinates[[spatial_z]] else NULL, 
                                  distance = 0.75, mod = 'B')
  seurat_obj <- call_monocle3(seurat_obj, 
                                            group.by = group.by, 
                                            root.cells = root_result$barcode, 
                                            reduction = 'umap', 
                                            assay = assay)  # Can change assay if needed
  } else {
    stop("Invalid use_root parameter. Use 'min_time', '2D', or '3D'.")
  }
  
# Step 3: Calculate spatial distances
if (is.null(spatial_z)) {
  # 2D case
# Use center point coordinates as starting point
seurat_obj <- get_spatial_distance(seurat_obj, spatial_x = "x_coord", spatial_y = "y_coord", 
                                   root.x = root_result$coordinates["x_coord"], 
                                   root.y = root_result$coordinates["y_coord"]
                                   )
} else {
  # 3D case
  # Use center point coordinates as starting point
seurat_obj <- get_spatial_distance(seurat_obj, spatial_x = "x_coord", spatial_y = "y_coord", spatial_z = "z_coord",
                                   root.x = root_result$coordinates["x_coord"], 
                                   root.y = root_result$coordinates["y_coord"],root.z = root_result$coordinates["z_coord"]
                                   )
}
  
  # Step 4: Scale distances and pseudotime
  seurat_obj@meta.data$scaled_spatial_distance <- scale_to_01(seurat_obj@meta.data$spatial_distance)
  seurat_obj@meta.data$scaled_pseudotime <- scale_to_01(seurat_obj@meta.data[[time_col]])
  
  # Step 5: Remove outliers if specified
  if (remove_out_point) {
    seurat_obj <- remove_outliers(seurat_obj, group.by = group.by)
    # Re-scale after removing outliers
    seurat_obj@meta.data$scaled_spatial_distance <- scale_to_01(seurat_obj@meta.data$spatial_distance)
    seurat_obj@meta.data$scaled_pseudotime <- scale_to_01(seurat_obj@meta.data[[time_col]])
  }
  
  # Step 6: Generate scatter plots
  result_seurat <- spatial_pseudo_scatter(seurat_obj, group.by = group.by, 
                                          root = root_cluster, end = end_cluster, 
                                          cols = cols)
  message('Visualization complete. You can use the returned seurat object to execute spatial_pseudo_scatter function separately for further visualization refinement without rerunning this function.')
  return(list(seurat_obj = seurat_obj , root_coords = root_result$coordinates, 
              end_coords = result_end$coordinates))
}
                                  
# Read color file
mycols <- readRDS('./data/fig1_cols_617.rds')                             
seurat_obj <- readRDS('./data/seurat_g1_slingshot_629.rds')
unique(seurat_obj$celltype)
seurat_obj <- subset(seurat_obj,celltype %in%c('Inter.Epi','Anterior.Epi','PrCP','Posterior.Epi') )
seurat_obj$z_coord <- seurat_obj$z_coord*100
unique(seurat_obj$celltype)
                                  
result <- Spatial_move(seurat_obj, 
                       group.by = "celltype", 
                       root_cluster = "Inter.Epi", 
                       end_cluster = c("PrCP","Anterior.Epi"), #, "Anterior.Epi"
                       spatial_x = "x_coord", 
                       spatial_y = "y_coord", 
                       #spatial_z = "z_coord", #If using 2D root, don't input z
                       fake_point = TRUE,
                       time_col = 'pseudotime',
                       cols = mycols[c( "Inter.Epi","PrCP","Anterior.Epi")], #, "Anterior.Epi"
                       remove_out_point = TRUE,
                       use_root = '2D')

pl1 <- spatial_pseudo_scatter(result$seurat_obj, group.by = "celltype", root = "Inter.Epi", end = "Posterior.Epi", cols = mycols[c('Inter.Epi', 'Posterior.Epi')])
pl2 <- spatial_pseudo_scatter(result$seurat_obj, group.by = "celltype", root = "Inter.Epi", end = c("Anterior.Epi",'PrCP'), cols = mycols[c('Inter.Epi', "Anterior.Epi",'PrCP')])
