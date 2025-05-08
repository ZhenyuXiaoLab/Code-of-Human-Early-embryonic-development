#code for Figure1
# Data integration and clustering

> 2025-5-7

## Integration and Formatting of Chip Data

library(Seurat)
#library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)
color_feature <- mako(9999)
read_many_rds <- function(path = './', pattern = '.rds') {
  # Load required libraries
  library(tools)
  
  # Get all filenames matching the pattern in the specified path
  file_paths <- list.files(path = path, pattern = paste0("\\", pattern, "$"), full.names = TRUE)
  
  # Initialize a list to store read objects
  obj.list <- list()
  
  # Iterate through file paths and read each .rds file
  for (file_path in file_paths) {
    # Extract filename without extension from file path
    file_name <- basename(file_path)
    file_name <- sub(pattern = paste0(pattern, "$"), replacement = "", x = file_name)
    message(paste('Reading', file_path))
    obj.list[[file_name]] <- readRDS(file_path)
  }
  
  # Return the list containing all read objects
  return(obj.list)
}

obj.list2 <- read_many_rds('/home/xlyang/align_result/human_saw8_cs6/saw8_result')

# Get names of each object in the list using names function
obj_names <- names(obj.list2)

# Iterate through each object name
for (name in obj_names) {
  # Remove trailing underscore from object name
  clean_name <- sub("_$", "", name)
  
  # Assign value to orig.ident field in meta.data for each object
  obj.list2[[name]]@meta.data$orig.ident <- rep(clean_name)
}

# Create a vector mapping chip IDs to sample names
chip_to_sample <- c("C03650F1"="EV1-62", "C03647E2"="EV1-135", "C03650A5"="EV1-105",
                    "C03647E6"="EV1-31", "C03650E6"="EV1-116", "C03649A6"="EV1-128",
                    "A03398C2"="EV1-68", "A03398E3"="EV1-24")

# Iterate through each object in the list
for (name in names(obj.list2)) {
  # Get object name without underscore
  clean_name <- sub("_$", "", name)
  # Find corresponding sample name using mapping vector
  sample_name <- chip_to_sample[clean_name]
  
  # Assign value to sample column for each object
  obj.list2[[name]]@meta.data$sample <- sample_name
  
  # To avoid duplicate names, check if new name already exists in list
  if (sample_name %in% names(obj.list2)) {
    stop(paste("Error: the name", sample_name, "already exists in obj.list!"))
  }
  
  # Update name in list
  names(obj.list2)[names(obj.list2) == name] <- sample_name
}

# Iterate through each object in the list
for (sample_name in names(obj.list2)) {
  # Get current object's meta.data
  meta_data <- obj.list2[[sample_name]]@meta.data
  # Assign modified meta.data back to object
  obj.list2[[sample_name]]@meta.data$`_index` = gsub("sample:", "", obj.list2[[sample_name]]@meta.data$`_index`)
  obj.list2[[sample_name]]@meta.data$barcode <- paste(meta_data$sample, meta_data$'_index', sep = '_')
  obj.list2[[sample_name]] <- RenameCells(obj.list2[[sample_name]], new.names = paste(meta_data$sample, x = obj.list2[[sample_name]]$'_index', sep = '_'))
}

for (i in names(obj.list2)) {
  # Set default assay
  DefaultAssay(obj.list2[[i]]) <- 'Spatial'
}

objlistbak <- obj.list2

# Extract numeric part from object names
nums <- sapply(names(obj.list2), function(name) {
  as.numeric(gsub("EV1-", "", name))
})

# Sort based on extracted numbers
sorted_names <- names(obj.list2)[order(nums)]

# Reorder the list
obj.list2 <- obj.list2[sorted_names]

# Initialize empty matrix for merging cell.embeddings
all_embeddings <- matrix(nrow = 0, ncol = 2)

# Iterate through obj.list2 and merge cell.embeddings
for (sample_name in names(obj.list2)) {
    spatial_emb <- obj.list2[[sample_name]]@reductions$spatial@cell.embeddings
    colnames(spatial_emb) <- c('spatial_1', 'spatial_2')
    rownames(spatial_emb) <- obj.list2[[sample_name]]@meta.data$barcode
    all_embeddings <- rbind(all_embeddings, spatial_emb)
}

### Integration
seurat_obj <- obj.list2[[1]]
for(i in 2:length(obj.list2)) {
  seurat_obj <- merge(seurat_obj, y = obj.list2[[i]])
}

Idents(seurat_obj) <- seurat_obj@meta.data$sample
seurat_obj@meta.data$barcode = rownames(seurat_obj@meta.data)


### Try more batch correction methods and perform LISI scoring

library(Seurat)
#library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)
library(scplotter)
color_feature <- mako(9999)
library(ggsci)
seurat_obj <- LoadH5Seurat('./data/saw8_bin30_obj4.h5seurat')
mycols = readRDS('./data/CS6_color_2412.rds')
simpsons_colors <- pal_simpsons("springfield")(16)
sample_cols = setNames(simpsons_colors[1:length(unique(seurat_obj$sample))], unique(seurat_obj$sample))

p1 = CellDimPlot(seurat_obj, group_by = "celltype", reduction = "umap",
            label = FALSE,palcolor  = mycols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
ggsave('./plot/qupici_normal_p1.pdf',p1,width = 8,height = 7)
p2 = CellDimPlot(seurat_obj, group_by = "sample", reduction = "umap",
            label = FALSE,palcolor  = sample_cols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
ggsave('./plot/qupici_normal_p2.pdf',p2,width = 8,height = 7)
p3 = CellDimPlot(seurat_obj, group_by = "celltype", reduction = "umap",
            label = FALSE,palcolor  = c(sample_cols,mycols),pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE,split_by = 'sample')
ggsave('./plot/qupici_normal_p3.pdf',p3,width = 24,height = 21)

# seurat_obj <- SCTransform(seurat_obj,vars.to.regress = c('sample','percent.mito'),assay = 'Spatial')
# seurat_obj <- RunPCA(seurat_obj,npcs = 50,verbose = T)
# # Plot PCA
# #DimPlot(seurat_obj, reduction = "pca")
# VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
# ElbowPlot(seurat_obj, ndims = 50)
# #From the plot, the first 30 PCs are sufficient
# seurat_obj <- RunUMAP(seurat_obj,reduction = "pca",dims = 1:30,return.model = TRUE ) #，min.dist= ,spread=These parameters are the core of iterative optimization, you can write a loop function to adjust UMAP.
# #seurat_obj <- FindNeighbors(seurat_obj,reduction = "pca",dims = 1:10)#Default is 20
# # saveRDS(seurat_obj,'./data/SCT_obj5.rds')
# # seurat_obj <- readRDS('./data/SCT_obj5.rds')
# p1 = CellDimPlot(seurat_obj, group_by = "celltype", reduction = "umap",
#             label = FALSE,palcolor  = mycols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
# ggsave('./plot/qupici_sct_p1.pdf',p1,width = 8,height = 7)
# p2 = CellDimPlot(seurat_obj, group_by = "sample", reduction = "umap",
#             label = FALSE,palcolor  = sample_cols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
# ggsave('./plot/qupici_sct_p2.pdf',p2,width = 8,height = 7)
# p3 = CellDimPlot(seurat_obj, group_by = "celltype", reduction = "umap",
#             label = FALSE,palcolor  = c(sample_cols,mycols),pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE,split_by = 'sample')
# ggsave('./plot/qupici_sct_p3.pdf',p3,width = 24,height = 21)

######
# ##CCA integration
# DefaultAssay(seurat_obj) = 'Spatial'
# seurat_list <- SplitObject(seurat_obj, split.by = "sample")

# #Perform normalization and find highly variable genes separately
# seurat_list <- lapply(X = seurat_list, FUN = function(x) {
#     x <- NormalizeData(x)
#     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
# })
# # Find integration anchors (batch correction)
# features <- SelectIntegrationFeatures(object.list = seurat_list)
# immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
# #immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features,dims = 1:10,k.filter = 10,k.anchor = 10,k.score = 10)
# # this command creates an 'integrated' data assay
# immune.combined <- IntegrateData(anchorset = immune.anchors) #,k.weight = 50
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
# #immune.combined <- IntegrateData(anchorset = immune.anchors,features.to.integrate = rownames(immune.combined)) #,k.weight = 50
# DefaultAssay(immune.combined) <- "integrated" #The purpose of this step is unclear, it might be to only use common highly variable genes to speed up subsequent steps, but this would make quality control and volcano plots impossible
# #Normalization, PCA, UMAP plotting
# seurat_cca = immune.combined
# rm(immune.combined)
# seurat_cca <- ScaleData(seurat_cca, verbose = T)
# seurat_cca <- RunPCA(seurat_cca,npcs = 30,verbose = T,assay = 'integrated')
# # Plot PCA
# DimPlot(seurat_cca, reduction = "pca")
# VizDimLoadings(seurat_cca, dims = 1:2, reduction = "pca")
# ElbowPlot(seurat_cca, ndims = 30)
# #From the plot, the first 30 PCs are sufficient
# seurat_cca <- RunUMAP(seurat_cca,reduction = "pca",dims = 1:15 ) #，min.dist= 
# # seurat_cca <- seurat_harmony
# # rm(seurat_harmony)
# # rm(seurat_list)
# # gc()
# p1 = CellDimPlot(seurat_cca, group_by = "celltype", reduction = "umap",
#             label = FALSE,palcolor  = mycols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
# ggsave('./plot/qupici_cca_p1.pdf',p1,width = 8,height = 7)
# p2 = CellDimPlot(seurat_cca, group_by = "sample", reduction = "umap",
#             label = FALSE,palcolor  = sample_cols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
# ggsave('./plot/qupici_cca_p2.pdf',p2,width = 8,height = 7)
# p3 = CellDimPlot(seurat_cca, group_by = "celltype", reduction = "umap",
#             label = FALSE,palcolor  = c(sample_cols,mycols),pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE,split_by = 'sample')
# ggsave('./plot/qupici_cca_p3.pdf',p3,width = 24,height = 21)
##Try harmony integration
DefaultAssay(seurat_obj) = 'Spatial'
library(dplyr)
scRNA_harmony <- seurat_obj
# scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
library(harmony)
#Start integration
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "sample")})
#Dimensionality reduction and clustering
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15)# %>% FindClusters(resolution = 0.1)
#Dimensionality reduction visualization
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15)
p1 = CellDimPlot(scRNA_harmony, group_by = "celltype", reduction = "umap",
            label = FALSE,palcolor  = mycols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
ggsave('./plot/qupici_har_p1.pdf',p1,width = 8,height = 7)
p2 = CellDimPlot(scRNA_harmony, group_by = "sample", reduction = "umap",
            label = FALSE,palcolor  = sample_cols,pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE)
ggsave('./plot/qupici_har_p2.pdf',p2,width = 8,height = 7)
p3 = CellDimPlot(scRNA_harmony, group_by = "celltype", reduction = "umap",
            label = FALSE,palcolor  = c(sample_cols,mycols),pt.size = 0.01,theme = "theme_blank", legend.position = "right",raster = FALSE,split_by = 'sample')
ggsave('./plot/qupici_har_p3.pdf',p3,width = 24,height = 21)

#######Evaluate batch effect removal and biological variation preservation objectively through scoring
seurat_obj <- AddScoreLISI(seurat_obj, integration = "unintegrated",
                            batch.var = 'sample', cell.var = 'celltype',
                            reduction = "pca")
mean(seurat_obj@meta.data$iLISI_sample)
mean(seurat_obj@meta.data$cLISI_celltype)
seurat_cca <- AddScoreLISI(seurat_cca, integration = "integrated",
                            batch.var = 'sample', cell.var = 'celltype',
                            reduction = "pca")
mean(seurat_cca@meta.data$iLISI_sample)
mean(seurat_cca@meta.data$cLISI_celltype)
DefaultAssay(scRNA_harmony) <- 'Spatial'
scRNA_harmony <- AddScoreLISI(scRNA_harmony, integration = "integrated",
                            batch.var = 'sample', cell.var = 'celltype',
                            reduction = "harmony")
mean(scRNA_harmony@meta.data$iLISI_sample)
mean(scRNA_harmony@meta.data$cLISI_celltype)
DefaultAssay(scRNA_harmony) <- 'SCT'
scRNA_sct = RunPCA(scRNA_harmony)
scRNA_sct <- AddScoreLISI(scRNA_sct, integration = "integrated",
                            batch.var = 'sample', cell.var = 'celltype',
                            reduction = "pca")
mean(scRNA_sct@meta.data$iLISI_sample)
mean(scRNA_sct@meta.data$cLISI_celltype)

seurat_obj <- subset(seurat_obj,barcode%in%seurat_cca$barcode)

library(dplyr)
# Extract iLISI and cLISI from each object
df_harmony <- data.frame(barcode = scRNA_harmony$barcode,
                         iLISI_sample = scRNA_harmony$iLISI_sample,
                         cLISI_celltype = scRNA_harmony$cLISI_celltype)

df_cca <- data.frame(barcode = seurat_cca$barcode,
                     iLISI_sample = seurat_cca$iLISI_sample,
                     cLISI_celltype = seurat_cca$cLISI_celltype)

df_obj <- data.frame(barcode = seurat_obj$barcode,
                     celltype = seurat_obj$celltype,
                     sample = seurat_obj$sample,iLISI_sample = seurat_cca$iLISI_sample,
                     cLISI_celltype = seurat_cca$cLISI_celltype)

df_sct <- data.frame(barcode = scRNA_sct$barcode,
                     iLISI_sample = scRNA_sct$iLISI_sample,
                     cLISI_celltype = scRNA_sct$cLISI_celltype)

# Ensure each source has unique suffixes when merging to avoid column name conflicts
df_merged <- df_obj %>%
  left_join(df_harmony, by = "barcode", suffix = c("", ".harmony")) %>%
  left_join(df_cca, by = "barcode", suffix = c("", ".cca")) %>%
  left_join(df_sct, by = "barcode", suffix = c("", ".sct"))

# Check results
head(df_merged)

# The result is a new dataframe containing all information aligned by barcode
print(df_merged)
library(tidyr)
library(dplyr)
library(ggplot2)

# Convert iLISI and cLISI values to long format for plotting
df_long <- df_merged %>%
  pivot_longer(
    cols = c(iLISI_sample,cLISI_celltype,
        iLISI_sample.harmony, cLISI_celltype.harmony,
             iLISI_sample.cca, cLISI_celltype.cca,
             iLISI_sample.sct, cLISI_celltype.sct),
    names_to = c("Metric", "Method"),
    names_sep = "\\.",
    values_to = "Score"
  )
library(dplyr)

# Replace NA with "normal"
df_long <- df_long %>%
  mutate(Method = if_else(is.na(Method), "normal", Method))

# Check results
table(df_long$Method)

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
        annotate("text", x = Inf, y = Inf, label = col, hjust = 1.1, vjust = 2, size = 5, angle = 0) # Add y-axis title at the top
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
df_long$
p = vln_plot_beauty(meta_data = subset(df_long,Metric == 'iLISI_sample'), cols = 'Score', color = c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0'), group.by = 'Method', mod = "A")
ggsave('./plot/R2-32_iLISI_ca6.pdf',p,width = 5,height = 4)
p = vln_plot_beauty(meta_data = subset(df_long,Metric == 'cLISI_celltype'), cols = 'Score', color = c('#E5D2DD','#53A85F','#F1BB72','#F3B1A0'), group.by = 'Method', mod = "A")
ggsave('./plot/R2-32_cLISI_ca6.pdf',p,width = 5,height = 4)


## Dimensionality Reduction, Clustering

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj,nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj,npcs = 50,verbose = T)
# Plot PCA
#DimPlot(seurat_obj, reduction = "pca")
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
ElbowPlot(seurat_obj, ndims = 50)
seurat_obj <- RunUMAP(seurat_obj,reduction = "pca",dims = 1:10,return.model = TRUE ) 
seurat_obj <- FindNeighbors(seurat_obj,reduction = "pca",dims = 1:10)
seurat_obj <- FindClusters(seurat_obj,resolution = 0.5)


## Visualization

### dotplot

p = DotPlot(seurat_obj,features=hox_list,group.by = 'celltype') + scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#f6f6f6', '#f5f2f3', '#f0ebeb', '#F1E6E6', '#f8e7eb', '#dc8a9a', '#be3b56', '#d15c6c', '#8d192b'))+theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.border = element_rect(color="black"), panel.spacing = unit(1, "mm"))

