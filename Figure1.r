#code for Figure1
# Data integration and clustering
#       Seurat      5.1.0
#   SeuratDisk 0.0.0.9021
#      ggplot2      3.4.4
#    patchwork 1.2.0.9000
#        dplyr      1.1.4
#     magrittr      2.0.3
#      viridis      0.6.5
#  scCustomize      2.1.2
#           qs     0.25.7

## Integration and Formatting of Chip Data

library(Seurat)
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

## slingshot
library(Seurat)
library(SCP)
library(ggplot2)

# Load Seurat object with slingshot trajectory
seurat_g1 <- readRDS('./data/seurat_g1_slingshot_629.rds')

# Run Slingshot for trajectory analysis
seurat_g1 <- RunSlingshot(srt = seurat_g1, group.by = "order_617", reduction = "umap", start = 'Inter.Epi')

# Plot trajectories on UMAP
p <- CellDimPlot(
  seurat_g1, 
  group.by = "order_617", 
  reduction = "umap", 
  lineages_span = 0.5, 
  lineages = paste0("Lineage", 1:3),
  palcolor = mycols[levels(seurat_my$order_617)],
  lineages_palcolor = c("#fe9929","#54278f", "#1c9099")
) 
ggsave('./plot/cs6_epi_slingshot_1.pdf', p, width = 7, height = 4.75)

# Plot lineage features on UMAP
p <- FeatureDimPlot(
  seurat_g1, 
  features = paste0("Lineage", 1:3), 
  reduction = "UMAP", 
  theme_use = "theme_blank"
)
ggsave('./plot/cs6_epi_slingshot_2.pdf', p, width = 17, height = 5)

# Identify dynamic features along trajectories
seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage1", "Lineage2", "Lineage3"), 
  n_candidates = 50,
  BPPARAM = BPPARAM
)

# Function to extract and sort top 50 genes
get_top_50_genes <- function(lineage_data) {
  # Sort by padjust value
  sorted_genes <- lineage_data[order(lineage_data$padjust), ]
  # Extract top 50 gene names
  top_50_genes <- rownames(sorted_genes)[1:50]
  return(top_50_genes)
}

# Run DynamicFeatures with custom gene sets
seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage1", "Lineage2", "Lineage3"),
  features = unique(amnion_markers)
)

seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage3"),
  features = unique(lineage3_genes),
  BPPARAM = BPPARAM
)

# Extract top dynamic feature genes for each lineage
lineage1_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage1$DynamicFeatures)
lineage2_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage2$DynamicFeatures)
lineage3_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage3$DynamicFeatures)
lineage3_genes <- c(lineage3_genes, 'CGB', 'OTX2', 'UCHL1', 'TBXT', 'POU5F1')

# Define gene sets for analysis
amnion_markers <- c("ISL1", "TFAP2B", "TFAP2A", "WNT6", "GABRP", "HEY1", "BAMBI", "DLX5", "SOX4", "GRHL1", 
                    "MEIS1", "BMP4", "PRTG", "DSP", "TGFBI", "MEST", "MSX2", "MITF", "VTCN1", "IGFBP3", 
                    "PRKD1", "KCNMA1", "STC1", "TCF4", "HAND1", "WNT11", "TBX3", "CDX2", "FURIN", "GATA6", 
                    "SALL4", "KRT19", "AQP3", "AMOTL1", "DAB2", "EMP2", "CA12", "CITED2", "TEAD1", "FOLR1", 
                    "PTGES")

lienage2_genes <- c(
  "CDX2", "HAND1", "PTN", "PGF", "CD55", "TBX3", "GATA6", "SALL4", "CITED2", "WNT11", "FURIN", 
  "KRT19", "AQP3", "AMOTL1", "DAB2", "EMP2", "CA12", "TEAD1", "FOLR1", "PTGES", "GABRP", "VTCN1", 
  "TFAP2A", "TFAP2B", "IGFBP3", "PRKD1", "KCNMA1", "STC1", "TCF4", "ISL1", "DLX5", "PRTG", "TGFB1", 
  "WNT6", "BMP4", "MSX2", "DSP", "MEIS1", "SOX4", "BAMBI", "HEY1", "MEST", "GRHL1", "MITF"
)

EMT_genes <- c("AXL", "BMI1", "CDH1", "CDKN2A", "E2-2", "E47", "EHMT2", "EPAS1", "ETS1", "EZH2", 
               "FOXC2", "G9A", "GSC", "HDAC1", "HDAC2", "HDAC3", "HIF1A", "KLF8", "LOX", "LOXL2", 
               "LSD1", "MCT1", "PIK3CA", "PRRX1", "PTEN", "RB1", "SIP1", "SIRT1", "SIX1", "SLUG", 
               "SMAD2", "SNAI1", "SNAI2", "SNAIL", "SNAIL1", "SNAIL2", "SUV39H1", "SUZ12", "T", 
               "TCF4", "TP53", "TWIST", "TWIST1", "ZEB1", "ZEB2")

# PRCP genes analysis
prcp_genes <- c("GDF3", "WNT3A", "FGF17", "GSC", "SAT1", "LHX1", "CHRD", "POU5F1", "NOTO", "HHEX", 
                "NODAL", "FOXA2", "MIXL1", "SP5", "TBXT", "SHISA2", "CDH2", "FZD8", 
                "MT1H", "SIX3", "OTX2", "SOX17", "SFRP1", "CER1", "DKK1", "MESP1", "FOXJ1", "DKK4")

prcp_genes <- intersect(prcp_genes, rownames(seurat_g1))

library(BiocParallel)
# Set BPPARAM for single-thread processing
BPPARAM <- SerialParam()

seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage1"),
  features = unique(prcp_genes),
  BPPARAM = BPPARAM
) 

seurat_l1 <- subset(seurat_g1, Lineage1 >= 0)

# Plot dynamic features for PRCP lineage
p <- DynamicPlot(
    srt = seurat_l1, 
    lineages = c("Lineage1"), 
    group.by = "celltype", # Don't use 'order' as it may cause errors
    features = prcp_genes,
    compare_lineages = TRUE, 
    compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#fe9929",  # Line color: #1c9099 for Lineage3, #54278f for Lineage2, #fe9929 for Lineage1
    point_palcolor = list(order_617 = mycols[levels(seurat_l1$order_617)][c('Inter.Epi', 'Anterior.Epi', 'PrCP')]),
    add_line = TRUE,
    add_interval = TRUE,
    line.size = 1,
    add_point = TRUE,
    add_rug = TRUE,
    flip = FALSE,
    reverse = FALSE,
    x_order = "value",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scp",
    ncol = 8
)
ggsave('./plot/lingage1_genes_prcp_1.pdf', p, width = 24, height = 12)

# EMT genes analysis
seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage3"),
  features = unique(emt_genes),
  BPPARAM = BPPARAM
) 

p <- DynamicPlot(
    srt = seurat_g1, 
    lineages = c("Lineage3"), 
    group.by = "order_617",
    features = emt_genes,
    compare_lineages = TRUE, 
    compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#1c9099",  # Line color: #1c9099 for Lineage3
    point_palcolor = mycols[levels(seurat_g1$order_617)][c('PrCP', 'Anterior.Epi', 'Inter.Epi', 'Posterior.Epi', 'AM.Ecto', 'AM')],
    add_line = TRUE,
    add_interval = TRUE,
    line.size = 1,
    add_point = TRUE,
    add_rug = TRUE,
    flip = FALSE,
    reverse = FALSE,
    x_order = "value",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scp",
    ncol = 8
)
ggsave('./plot/lingage3_genes_emt_test.pdf', p, width = 24, height = 14)

# Combine gene sets for Lineage2 analysis
genes <- unique(c(lienage2_genes, amnion_markers))

seurat_g1_s1 <- RunDynamicFeatures(
  srt = seurat_g1_s1, 
  lineages = c("Lineage2"),
  features = unique(genes),
  BPPARAM = BPPARAM
)

# Plot dynamic features for Lineage2
p <- DynamicPlot(
    srt = seurat_g1_s1, 
    lineages = c("Lineage2"), 
    group.by = "celltype",
    features = genes,
    compare_lineages = TRUE, 
    compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#54278f",  # Line color: #54278f for Lineage2
    point_palcolor = mycols[levels(seurat_g1$order_617)][c('PrCP', 'Anterior.Epi', 'Inter.Epi', 'Posterior.Epi', 'AM.Ecto', 'AM')],
    add_line = TRUE,
    add_interval = TRUE,
    line.size = 1,
    add_point = TRUE,
    add_rug = TRUE,
    flip = FALSE,
    reverse = FALSE,
    x_order = "value",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scp",
    ncol = 8,
    heatmap_palcolor = heatcolor
)
ggsave('./plot/lingage2_genes_2.pdf', p, width = 24, height = 17)


## pyScenic

### create_loom_input.R

library(optparse)
op_list <- list(
make_option(c("-i", "--inrds"), type = "character", default = NULL, action = "store", help = "The input of Seurat RDS",metavar="rds"),
make_option(c("-d", "--ident"), type = "character", default = NULL, action = "store", help = "The sample Ident of Seurat object",metavar="idents"),
make_option(c("-s", "--size"),  type = "integer", default = NULL, action = "store", help = "The sample size of Seurat object",metavar="size"),
make_option(c("-l", "--label"), type = "character", default = "out", action = "store", help = "The label of output file",metavar="label"),
make_option(c("-a", "--assay"), type = "character", default = "Spatial", action = "store", help = "The assay of input file",metavar="assay")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

assay <- opt$assay

library(Seurat)
obj <- readRDS(opt$inrds)
if (!is.null(opt$ident)) {
Idents(obj) <-  opt$ident
size=opt$size
if (!is.null(size)) {
obj <- subset(x = obj, downsample = opt$size)
}
saveRDS(obj,"subset.rds")
}
if (is.null(opt$label)) {
label1 <- 'out'
}else{
label1 <- opt$label
}

library(SCopeLoomR)
outloom <- paste0(label1,".loom")
build_loom(file.name = outloom,dgem = obj@assays[[assay]]@counts)
write.table(obj@meta.data,'metadata_subset.xls',sep='\t',quote=F)



### pyscenic_from_loom.sh

``` sh
input_loom=out.loom
n_workers=20
#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-h|--help:\t Usage information"
exit 1
}
#get value
while getopts :i:n:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done
#database path
tfs=/home/xlyang/python_work/pyscenic_pipeline/01.database/human_base/hs_hgnc_tfs.txt
feather=/home/xlyang/python_work/pyscenic_pipeline/01.database/human_base/*.feather
tbl=/home/xlyang/python_work/pyscenic_pipeline/01.database/human_base/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
pyscenic=/sdc/xlyang/software/anaconda3/envs/pyscenic/bin/pyscenic

# grn
 $pyscenic grn \
 --num_workers $n_workers \
 --output grn.tsv \
 --method grnboost2 \
 $input_loom  $tfs

# cistarget
$pyscenic ctx \
grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
ctx.csv \
--output aucell.loom \
--num_workers $n_workers
```

### calcRSS_by_scenic.R

library(optparse)
op_list <- list(
make_option(c("-l", "--input_loom"), type = "character", default = NULL, action = "store", help = "The input of aucell loom file",metavar="rds"),
make_option(c("-m", "--input_meta"), type = "character", default = NULL, action = "store", help = "The metadata of Seurat object",metavar="idents"),
make_option(c("-a", "--assay"), type = "character", default = 'Spatial', action = "store", help = "The assay of Seurat object",metavar="assay"),
make_option(c("-c", "--celltype"), type = "character", default = NULL, action = "store", help = "The colname of metadata to calculate RSS",metavar="label")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
celltype <- opt$celltype
message(paste0('celltype group: ',celltype))
assay <- opt$assay
loom <- open_loom(opt$input_loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

meta <- read.table(opt$input_meta,sep='\t',header=T,stringsAsFactor=F)
cellinfo <- meta[,c(opt$celltype,paste0("nFeature_",assay),paste0("nCount_",assay))]
colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"

sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
try({
rssPlot <- plotRSS(rss)
save(regulonAUC,rssPlot,regulons,file='regulon_RSS.Rdata')
})

saveRDS(rss,paste0(celltype,"_rss.rds"))

source('/home/xlyang/python_work/pyscenic_pipeline/00.scripts/function_pyscenic_visualize.R')
plot_pyscenic(inloom='aucell.loom',incolor=incolor,inrss=paste0(celltype,"_rss.rds"),inrds='subset.rds',infun='median', ct.col=celltype,inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50)

### Execute the above script to run pyscenic

``` R
inscp='../00.scripts/'
inrds='./seurat_g1_slingshot_629.rds'
ingroup='celltype'

# Step 1: Activate snapatacv1 environment and run the first R script
source /home/xlyang/software/anaconda3/bin/activate snapatacv1
Rscript ${inscp}/create_loom_input.R -i ${inrds} -d ${ingroup} -l out -a Spatial
sleep 30

# Step 2: Activate pyscenic environment and run the shell script
source /home/xlyang/software/anaconda3/bin/activate pyscenic
sh ${inscp}/pyscenic_from_loom.sh -i out.loom -n 10

# Step 3: Reactivate snapatacv1 environment and run the second R script
source /home/xlyang/software/anaconda3/bin/activate snapatacv1
Rscript ${inscp}/calcRSS_by_scenic.R -l aucell.loom -m metadata_subset.xls -c ${ingroup} -a Spatial

# Optionally, deactivate the environment at the end of the script
conda deactivate
```

### Visualize the RSS

``` R
library(Seurat)
library(pheatmap)

df <- readRDS('../celltype2_rss.rds')

# Define Min-Max normalization function
min_max_normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

# Perform Min-Max normalization for each row (each gene)
#htdf <- t(apply(df, 1, min_max_normalize))
htdf <- scale(df, scale = TRUE, center = TRUE)
# Calculate row standard deviation and filter
row_sds <- apply(htdf, 1, sd)
filtered_mat <- htdf[row_sds > 0.95, ]
htdf <- filtered_mat

# Create column annotation
Cell_Type <- colnames(df)
annotation_col_df <- data.frame(Cell_Type = Cell_Type)
rownames(annotation_col_df) <- Cell_Type

# Define column order
#custom_order <- c("CTB_pri", "CTB_ccc", "iEVT_naïve", "iEVT_mature", "GC")
htdf <- htdf[, custom_order]
annotation_col_df <- annotation_col_df[custom_order, , drop = FALSE]

# Define colors for cell types
#mycols <- c('CTB_pri'='#35648F', 'CTB_ccc'='#BF3F45', 'iEVT_naïve'='#E29F46', 'iEVT_mature'='#BAC76B', 'GC'='#8C29C9')
colors_list <- list('Cell_Type' = mycols)

# Create heatmap
p <- pheatmap(htdf, 
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              annotation_col = annotation_col_df,
              annotation_colors = colors_list,
              show_rownames = TRUE,
              show_colnames = TRUE,
              color = colorRampPalette(c('#253494', '#2c7fb8', '#c7e9b4', '#ffffcc'))(100),
              main = "RSS Score (Z-Score Normalized)")

print(p)
pdf('./circular_heatmap.pdf',width = 5,height = 7)
print(p)
dev.off()

library(dendextend) # Add this line
# Load necessary packages
library(circlize)
library(ComplexHeatmap)

# Define color mapping
mycol2 <- colorRamp2(c(-1.7, 0.3, 2.3), c("#57ab81", "white", "#ff9600"))

# Draw circular heatmap
circos.clear()
circos.par(gap.after = c(22))
circos.heatmap(hfdf, col = mycol2, dend.side = "inside", rownames.side = "outside", 
               track.height = 0.38, rownames.col = "black", rownames.cex = 0.9, 
               rownames.font = 1, cluster = TRUE, dend.track.height = 0.18,
               dend.callback = function(dend, m, si) {
                 color_branches(dend, k = length(mycols), col = mycols)
               })

# Add legend
lg <- Legend(title = "Exp", col_fun = mycol2, direction = "vertical")
grid.draw(lg)

# Add column names
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 1) {
    cn <- colnames(hfdf)
    n <- length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.8, "mm"),
                7.8 + (1:n) * 1.1,
                cn, cex = 0.8, adj = c(0, 1), facing = "inside")
  }
}, bg.border = NA)

circos.clear()
```



## FASTMNN

``` R
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(SeuratDisk)

# Load Seurat object and metadata
seurat_obj <- LoadH5Seurat('./data/cs6_all.h5seurat')
meta_cs6 <- readRDS('./data/CS6_all_meta_250402.rds')
head(meta_cs6)

# Add metadata and subset cells
seurat_obj <- AddMetaData(seurat_obj, meta_cs6)
unique(seurat_obj@meta.data$celltype_new)
seurat_obj <- subset(seurat_obj, celltype_new %in% c('Posterior.Epi', 'Hypo.2', 'Hypo.1', 'ave', 'Anterior.Epi', 'Inter.Epi', 'Anterior pole'))

# Set RNA assay as default and remove Spatial assay
seurat_obj[['RNA']] = seurat_obj[['Spatial']]
DefaultAssay(seurat_obj) = 'RNA'
seurat_obj[['Spatial']] = NULL
seurat_obj$Day = 'CS6'
seurat_obj$celltype = seurat_obj$celltype_new

# Load and process LTQ data
seurat_ltq <- readRDS('/home/xlyang/download/refer_data/ltx_14d/ltq_14.rds')
unique(seurat_ltq@meta.data$Day)
unique(seurat_ltq@meta.data$Group)
seurat_ltq$celltype <- seurat_ltq@meta.data$Group
seurat_ltq <- subset(seurat_ltq, Day %in% c('D12', 'D14') & Group %in% c('EPI', 'PrE', 'PSA-EPI'))

# Define color palettes
ltq_col <- c('EPI' = '#813a90', 'PSA-EPI' = '#904410', 'PrE' = '#62509d')
CS6_color <- c('Anterior pole' = '#E86c79',  # "PrCP"
               "Anterior.Epi" = '#253494',
               "Inter.Epi" = '#006837', "Posterior.Epi" = '#1d91c0',
               "AM.Ecto" = '#4682B4', "AM" = '#bcbddc',
               "AM.EXMC" = '#9e9ac8', "Connecting Stalk" = '#fdd0a2',
               "Hypo.1" = '#6a51a3', "Hypo.2" = '#ae017e',
               "PL.EXMC" = '#ccebc5', "CTB" = '#2b8cbe',
               "CTB.Fusion" = '#7bccc4', "STB" = '#58BCE8',
               "MTB" = '#084081', "YS.Endo_1" = '#FFAA92',
               "YS.Endo_2" = '#8FB0FF', "YS.Endo_3" = '#00C2A0',
               "YS.EXMC_1" = '#6F0062', "YS.EXMC_2" = '#EEC3FF',
               "Myeloid Progenitor" = '#D16100', "Primitive_Ery1" = 'magenta',
               "Primitive_Ery2" = '#B79762', "Primitive_Mk1" = 'red3',
               "Primitive_Mk2" = 'purple', "YS.Endo" = '#ec7014',
               "YS.EXMC" = '#993404', "Blood" = '#e31a1c',
               "grey" = '#EBEBEB', "Epi" = '#006837', "Hypo" = '#6a51a3',
               'ave' = 'green')
mycols <- c(ltq_col, CS6_color)

# Prepare datasets for integration
seuratcs6 = seurat_obj
seuratcs6 <- subset(seuratcs6, features = rownames(seurat_ltq))
seurat_ltq <- subset(seurat_ltq, features = rownames(seuratcs6))
seurat_ltq$sample <- 'E12-E14'
seuratcs6$sample <- 'cs6'

# Merge datasets and prepare for integration
seurat_obj = merge(seuratcs6, seurat_ltq)
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")
seurat_obj <- JoinLayers(seurat_obj)

# Normalize data and select integration features
seurat_obj <- NormalizeData(seurat_obj)
features <- SelectIntegrationFeatures(object.list = SplitObject(seurat_obj, 'sample'))
VariableFeatures(seurat_obj) = features

# Split objects for integration
split_objects <- SplitObject(seurat_obj, split.by = "sample")
print(names(split_objects))

# Reorder objects to ensure cs6 is first
split_objects <- split_objects[c("cs6", "E12-E14")]

# Run FastMNN integration
seurat_obj <- RunFastMNN(
  object.list = split_objects,
  k = 8,
  merge.order = c("cs6", "E12-E14"),
  auto.merge = FALSE
)

# Run UMAP and visualize results
seurat_obj <- RunUMAP(seurat_obj, reduction = "mnn", dims = 1:10)

library(scplotter)
p1 = CellDimPlot(
  seurat_obj,
  group_by = "celltype",
  reduction = "umap",
  label = TRUE,
  palcolor = mycols,
  pt.size = 1.5,
  sizes.highlight = 2,
  theme = "theme_blank",
  legend.position = "right",
  raster = FALSE,
  highlight = 'sample == "E12-E14"'
)

# Save plot
library(ggplot2)
ggsave('./plot/e1214_cs6_fmnn_p1.pdf', p1, width = 6.5, height = 5.75)
```