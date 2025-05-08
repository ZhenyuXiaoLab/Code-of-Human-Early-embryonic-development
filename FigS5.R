#code for Fig.S5

### integrate
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



# Mapping prediction
anchors <- FindTransferAnchors(reference = seurat_obj_ref, query = seurat_obj_my, reference.reduction = "pca",k.filter = NA,mapping.score.k = 50,dims = 1:10,k.anchor = 10,k.score = 10) #,k.filter = NA
seurat_obj_my <- MapQuery(anchorset = anchors, reference = seurat_obj_ref, query = seurat_obj_my,
                          refdata = "celltype", reference.reduction = "pca", 
                          reduction.model = ifelse(use_raw_umap, "umapnew", "umap"))
    
meta_my <- cbind(seurat_obj_my@meta.data, as.data.frame(seurat_obj_my@reductions$ref.umap@cell.embeddings))
meta_ref <- cbind(seurat_obj_ref@meta.data, as.data.frame(seurat_obj_ref@reductions[[ifelse(use_raw_umap, "umapnew", "umap")]]@cell.embeddings))
    
options(repr.plot.width=7, repr.plot.height=6)
p1 <- Dimplot_sort(seurat_obj_ref, group_name = c('celltype'), order_list = 'auto', save_name='', color_palette = mycols, reduction = "umap")
    
library(ggplot2)
p2 <- ggplot(meta_ref, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey", size = 0.01, alpha = 0.5) +
  geom_point(data = meta_my, aes(x = refUMAP_1, y = refUMAP_2, color = celltype), size = 0.75, alpha = 1, shape = 17) +
  scale_color_manual(values = mycols) +
  theme_minimal() +
  labs(x = "", y = "", color = "Cell Type") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
p <- wrap_plots(plotlist = list(p1,p2), nrow = 1)
}
else {
  message('The plot_umap parameter must be either onmy, reverse or project')
  return(NULL)
}


### integrate

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