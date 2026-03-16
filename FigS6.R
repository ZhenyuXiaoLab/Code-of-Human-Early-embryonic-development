## code for Fig.S6
# Seurat	5.1.0
# tidyr	1.3.1
# dplyr	1.1.4
# tibble	3.2.1.9016
# ggplot2	3.4.4
### integrate
######integrative analysis between CS6 and CS7 ######
#load CS6 data
load('cs6_analysis_0527.Rdata')
figure1_anno <- readRDS('fig1_metadata_823.rds')
ana1$Cluster <- as.character(figure1_anno[colnames(ana1), 'celltype'])
ana1$Cluster[colnames(ana1) %in% colnames(cs6_hema2) ] <- as.character(na.omit(as.character(cs6_hema2$cluster_final)[match(colnames(ana1), colnames(cs6_hema2))]))

cs6_hypo_epi2 <- ana1[, ana1$Cluster %in% c('Hypo1','Hypo.2',  "PrCP",  "Anterior.Epi", "Inter.Epi", "Posterior.Epi", 'AM.EXMC','PL.EXMC','CS',"YS.Endo1", "YS.Endo2", "YS.Endo3", "YS.EXMC1" ,"YS.EXMC2", "pMeg1", "pMeg2", "pEry1", 'pEry2','pMP')]
cs6_hypo_epi2$stage <- 'CS6'
cs6_hypo_epi2$batch <- 'CS6'
cs6_hypo_epi2$cluster <- as.character(cs6_hypo_epi2$Cluster)

#downsample cells of CS6
cs6_hypo_epi_cells2 <- Reduce(union, purrr::map(unique(cs6_hypo_epi2$cluster), function(x){
  
  cells <- colnames(cs6_hypo_epi2)[cs6_hypo_epi2$cluster==x]
  d
  if(length(cells) >200){
    
    cells <- sample(cells, 200)
  }
  
  return(cells)
  
}))

#CS7 data Tyser
cs7_ps2 <- n2021[, as.character(n2021$sub_cluster) %in% c('Epiblast', 'Hypoblast','YS Endoderm', 'Primitive Streak', 'Nascent Mesoderm', 'Emergent Mesoderm','Axial Mesoderm', 'Advanced Mesoderm',
                                            'YS.EXMC',"Primitive_Mk", "pMP", 'YSMP', 'Primitive_Ery', 'Mac')]
cs7_ps2$stage <- 'CS7'
cs7_ps2$batch <- 'CS7'
cs7_ps2$cluster <- as.character(cs7_ps2$sub_cluster)


DefaultAssay(cs6_hypo_epi2) <- 'RNA'
DefaultAssay(cs7_ps2) <- 'RNA'

#integrate data
co_genes6 <- intersect(rownames(cs6_hypo_epi2), rownames(cs7_ps2))

cs6_cs7_comv2 <- merge(cs6_hypo_epi2[co_genes6, cs6_hypo_epi_cells2],
                       cs7_ps2[co_genes6, ])

cs6_cs7_comv2 <- NormalizeData(cs6_cs7_comv2)
cs6_cs7_comv2 <- FindVariableFeatures(cs6_cs7_comv2)
cs6_cs7_comv2@assays$RNA@scale.data <- as.matrix(0)
cs6_cs7_comv2 <- RunFastMNN(SplitObject(cs6_cs7_comv2, split.by = 'batch')[c('CS7', 'CS6')], k=50)
cs6_cs7_comv2 <- RunUMAP(cs6_cs7_comv2, reduction = 'mnn', dims = 1:15,
                        local.connectivity = 1, n.neighbors = 15)

cs6_cs7_comv2$Cluster <- cs6_cs7_comv2$cluster
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('Hypo1', 'Hypo.2')] <- 'Hypoblast'
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('YS.EXMC1', 'YS.EXMC2')] <- 'YS.EXMC'
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('PrCP')] <- 'Anterior Pole'
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('Posterior.Epi')] <- 'Gast-primed.Epi'
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('Inter.Epi')] <- 'Posterior.Epi'
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('Myeloid Progenitor','pMP', 'Primitive_Mk', 'Primitive_Ery','pEry1','pEry2','pMeg1','pMeg2',
                                                   'YSMP','Mac',"Primitive_Mk", "Primitive_Ery")] <- 'Blood'
cs6_cs7_comv2$cluster[cs6_cs7_comv2$cluster %in% c('YS.Endo1', 'YS.Endo2', 'YS.Endo3', 'YS Endoderm')] <- 'YS.Endo'
cs6_cs7_comv2$com <- paste(cs6_cs7_comv2$cluster, cs6_cs7_comv2$stage, sep="_")
cs6_cs7_comv2$com <- factor(cs6_cs7_comv2$com, levels = c('Anterior Pole_CS6', "Anterior.Epi_CS6", "Posterior.Epi_CS6", "Gast-primed.Epi_CS6", "Epiblast_CS7", "Axial Mesoderm_CS7",
                                                           "Hypoblast_CS6", "Hypoblast_CS7","YS.Endo_CS6", "YS.Endo_CS7", "Primitive Streak_CS7",
                                                           "Nascent Mesoderm_CS7", "Emergent Mesoderm_CS7", "Advanced Mesoderm_CS7",
                                                           'YS.EXMC_CS6','YS.EXMC_CS7', 'AM.EXMC_CS6','PL.EXMC_CS6','CS_CS6', 'Blood_CS6', 'Blood_CS7'))
scales::show_col(col_cs6_cs7_com)
col_cs6_cs7_com <- c(col_cs6_cs7_com[1:19], 'PL.EXMC_CS6' = 'orange3', 'CS_CS6' = '#6B7900')

cs6_cs7_comv2@reductions$umap@cell.embeddings[ ,1] <- 0-cs6_cs7_comv2@reductions$umap@cell.embeddings[ ,1]
DimPlot(cs6_cs7_comv2, group.by = 'com', cols = col_cs6_cs7_com, label = F, repel = T)

cs6_cs7_comv2$com_cs6 <- as.character(cs6_cs7_comv2$com)
cs6_cs7_comv2$com_cs6[cs6_cs7_comv2$stage=='CS7'] <- 'others'
cs6_cs7_comv2$com_cs6 <- factor(cs6_cs7_comv2$com_cs6, levels = c('Anterior Pole_CS6', "Anterior.Epi_CS6", "Posterior.Epi_CS6", "Gast-primed.Epi_CS6",
                                "Hypoblast_CS6", "YS.Endo_CS6", 
                                'YS.EXMC_CS6','AM.EXMC_CS6','PL.EXMC_CS6','CS_CS6', 'Blood_CS6', 'others'))

cs6_cs7_comv2$com_cs7 <- as.character(cs6_cs7_comv2$com)
cs6_cs7_comv2$com_cs7[cs6_cs7_comv2$stage=='CS6'] <- 'others'
cs6_cs7_comv2$com_cs7 <- factor(cs6_cs7_comv2$com_cs7, levels = c("Epiblast_CS7", "Axial Mesoderm_CS7",
                                                                  "Hypoblast_CS7", "YS.Endo_CS7", "Primitive Streak_CS7",
                                                                  "Nascent Mesoderm_CS7", "Emergent Mesoderm_CS7", "Advanced Mesoderm_CS7",
                                                                  'YS.EXMC_CS7', 'Blood_CS7','others'))


DimPlot(cs6_cs7_comv2, group.by = 'com', cols = c(col_cs6_cs7_com))
DimPlot(cs6_cs7_comv2, group.by = 'com_cs6', cols = c('others' = '#E6E6E6', col_cs6_cs7_com))+
DimPlot(cs6_cs7_comv2, group.by = 'com_cs7', cols = c('others' = '#E6E6E6', col_cs6_cs7_com))


library(Seurat)
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
	### vlnplot
big_vln_plot <- function (seurat_obj, 
                          genes = c(), 
                          assay = "RNA", 
                          slot = "data",
                          reverse = TRUE, 
                          cols = NULL, 
                          add_line = TRUE, 
                          group.by = NULL,
                          cluster_col = FALSE, 
                          cluster_row = FALSE, 
                          group.level = NULL,
                          hide_triangle = FALSE) {
  
  require(Seurat)
  require(tidyr)
  require(dplyr)
  require(tibble)
  require(ggplot2)
  
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) {
    stop("No valid genes provided or found in the Seurat object.")
  }
  
  if (slot == "scale.data") {
    seurat_obj <- ScaleData(seurat_obj, features = genes, assay = assay)
  }
  
  vars_to_fetch <- c(genes, if (is.null(group.by)) "ident" else group.by)
  dt <- FetchData(seurat_obj, vars = vars_to_fetch, slot = slot, assays = assay)
  
  dt <- dt %>% 
    rownames_to_column("cell") %>% 
    pivot_longer(cols = all_of(genes), names_to = "genes", values_to = "expressions") %>% 
    rename(cluster = if (is.null(group.by)) "ident" else group.by)
  
  dt$genes <- factor(dt$genes, levels = genes)
  
  if (!is.null(group.level)) {
    dt$cluster <- factor(dt$cluster, levels = group.level)
  }
  
  if (is.null(cols)) {
    cols <- (scales::hue_pal())(length(unique(dt$cluster)))
  }
  
  if (cluster_col || cluster_row) {
    require(ggdendro)
    
    if (reverse) {
      if (cluster_col) {
        gene_dist <- dist(t(dt %>% pivot_wider(names_from = genes, values_from = expressions) %>% select(-cell, -cluster)))
        gene_hclust <- hclust(gene_dist)
        gene_order <- gene_hclust$labels[gene_hclust$order]
        dt$genes <- factor(dt$genes, levels = gene_order)
      }
    } else {
      if (cluster_row) {
        cluster_data <- dt %>% 
          group_by(cluster, genes) %>% 
          summarise(mean_expression = mean(expressions, na.rm = TRUE), .groups = "drop") %>% 
          pivot_wider(names_from = genes, values_from = mean_expression) %>% 
          column_to_rownames("cluster")
        
        if (any(is.na(cluster_data))) {
          warning("NAs found in cluster data. Removing NAs for clustering.")
          cluster_data <- na.omit(cluster_data)
        }
        
        if (nrow(cluster_data) > 1) {
          cluster_dist <- dist(cluster_data)
          cluster_hclust <- hclust(cluster_dist)
          cluster_order <- rownames(cluster_data)[cluster_hclust$order]
          dt$cluster <- factor(dt$cluster, levels = cluster_order)
        } else {
          warning("Not enough data for clustering after removing NAs.")
        }
      }
    }
  }
  
  max_expr_clusters <- dt %>% 
    group_by(genes, cluster) %>% 
    summarise(mean_expression = mean(expressions, na.rm = TRUE), .groups = "drop") %>% 
    group_by(genes) %>% 
    summarise(max_cluster = cluster[which.max(mean_expression)], max_expression = max(mean_expression))
  
  p <- ggplot(data = dt, aes(x = if (reverse) expressions else genes, 
                             y = if (reverse) cluster else expressions, 
                             fill = cluster)) +
    geom_violin(scale = "width", draw_quantiles = if (add_line) c(0.25, 0.5, 0.75) else NULL, 
                color = "black", size = 0.45, alpha = 0.8) +
    scale_fill_manual(values = cols, limits = rev(levels(dt$cluster))) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(size = 14, angle = if (reverse) 60 else 45, hjust = if (reverse) 1 else 1), 
          axis.text.y = element_text(size = 16), 
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16), 
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 16, angle = 60), 
          legend.title = element_text(size = 16), 
          legend.text = element_text(size = 15))
  
  custom_breaks <- function(x) { c(0, max(x)) }
  custom_labels <- function(x) { c("0", format(round(max(x), 1), nsmall = 1)) }
  
  if (reverse) {
    p <- p + 
      facet_grid(cols = vars(genes), scales = "free_x") +
      labs(x = "Log Normalized Expression", y = "") +
      scale_x_continuous(breaks = custom_breaks, labels = custom_labels)
  } else {
    p <- p + 
      facet_grid(rows = vars(cluster), scales = "free_y") +
      labs(x = "", y = "Log Normalized Expression") +
      scale_y_continuous(breaks = custom_breaks, labels = custom_labels)
  }
  
  if (!add_line && !hide_triangle) {
    star_data <- dt %>% 
      group_by(genes, cluster) %>% 
      summarise(mean_expression = mean(expressions), .groups = "drop") %>% 
      left_join(max_expr_clusters, by = "genes") %>% 
      filter(cluster == max_cluster)
    
    if (reverse) {
      p <- p + 
        geom_point(data = star_data, aes(x = mean_expression, y = cluster), 
                   shape = 24, size = 2.5, color = "black", fill = "black")
    } else {
      p <- p + 
        geom_point(data = star_data, aes(x = genes, y = mean_expression), 
                   shape = 24, size = 2.5, color = "black", fill = "black")
    }
  }
  
  options(repr.plot.width = 12, repr.plot.height = if (reverse) 4 else 12)
  return(p)
}

p <- big_vln_plot(seurat_obj, 
                  genes = unique(c("GATA4","GATA6","MYL7","LUM", "ISM2","GATA3","GJA5", "KRT8","CGA", 
                                   "KRT18", "INSL4", "PRG2", "MMP12", "PPBP", "DKK1", "CXCL8", "PLAC8", 
                                   "CD82", "KRT7", "ERVFRD-1", "EPCAM", "SOX9","LGR5","PECAM1","ACKR1", 
                                   "RGS5", "KDR", "PTPRC", "PRAP1", "IGF1","HAND2","PAEP", 'ACTA2',
                                   'MYH11','VIM','DCN',"PCNA","MKI67")), 
                  slot = 'scale.data', 
                  reverse = TRUE, 
                  add_line = FALSE, 
                  cols = rev(mycols[unique(seurat_obj$celltype)]), 
                  group.by = 'celltype', 
                  cluster_row = FALSE, 
                  cluster_col = FALSE, 
                  group.level = unique(seurat_obj$celltype), 
                  assay = 'Spatial')
```

### GO

same as the code provided in EXT.Fig3

### project 

same as the code provided in EXT.Fig5, but we transform the mouse/monkey gene name to the same as human genes first.
