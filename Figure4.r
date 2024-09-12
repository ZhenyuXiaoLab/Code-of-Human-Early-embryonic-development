#######subset hematopoietic related clusters #############
cs6_hema <- ana1[, ana1$Main_CellType %in% c('VE', 'YS.Endo', 'YS.EXMC', 'EXMC_Prog', "Myeloid Progenitor",
                                             'Primitive Ery1', 'Primitive Ery2', 'Primitive Mk1','Primitive Mk2')]


cs6_hema <- NormalizeData(cs6_hema)
cs6_hema <- FindVariableFeatures(cs6_hema)
cs6_hema <- ScaleData(cs6_hema)
cs6_hema <- RunPCA(cs6_hema)

ElbowPlot(cs6_hema, ndims = 30)
cs6_hema <- RunUMAP(cs6_hema, dims = 1:20, min.dist = 0.3, spread = 0.5)
cs6_hema$CellType_V2[cs6_hema$CellType_V2=='VE1'] <- 'VE'

col_hema <- setNames(col.use[1:12], c('VE', 'YS.Endo_1','YS.Endo_2','YS.Endo_3', 'YS.EXMC_1','YS.EXMC_2', 'EXMC_Prog', "Myeloid Progenitor",
                                      'Primitive_Ery1', 'Primitive_Ery2', 'Primitive_Mk1','Primitive_Mk2'))

cs6_hema$CellType_V2 <- factor(cs6_hema$CellType_V2, levels = c('VE', 'YS.Endo_1','YS.Endo_2','YS.Endo_3', 'YS.EXMC_1','YS.EXMC_2', 'EXMC_Prog', "Myeloid Progenitor",
                                                                'Primitive_Ery1', 'Primitive_Ery2', 'Primitive_Mk1','Primitive_Mk2'))

DimPlot(cs6_hema, group.by = 'CellType_V2', cols = col.use, label = T, repel = T, pt.size = 1)

######Feature Genes 
cs6_hema <- SetIdent(cs6_hema, value = 'CellType_V2')
degs_cs6_hema <- FindAllMarkers(cs6_hema, logfc.threshold = log(1.25), only.pos = T, max.cells.per.ident = 1000)
degs_cs6_hema <- degs_cs6_hema[degs_cs6_hema$p_val_adj<0.05, ]

DotPlot(cs6_hema, features = rev(sc_tl_topgene(degs_cs6_hema, 5)),
        group.by = 'CellType_V2')+coord_flip()+
        scale_color_gradientn('Pseudotime', colours = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))+
        theme(axis.text.x = element_text(angle = 60, hjust = 1))

sc_tl_average <- function(data, cluster){
  
  
  library(reshape2)
  tmp_data <- data.frame(data.frame(t(data), check.names = F), cluster = cluster,
                         check.names = F)
  
  mean_data <- aggregate(.~cluster, tmp_data, mean)
  clustername <- mean_data[,1]
  mean_data <- mean_data[,2:ncol(mean_data)]
  rownames(mean_data) <- clustername; mean_data <- t(mean_data)
  
  return(mean_data)
  
}


average_cs6_hema <- sc_tl_average(cs6_hema@assays$RNA@data, cluster = cs6_hema$CellType_V2)

pheatmap(average_cs6_hema[c('CST1','SOX17','APOE', 'APOC1',
                            'AFP','TTR', 'KDR', 'DCN','HAND1','MDK','COL1A1','MPO','PRTN3','HBE1', 'HBG1', 'GYPA', 'PF4', 'PPBP', 'GP9'), ], 
         cluster_rows = F, cluster_cols = F, scale = 'row')

FeaturePlot(cs6_hema, c('CST1','SOX17','APOE', 'APOC1',
                            'AFP','TTR', 'KDR', 'DCN','HAND1','MDK','COL1A1','MPO','PRTN3','HBE1', 'HBG1', 'GYPA', 'PF4', 'PPBP'), ncol = 6, cols = c('grey', 'red'), order = T)

DimPlot(cs6_hema, reduction = 'spatial', group.by = 'CellType_V2', cols = col_hema)
DimPlot(ana1, reduction = 'spatial', group.by = 'CellType_V2', cols = col.use, raster = F)

########spatial mapping
col_ana1_v2 <- setNames(col.use[1:28], unique(ana1$CellType_V2)) 
DimPlot(ana1[, ana1$slice_num==1], reduction = 'spatial', group.by = 'CellType_V2', cols = col_ana1_v2, raster = F)+
DimPlot(ana1[, ana1$slice_num==5], reduction = 'spatial', group.by = 'CellType_V2', cols = col_ana1_v2, raster = F)+
DimPlot(ana1[, ana1$slice_num==28], reduction = 'spatial', group.by = 'CellType_V2', cols = col_ana1_v2, raster = F)+
DimPlot(ana1[, ana1$slice_num==36], reduction = 'spatial', group.by = 'CellType_V2', cols = col_ana1_v2, raster = F)


pca_average_cs6_hema <- sc_tl_average(t(cs6_hema@reductions$pca@cell.embeddings[,1:20]), cs6_hema$CellType_V2)
pheatmap(cor(pca_average_cs6_hema), clustering_method = 'ward.D2')
