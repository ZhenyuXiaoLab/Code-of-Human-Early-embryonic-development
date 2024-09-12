#######subset hematopoietic related clusters #############
###loading cs6_hema2
cs6_hema2 <- NormalizeData(cs6_hema2)
cs6_hema2 <- FindVariableFeatures(cs6_hema2)
cs6_hema2 <- ScaleData(cs6_hema2)
cs6_hema2 <- RunPCA(cs6_hema2)

ElbowPlot(cs6_hema2, ndims = 30)
cs6_hema2 <- RunUMAP(cs6_hema2, dims = 1:20, min.dist = 0.3, spread = 0.5)

col_hema <- setNames(col.use[1:12], c('Hypo1', 'YS.Endo_1','YS.Endo_2','YS.Endo_3', 'YS.EXMC_1','YS.EXMC_2', 'EXMC_Prog', "pMP",
                                      'pEry1', 'pEry2', 'pMega1','pMega2'))

cs6_hema2$cluster <- factor(cs6_hema2$cluster, levels = c('Hypo1', 'YS.Endo_1','YS.Endo_2','YS.Endo_3', 'YS.EXMC_1','YS.EXMC_2', 'EXMC_Prog', "pMP",
                                      'pEry1', 'pEry2', 'pMega1','pMega2'))

DimPlot(cs6_hema2, group.by = 'cluster', cols = col.use, label = T, repel = T, pt.size = 1)

######Feature Genes 
cs6_hema2 <- SetIdent(cs6_hema2, value = 'cluster')
degs_cs6_hema <- FindAllMarkers(cs6_hema2, logfc.threshold = log(1.25), only.pos = T, max.cells.per.ident = 1000)
degs_cs6_hema <- degs_cs6_hema[degs_cs6_hema$p_val_adj<0.05, ]

DotPlot(cs6_hema2, features = rev(sc_tl_topgene(degs_cs6_hema, 5)),
        group.by = 'cluster')+coord_flip()+
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


average_cs6_hema <- sc_tl_average(cs6_hema2@assays$RNA@data, cluster = cs6_hema2$cluster)

pheatmap(average_cs6_hema[c('CST1','SOX17','APOE', 'APOC1',
                            'AFP','TTR', 'KDR', 'DCN','HAND1','MDK','COL1A1','MPO','PRTN3','HBE1', 'HBG1', 'GYPA', 'PF4', 'PPBP', 'GP9'), ], 
         cluster_rows = F, cluster_cols = F, scale = 'row')

FeaturePlot(cs6_hema2, c('CST1','SOX17','APOE', 'APOC1',
                            'AFP','TTR', 'KDR', 'DCN','HAND1','MDK','COL1A1','MPO','PRTN3','HBE1', 'HBG1', 'GYPA', 'PF4', 'PPBP'), ncol = 6, cols = c('grey', 'red'), order = T)

DimPlot(cs6_hema2, reduction = 'spatial', group.by = 'cluster', cols = col_hema)
DimPlot(ana1, reduction = 'spatial', group.by = 'cluster', cols = col.use, raster = F)

########spatial mapping
col_ana1_v2 <- setNames(col.use[1:28], unique(ana1$cluster)) 
DimPlot(ana1[, ana1$slice_num==1], reduction = 'spatial', group.by = 'cluster', cols = col_ana1_v2, raster = F)+
DimPlot(ana1[, ana1$slice_num==5], reduction = 'spatial', group.by = 'cluster', cols = col_ana1_v2, raster = F)+
DimPlot(ana1[, ana1$slice_num==28], reduction = 'spatial', group.by = 'cluster', cols = col_ana1_v2, raster = F)+
DimPlot(ana1[, ana1$slice_num==36], reduction = 'spatial', group.by = 'cluster', cols = col_ana1_v2, raster = F)


pca_average_cs6_hema <- sc_tl_average(t(cs6_hema2@reductions$pca@cell.embeddings[,1:20]), cs6_hema2$cluster)
pheatmap(cor(pca_average_cs6_hema), clustering_method = 'ward.D2')
