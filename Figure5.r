####### code for Figure 5 ########
#######correlation analysis
pca_average_cs6_hema <- sc_tl_average(t(cs6_hema2@reductions$pca@cell.embeddings[,1:15]), cs6_hema2$CellType_V2)

pheatmap::pheatmap(cor(pca_average_cs6_hema), cluster_rows = T, cluster_cols = T, clustering_method = 'ward.D2',
         display_numbers = F, color = colors.use$BluewhiteRed, treeheight_row = 10, treeheight_col = 10,
         fontsize_number = 8, border_color = 'grey')

########spatial distance
average_dist <- sc_tl_average(t(cs6_hema2_3d_coord[, 1:3]), cs6_hema2_3d_coord$CellType)
heatmap_order3 <- c("Hypo.2", "YS.Endo_1", "YS.EXMC_1" ,"Myeloid Progenitor",
                    "YS.Endo_2", "YS.Endo_3",        
                    "YS.EXMC_2", "Primitive_Ery1", "Primitive_Mk1")

heatmap_order4 <- c("Hypo.2", "YS.Endo_1", "YS.EXMC_1" ,
                    "YS.Endo_2", "YS.Endo_3",        
                    "YS.EXMC_2")

pheatmap::pheatmap(as.matrix(dist(t(average_dist), upper = F))[heatmap_order3, heatmap_order3], 
                   color = rev(colors.use$cm), 
                   border_color = 'grey', treeheight_row = 10, treeheight_col = 10)

##########position between EXMC and Endo
####transcriptomic similarity between EXMC and Endo
exmc_endo <- cs6_hema2[, cs6_hema2$CellType_V2 %in% heatmap_order4]
average_exmc_endo <- sc_tl_average(t(exmc_endo@reductions$pca@cell.embeddings[, 1:10]), exmc_endo$CellType_V2)

pheatmap::pheatmap(cor(average_exmc_endo), 
                   color = colors.use$BluewhiteRed, 
                   border_color = 'grey', treeheight_row = 10, treeheight_col = 10, display_numbers = F)

#####spatial mapping
a <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==7], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)

b <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==8], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)

c <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==9], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)

d <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==10], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)

e <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==11], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)

f <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==12], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)

ggarrange(a, b, c, d, e, f, ncol = 3, nrow = 2, common.legend = T, legend = 'right')

#####
average_dist2 <- sc_tl_average(t(cs6_hema2_3d_coord[colnames(cs6_hema2)[cs6_hema2$slice_num %in% 7:12], 1:3]),
                               cs6_hema2_3d_coord[colnames(cs6_hema2)[cs6_hema2$slice_num %in% 7:12],]$CellType)

pheatmap::pheatmap(as.matrix(dist(t(average_dist2), upper = F))[heatmap_order4, heatmap_order4], 
                   color = rev(colors.use$BluewhiteRed), 
                   border_color = 'grey', treeheight_row = 10, treeheight_col = 10, display_numbers = F)

#####
VlnPlot(cs6_endo2[, cs6_endo2$CellType_V2 %in% c('YS.Endo_1', 'YS.Endo_2', 'YS.Endo_3')], group.by = 'CellType_V2',
        features = c('ID1','FGFR1','VCAN','MDK','KRT8', 'FGA', 'FGB', 'AFP', "TTR", 'GPC3'), cols = all_colors$cs6_hema_cluster_v1, pt.size = 0, ncol = 5)

#####YE maration
DotPlot(cs6_endo2[, cs6_endo2$CellType_V2 %in% c('YS.Endo_1', 'YS.Endo_2', 'YS.Endo_3')], group.by = 'CellType_V2', scale.by = 'size',
        features = rev(c('FOXA2', 'SOX17','DPYS','AKR1D1', 'VTN', 'TF', "GJB1", 'APOA2')))+coord_flip()+
   theme(axis.text.x = element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = colors.use$gradient)


##########DEG between 2 EXMCs#####
ys_exmc <- cs6_hema2[, cs6_hema2$CellType_V2 %in% c('YS.EXMC_1', 'YS.EXMC_2')]
ys_exmc <- SetIdent(ys_exmc, value = 'CellType_V2')

degs_ys_exmc <- FindAllMarkers(ys_exmc, logfc.threshold = log(1.25), only.pos = T)
degs_ys_exmc <- degs_ys_exmc[degs_ys_exmc$p_val_adj<0.05, ]

DefaultAssay(ys_exmc) <- 'RNA'
ys_exmc <- ScaleData(ys_exmc)
average_ys_exmc <- sc_tl_average(as.matrix(ys_exmc@assays$RNA@data), ys_exmc$CellType_V2)
average_ys_exmc2 <- sc_tl_average(as.matrix(ys_exmc@assays$RNA@scale.data), ys_exmc$CellType_V2)

DotPlot(ys_exmc, features = rev(sc_tl_topgene(degs_ys_exmc, 20)), scale.by = 'size', scale = F)+
  coord_flip()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_color_gradientn(colours = c(colorRampPalette(c('#E6E6E6','blue', 'yellow'))(15), colorRampPalette(c('yellow','red', 'brown'))(80)))
  

pheatmap::pheatmap(t(scale(t(average_ys_exmc[sc_tl_topgene(degs_ys_exmc, 20), ]))), cluster_rows = F,
                   cluster_cols = F, color = colors.use$BluewhiteRed)

pheatmap::pheatmap(cs6_hema2_average2[intersect(sc_tl_topgene(degs_ys_exmc, 20), rownames(average_ys_exmc2)), c('YS.EXMC_1', 'YS.EXMC_2')], cluster_rows = F,
                   cluster_cols = F, color = colors.use$BluewhiteRed)


ys_exmc$cluster <- as.character(ys_exmc$CellType_V2)

ysexmc_cellcycle <- sc_tl_cell_cycle_analyze(data_exprMat = ys_exmc@assays$RNA@data, anno = ys_exmc@meta.data,
                                         col = all_colors$cs6_hema_cluster_v1[c("YS.EXMC_1", "YS.EXMC_2")], thre = 0.2, org = 'human',
                                         label = c("YS.EXMC_1", "YS.EXMC_2"))

######hypixia gene  score
library(AnnotationDbi)
o2_genelist <- list('oxi_stress' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0006979"),
                    'oxi_level' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0036293"),
                    'hypoxia' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0001666"))

cs6_hema2 <- AddModuleScore(cs6_hema2, features = o2_genelist, name = names(o2_genelist))
names(cs6_hema2@meta.data)[98:100] <- names(o2_genelist)

load("/data2/public/dong_all_colors.RData")
cs6_hema3 <- cs6_hema2[, !cs6_hema2$CellType_V2 %in% c('Connecting Stalk', 'Primitive_Ery2', 'Primitive_Mk2')]

my_comparisons <- list( c("YS.Endo_1", "YS.Endo_2"), c("YS.EXMC_1", "YS.EXMC_2"))

ggviolin(cs6_hema3@meta.data, x = 'CellType_V2', y= "oxi_level", fill = 'CellType_V2')+
      geom_hline(yintercept = mean(cs6_hema3$oxi_level), lty=2)+geom_boxplot(width=0.3)+scale_fill_manual(values = all_colors$cs6_hema_cluster_v1)+
       stat_compare_means(method = 't', comparisons = list(c(3, 4), c(5,6)))+

ggviolin(cs6_hema3@meta.data, x = 'CellType_V2', y= "hypoxia", fill = 'CellType_V2')+
      geom_hline(yintercept = mean(cs6_hema3$hypoxia), lty=2)+geom_boxplot(width=0.3)+scale_fill_manual(values = all_colors$cs6_hema_cluster_v1)+
       stat_compare_means(method = 't', comparisons = list(c(3, 4), c(5,6)))

a <- VlnPlot(cs6_hema3, group.by = 'CellType_V2', cols = all_colors$cs6_hema_cluster_v1,
        features = names(o2_genelist)[2], pt.size = 0)+geom_boxplot(width = 0.3, outlier.size = 0)+
      geom_hline(yintercept = mean(cs6_hema3$oxi_level), lty=2)+stat_compare_means(method = 't', comparisons = my_comparisons, label.y = c(0.2, 0.2))

b <- VlnPlot(cs6_hema3, group.by = 'CellType_V2', cols = all_colors$cs6_hema_cluster_v1,
        features = names(o2_genelist)[3], pt.size = 0)+geom_boxplot(width = 0.3, outlier.size = 0)+
        geom_hline(yintercept = mean(cs6_hema3$hypoxia), lty=2)+stat_compare_means(comparisons = my_comparisons)

ggarrange(a, b, common.legend = T, align = 'hv', legend = 'right')

######
a <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==7], group.by = 'CellType_V2',
        reduction = 'spatial', cols = all_colors$cs6_hema_cluster_v1)+labs(title = 'slice 7')

b <- FeaturePlot(cs6_hema2[, cs6_hema2$slice_num==7], features = names(o2_genelist)[2], reduction = 'spatial')+
  scale_color_gradientn(colours = colors.use$BluewhiteRed)

c <- FeaturePlot(cs6_hema2[, cs6_hema2$slice_num==7], features = names(o2_genelist)[3], reduction = 'spatial')+
  scale_color_gradientn(colours = colors.use$BluewhiteRed)

ggarrange(a, b, c, ncol = 3, nrow = 1, align = 'hv')


p_5 <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==7], group.by = 'CellType_V2',
        reduction = 'spatial', cols = hema_cluster_col)

p_6 <- DimPlot(cs6_hema2[, cs6_hema2$slice_num==6], group.by = 'CellType_V2',
        reduction = 'spatial', cols = hema_cluster_col)

ggarrange(p_5, p_6, ncol = 1, nrow = 2, align = 'hv', common.legend = T, legend = 'right')


############lineange gene expression#####
cs6_hema2_average <- sc_tl_average(cs6_hema2@assays$RNA@data, cs6_hema2$CellType_V2)

cs6_hema2 <- ScaleData(cs6_hema2)
cs6_hema2_average2 <- sc_tl_average(cs6_hema2@assays$RNA@scale.data, cs6_hema2$CellType_V2)
range(cs6_hema2_average)

genes.use <- c('ANXA1','POSTN','KDR','GATA2',
               'LMO2', 'TAL1', 'LYL1', 'FLI1', 'MPO', 'PRTN3', 'LYZ','APOA1','TF','AFP','TTR','HBE1',"HBZ", 'PF4','PLEK')

pheatmap::pheatmap(cs6_hema2_average[intersect(c('KDR','MPO', 'PRTN3','SRGN',
                                                 'LMO2', 'TAL1', 'LYL1', 'GATA2', 'NFE2', 'HBZ','HBG','HBE1', 'PF4','PLEK'),
                            rownames(cs6_hema2)), c('YS.EXMC_1', 'YS.EXMC_2')], cluster_rows = F, cluster_cols = F)


pheatmap::pheatmap(t(cs6_hema2_average[genes.use, c('YS.EXMC_1', 'YS.EXMC_2')]), cluster_rows = F, cluster_cols = F,
                   
                            color = c(colorRampPalette(c('#E6E6E6', 'yellow'))(80), colorRampPalette(c('orange', 'red', 'brown'))(2000)))

