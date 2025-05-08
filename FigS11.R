#code for FigS11

#### perform developmental trajectory analysis using clusters located in yolk sac ####
#slingshot analysis
cs6_sub_slingshot <- RunSlingshot(cs6_hema2[, !cs6_hema2$cluster_final %in% c('CS', 'pEry2', 'pMeg2')],
                                  group.by = 'cluster_final', reduction = 'pca', dims = 1:10,  
                                  align_start = F, end = c( 'pMP', 'pEry1', 'pMeg1','YS.Endo3'), start = c('Hypo.1'))

cs6_sub_slingshot$cluster_final <- factor(cs6_sub_slingshot$cluster_final, levels = names(all_colors$cs6_hema_cluster_v1)[1:9])

CellDimPlot(cs6_sub_slingshot, group.by = "cluster_final", reduction = "umap",
            lineages = paste0("Lineage", 1:4), lineages_span = 0.5, lineages_trim = c(0.01, 0.99),
            palcolor = all_colors$cs6_hema_cluster_v1)

plot(cs6_sub_slingshot@tools$Slingshot_cluster_final_pca@metadata$lineages, col = col_cluster2[cs6_sub_slingshot@meta.data$cluster_final], pch=16, asp=1)
plot(cs6_sub_slingshot@tools$Slingshot_cluster_final_pca@metadata$mst)

#monocle3 analysis
library(SeuratWrappers)
library(monocle3)


cs6_sub_monocle3 <- RunMonocle3(cs6_hema2[, !cs6_hema2$cluster_final %in% c('CS', 'pEry2', 'pMeg2')],
                                  clusters = 'cluster_final', reduction = 'umap', k = 100,
                                close_loop = F, root_cells = colnames(cs6_hema2)[cs6_hema2$cluster_final=='Hypo.1'], resolution = 0.1)

trajectory <- cs6_sub_monocle3@tools$Monocle3$trajectory
CellDimPlot(cs6_sub_monocle3 , group.by = "cluster_final", reduction = "umap",
            label = TRUE, theme_use = "theme_blank") + trajectory



####DEGs between 2 EXMC  clusters #####
ys_exmc <- cs6_hema2[, cs6_hema2$CellType_V2 %in% c('YS.EXMC_1', 'YS.EXMC_2')]
ys_exmc <- SetIdent(ys_exmc, value = 'CellType_V2')

degs_ys_exmc <- FindAllMarkers(ys_exmc, logfc.threshold = log(1.25), only.pos = T)
degs_ys_exmc <- degs_ys_exmc[degs_ys_exmc$p_val_adj<0.05, ]

DefaultAssay(ys_exmc) <- 'RNA'
ys_exmc <- ScaleData(ys_exmc)
average_ys_exmc <- sc_tl_average(as.matrix(ys_exmc@assays$RNA@data), ys_exmc$CellType_V2)

DotPlot(ys_exmc, features = rev(sc_tl_topgene(degs_ys_exmc, 20)), scale.by = 'size', scale = F)+
  coord_flip()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_color_gradientn(colours = c(colorRampPalette(c('#E6E6E6','blue', 'yellow'))(15), colorRampPalette(c('yellow','red', 'brown'))(80)))
  
##### FateID prediction #####
library(FateID)
data.use <- cs6_hema3@assays$RNA@data

cs6_hema3$cluster.use <- as.character(cs6_hema3$CellType_V2)
cs6_hema3$cluster.use[cs6_hema3$cluster.use %in% c('Primitive_Ery1', 'Primitive_Mk1')] <- 'Ery_Meg'

fb  <- fateBias(as.matrix(cs6_hema3@assays$RNA@data[unique(degs_hema3$gene[degs_hema3$cluster %in% c('Myeloid Progenitor','Primitive_Ery1', 'Primitive_Mk1')]), ]) , cs6_hema3$cluster.use,
                c('Myeloid Progenitor', 'Ery_Meg'), z=NULL, minnr=5, minnrh=10, adapt=TRUE,
                confidence=0.75, nbfactor=5, use.dist=FALSE, seed=12345, nbtree=NULL)

fb$votes <- fb$votes/500 
fb$votes$cluster <- cs6_hema3@meta.data[rownames(fb$votes), 'CellType_V2']

fb$votes$cluster <- as.character(fb$votes$cluster)
fb$votes$cluster[fb$votes$cluster=='YS.EXMC_1'] <- 'YS.EXMC1'
fb$votes$cluster[fb$votes$cluster=='YS.EXMC_2'] <- 'YS.EXMC2'

a <- fb$votes %>%  filter(cluster %in% c("YS.EXMC1", "YS.EXMC2")) %>%  ggplot(aes(x = cluster, y= `tMyeloid Progenitor`, fill=cluster))+geom_boxplot(outlier.size=0)+
  stat_compare_means()+theme_bw()+scale_fill_manual(values = all_colors$cs6_hema_cluster_v1)

b <- fb$votes %>%  filter(cluster %in% c("YS.EXMC1", "YS.EXMC2")) %>%  ggplot(aes(x = cluster, y= `tEry_Meg`, fill=cluster))+geom_boxplot(outlier.size=0)+
  stat_compare_means()+theme_bw()+scale_fill_manual(values = all_colors$cs6_hema_cluster_v1)

ggarrange(a,b, ncol = 2, common.legend = T, legend = 'right')
