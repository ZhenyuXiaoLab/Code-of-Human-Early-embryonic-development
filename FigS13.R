#code for FigS13
#software version: Seuratwrapper:V0.30; Seurat 4.3.0.1; clusterProfiler 4.2.2; pheatmap 1.0.12; ggplot2 3.5.1
cs7_umap <- readRDS('annot_umap.rds')
cs7_anno <- read.delim('E-MTAB-9388.sdrf.txt', header = T)

n2021 <- CreateSeuratObject(counts = t(cs7_rawdata), meta.data = cs7_umap)
n2021 <- SCTransform(n2021 )
n2021 <- RunPCA(n2021)
n2021 <- RunUMAP(n2021, dims = 1:30)
n2021@reductions$umap@cell.embeddings[] <- as.matrix(cs7_umap[, 2:3])
DimPlot(n2021, group.by = "sub_cluster", label = T)
FeaturePlot(n2021, features = 'ANXA1')

#re-define clusters for hema in CS7
n2021_hema <- n2021[, as.character(n2021$sub_cluster) %in% c('Blood Progenitors', 'Erythroblasts','Hemogenic Endothelium',
                                                             'Erythro-Myeloid Progenitors', 'Myeloid Progenitors')]

n2021_hema <- SCTransform(n2021_hema)
n2021_hema <- RunPCA(n2021_hema)
ElbowPlot(n2021_hema, ndims = 30)
n2021_hema <- RunUMAP(n2021_hema, dims = 1:20)
n2021_hema <- FindNeighbors(n2021_hema, dims = 1:10, k.param = 5)
n2021_hema <- FindClusters(n2021_hema, resolution = seq(0.2, 3, 0.2))

n2021_hema$cluster <- plyr::mapvalues(n2021_hema$SCT_snn_res.2, 0:12,
                                      c('Mac', 'Primitive_Ery', 'YSMP','Primitive_Ery', 'YSMP',
                                        'Primitive_Mk', 'YSMP', 'EC', 'pMP','EC', 'Mes',
                                        'pMP', 'YSMP'))




#####ery lineage comparision between CS6 and CS7, com_cs6_cs7:integrative data between CS6 and CS7 #####
DefaultAssay(com_cs6_cs7) <- 'RNA'
com_cs6_cs7 <- NormalizeData(com_cs6_cs7)

pery <- com_cs6_cs7[, com_cs6_cs7$cluster2=='Primitive_Ery']
pery <- SetIdent(pery, value = 'group')
degs_pery <- FindAllMarkers(pery, logfc.threshold = log(1.25), only.pos = T)
degs_pery <- degs_pery[degs_pery$p_val_adj<0.05, ]


DotPlot(pery, features = sc_tl_topgene(degs_pery, 15), group.by = 'group')+
  coord_flip()+scale_color_gradientn(colours = BlueAndRed(100))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

DotPlot(pery, features = c('HBE1','GYPA', "RHAG", "BPGM", "UROD", "FECH", "ALAS2"),
        group.by = 'group', scale = F)+
  coord_flip()+scale_color_gradientn(colours = colors.use$gradient)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

go_ery_cs6 <- sc_tl_goterms(gene = degs_pery$gene[degs_pery$cluster=='CS6(Pre-Gas)'][1:100],
                            org = 'human', title = 'GO terms of CS6 (Pre-Gas) Ery',
                            col = col_group[1], length = 5, on = 'BP', graph = F)

go_ery_cs6[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))


go_ery_cs7 <- sc_tl_goterms(gene = degs_pery$gene[degs_pery$cluster=='CS7(Peri-Gas)'][1:100],
                            org = 'human', title = 'GO terms of CS7(Peri-Gas) Ery',
                            col = col_group[2], length = 5, on = 'BP', graph = F)

sc_tl_select_go(go_ery_cs7, select_GOID = c("GO:0048821", "GO:0006783", "GO:0030218", "GO:0033014", "GO:0034101"),
                color_use = col_group[2], title = 'GO terms of CS7(Peri-Gas) Ery')

sc_tl_select_go(go_ery_cs7, select_GOID = c("GO:0048821", "GO:0006783", "GO:0030218", "GO:0033014", "GO:0034101"),
                color_use = col_group[2], title = 'GO terms of CS7(Peri-Gas) Ery')


library(dplyr)
ery_genelist <- list('gas_transport' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0015669"),
                     'Ery Differentiation' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0030218"))

pery <- AddModuleScore(pery, features = ery_genelist)
names(pery@meta.data)[64:65] <- names(ery_genelist)

library(ggpubr)

a <- VlnPlot(pery, features = names(ery_genelist)[1], 
        group.by = 'group', cols = col_group, ncol = 1)+stat_compare_means()

b <- VlnPlot(pery, features = names(ery_genelist)[2], 
        group.by = 'group', cols = col_group, ncol = 1)+stat_compare_means()

ggarrange(a, b, ncol = 1, nrow = 2, common.legend = T)

##### mega lineage comparision between CS6 and CS7, com_cs6_cs7:integrative data between CS6 and CS7 ####
pmk <- com_cs6_cs7[, com_cs6_cs7$cluster2=='Primitive_Mk']
pmk <- SetIdent(pmk, value = 'group')

degs_pmk <- FindAllMarkers(pmk, logfc.threshold = log(1.25), only.pos = T)
degs_pmk <- degs_pmk[degs_pmk$p_val_adj<0.05, ]
degs_pmk <- degs_pmk[!degs_pmk$gene %in% rownames(cs6_hema2)[rowMeans(cs6_hema2@assays$RNA@counts)==0],]

DotPlot(pmk, features = sc_tl_topgene(degs_pmk, 15), group.by = 'group')+
  coord_flip()+scale_color_gradientn(colours = BlueAndRed(100))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

DotPlot(pmk, features = c('ITGA2B', 'ITGB3','GP1BA', 'GP6', 'GP9', 'F11R', 'F2RL3', 'F2RL2', 'ACTG1', 'PLEK'), group.by = 'group', scale = F)+
  coord_flip()+scale_color_gradientn(colours =colors.use$gradient)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

go_mk_cs6 <- sc_tl_goterms(gene = degs_pmk$gene[degs_pmk$cluster=='CS6(Pre-Gas)'],
                           org = 'human', title = 'GO terms of CS6(Pre-Gas) Mk',
                           col = col_group[1], length = 5, on = 'BP', graph = F)

go_mk_cs6[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))

go_mk_cs7 <- sc_tl_goterms(gene = degs_pmk$gene[degs_pmk$cluster=='CS7(Peri-Gas)'],
                           org = 'human', title = 'GO terms of CS7(Peri-Gas) mk',
                           col = col_group[2], length = 5, on = 'BP', graph = F)

go_mk_cs7[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))
View(go_mk_cs7[[2]]@result)

sc_tl_select_go(GO_result = go_mk_cs7, select_GOID = c("GO:0030168",
                                                       "GO:0050817",
                                                       "GO:0007596",
                                                       "GO:0030220",
                                                       "GO:0036344"), color_use = col_group[2], title = 'GO terms of CS7(Peri-Gas) Mk')

mk_genelist <- list('coagulation' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0050817"),
                    'MK development' = sc_tl_get_genes_from_GO(org = 'Hsa', GOID = "GO:0035855"))

pmk <- AddModuleScore(pmk, features = mk_genelist)
names(pmk@meta.data)[59:60] <- names(mk_genelist)

p1 <- VlnPlot(pmk, features = names(mk_genelist)[1], group.by = 'group', cols = col_group, ncol = 1)+stat_compare_means()
p2 <- VlnPlot(pmk, features = names(mk_genelist)[2], group.by = 'group', cols = col_group, ncol = 1)+stat_compare_means()

ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = T)

VlnPlot(pmk, features = '', group.by = 'group', cols = col_group, ncol = 1)+geom_boxplot()

####comparision among CS6 pMP, CS7 pMP and YSMP ####
ysmp <- com_cs6_cs7[, com_cs6_cs7$cluster2 %in% c("pMP(CS6)", "pMP(CS7)", "YSMP(CS7)")]
ysmp <- SetIdent(ysmp, value = 'cluster2')
degs_ysmp <- FindAllMarkers(ysmp, logfc.threshold = log(1.25), only.pos = T)
degs_ysmp <- degs_ysmp[degs_ysmp$p_val_adj<0.05, ]
degs_ysmp <- degs_ysmp[grep('RPL|RPS', degs_ysmp$gene, invert = T), ] #filter ribo-related genes


DotPlot(ysmp, group.by = 'cluster', features = rev(sc_tl_topgene(degs_ysmp2, 15)))+
  coord_flip()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = colors.use$BluewhiteRed)

