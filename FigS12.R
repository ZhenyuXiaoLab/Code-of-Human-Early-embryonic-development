#code for Fig.S12
#software version: Seuratwrapper:V0.30; Seurat 4.3.0.1; clusterProfiler 4.2.2; pheatmap 1.0.12; ggplot2 3.5.1
scs_data1 <- readRDS('scs_slice1-7_datalist.rds') #SCS expression 
scs1_loc <- data.frame(scs_data1$reduction$spatial@cell.embeddings[], check.names = F) #location

scs1_loc$slice_num <- scs_data1$meta[rownames(scs1_loc), 'slice_num']
scs_data1 <- CreateSeuratObject(counts = scs_data1$matrix, meta.data = scs_data1$meta)

scs_data2 <- readRDS('scs_ev1-31_obj4.rds')
scs_data2 <- scs_data2@assays$RNA@counts
colnames(scs_data2) <- paste('EV131_', colnames(scs_data2), sep="")
scs_data2 <- CreateSeuratObject(counts = scs_data2)
scs2_loc <- data.frame(readRDS('scs_ev131_trans_spatial.rds'), check.names = F)


scs <- merge(scs_data1[intersect(rownames(scs_data1), rownames(scs_data2)), ],
             scs_data2[intersect(rownames(scs_data1), rownames(scs_data2)), ])

#QC
mt_genes <- grep('^MT-', rownames(scs), value = T)
scs$mt_percent <- colSums(scs@assays$RNA@counts[mt_genes, ])/colSums(scs@assays$RNA@counts)
scs <- scs[, scs$nCount_RNA > 400 & scs$nFeature_RNA>100 & scs$mt_percent <0.05]

#standard pipline
scs <- SCTransform(scs)
scs <- RunPCA(scs, npcs = 50)
ElbowPlot(scs, ndims = 30)
scs <- RunUMAP(scs, dims = 1:10)
scs <- FindNeighbors(scs, dims = 1:10)
scs <- FindClusters(scs, resolution = seq(0.2, 4, 0.2))

DimPlot(scs, group.by = 'SCT_snn_res.2', label = T)
scs$celuster <- paste('c', scs$SCT_snn_res.2, sep="")


#re-cluster for cluster 26, which including Hypo and Epi
c26 <- scs[, scs$SCT_snn_res.2 %in% c(26)]
c26 <- SCTransform(c26)
c26 <- RunPCA(c26, npcs = 50)
ElbowPlot(c26, ndims = 30)
c26 <- RunUMAP(c26, dims = 1:10)
c26 <- FindNeighbors(c26, dims = 1:10)
c26 <- FindClusters(c26, resolution = seq(0.2, 4, 0.2))


c26$cluster <- 'Othercells'
c26$cluster[c26$SCT_snn_res.2 %in% c(0, 7)] <- 'Hypo'
c26$cluster[c26$SCT_snn_res.2 %in% c(4)] <- 'YS.Endo1'
c26$cluster[c26$SCT_snn_res.2 %in% c(8)] <- 'YS.EXMC1'

scs$celuster[colnames(scs) %in% colnames(c26)] <- 
  as.character(na.omit(c26$cluster[match(colnames(scs), colnames(c26))]))


#define clusters
FeaturePlot(scs, features = c('SOX2', 'POU5F1', 'SOX17', 'FOXA2', 'AFP',
                             'GABRP', 'VIM', 'MPO', 'TAL1', 'ALAS2','PF4', 'CGA', 'LEP','COL1A1','LUM','DCN','HLA-C','LAIR2'), ncol = 5)

DimPlot(scs, group.by = 'celuster', label = T)

scs$celuster[scs$celuster %in% c('c1', 'c8', 'c14', 'c15', 'c25')] <- 'PL.EXMC'
scs$celuster[scs$celuster %in% c('c0', 'c2', 'c3', 'c4','c7','c12','c13','c19','c20','c21', 'c25','c29')] <- 'CTB/STB'
scs$celuster[scs$celuster %in% c('c5','c9','c11', 'c16', 'c18','c23','c30','c34','c35')] <- 'MTB'
scs$celuster[scs$celuster %in% c('c33')] <- 'Anterior.pole'
scs$celuster[scs$celuster %in% c('c6','c22')] <- 'Epi'
scs$celuster[scs$celuster %in% c('c10')] <- 'AM/AM.Ecto/AM.EXMC'
scs$celuster[scs$celuster %in% c('c17')] <- 'YS.EXMC2'
scs$celuster[scs$celuster %in% c('c27','c28')] <- 'YS.Endo'
scs$celuster[scs$celuster %in% c('c24','c31','c32')] <- 'Blood'

scs2 <- scs[, !scs$celuster=='Othercells'] #filter low quality cells
scs2$celuster[scs2$celuster %in% c('YS.EXMC1', 'YS.EXMC2')] <- 'YS.EXMC'
scs2$celuster[scs2$celuster %in% c('YS.Endo1', 'YS.Endo')] <- 'YS.Endo'


scs_cluster_order <- c('Anterior.pole', 'Epi', 'Hypo', 'AM/AM.Ecto/AM.EXMC', 'PL.EXMC', 'CTB/STB', 'MTB', 'YS.Endo', 'YS.EXMC', 'Blood')
scs2$celuster <- factor(scs2$celuster, levels = scs_cluster_order)

col_scs_cluster <- setNames(rgb(c(232, 220, 143, 206, 23, 46, 197, 183, 156),
                       c(108, 111, 54, 43, 63, 129, 224, 186, 38),
                       c(121, 42, 31, 42, 116, 171, 193, 214, 116), maxColorValue = 255), 
                       c('Anterior.pole', 'YS.Endo', 'YS.EXMC', 'Blood', 'MTB', 'CTB/STB', 'PL.EXMC', 'AM/AM.Ecto/AM.EXMC', 'Hypo'))

col_scs_cluster <- c(col_scs_cluster, 'Epi' = 'darkgreen')[scs_cluster_order]
col_scs_cluster['Hypo'] <- 'magenta3'

DimPlot(scs2, group.by = 'celuster', label = F, label.size = 4, repel = F, cols = col_scs_cluster, raster = F)


genes2 <- 
  c("GSC", "NODAL", "FGF4", "FOXA2", "CHRD",
             "NOTO", 
             "CER1",
             "POU5F1",
             "UCHL1",
             "SOX2",
             "OTX2",
    "SOX17",
    "FOXA3",
    "GATA4",
             "GABRP",
             "ISL1",
             "MYL7",
             "DLK1",
             "VIM",
             "ERVW-1", 
             "CGA", 
             "SDC1",
             "HLA-G",
             "LAIR2",
             "TTR",
             "AFP",
             "TF",
    'KDR','DCN','POSTN',
             "PF4",
             "PPBP",
             "GYPC")

DefaultAssay(scs2) <- 'RNA'
scs2 <- NormalizeData(scs2)
scs2 <- ScaleData(scs2)

scs2 <- SetIdent(scs2, value = 'celuster')
scs2_average <- AverageExpression(scs2, slot = 'data', assays = 'RNA')

DotPlot(scs2, features = rev(genes2), group.by = 'celuster', scale = T)+coord_flip()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_color_gradientn(colours = BlueAndRed(100))

pheatmap(t(scale(t(scs2_average$RNA)))[genes2, ], cluster_rows = F, cluster_cols = F, color = colors.use$BluewhiteRed)

#foucus on yolk sac
ys <- scs[, scs$celuster %in% c('Hypo', 'YS.Endo1', 'YS.Endo','YS.EXMC1','YS.EXMC2','Blood')]
ys <- SCTransform(ys)
ys <- RunPCA(ys, npcs = 50)
ElbowPlot(ys, ndims = 30)

ys <- RunUMAP(ys, dims = 1:10, min.dist = 0.5, spread = 1, local.connectivity = 2, n.neighbors = 100)
ys <- FindNeighbors(ys, dims = 1:10)
ys <- FindClusters(ys, resolution = seq(0.2, 4, 0.2), k=10)
ys$celuster <- scs[, colnames(ys)]$celuster

ys$cluster <- as.character(ys$SCT_snn_res.4)
ys$cluster[ys$SCT_snn_res.4 %in% c(21,26)] <- 'Hypo'
ys$cluster[ys$SCT_snn_res.4 %in% c(19, 29)] <- 'YS.Endo1'
ys$cluster[ys$SCT_snn_res.4 %in% c(10)] <- 'YS.EXMC1'
ys$cluster[ys$SCT_snn_res.4 %in% c(17, 27, 31)] <- 'MTB'
ys$cluster[ys$SCT_snn_res.4 %in% c(0, 1, 8)] <- 'YS.EXMC2'
ys$cluster[ys$SCT_snn_res.4 %in% c(18, 30, 22, 4, 18)] <- 'YS.Endo2'
ys$cluster[ys$SCT_snn_res.4 %in% c(2, 3, 5, 6, 20, 28)] <- 'YS.Endo3'
ys$cluster[ys$SCT_snn_res.4 %in% c(11, 12, 25, 32)] <- 'pMeg'
ys$cluster[ys$SCT_snn_res.4 %in% c(7, 16, 23, 24, 14, 9, 13)] <- 'pEry'

#define pMP which locate in YS.EXMC1 and YS.Endo1
ys$tmp <- ys$SCT_snn_res.4
ysexmc1 <- ys[, ys$cluster %in% c('YS.EXMC1', 'YS.Endo1')]
ysexmc1 <- SCTransform(ysexmc1)
ysexmc1 <- RunPCA(ysexmc1, npcs = 50)
ElbowPlot(ysexmc1, ndims = 30)
ysexmc1 <- RunUMAP(ysexmc1, dims = 1:5)
ysexmc1 <- FindNeighbors(ysexmc1, dims = 1:5, k.param = 5)
ysexmc1 <- FindClusters(ysexmc1, resolution = seq(0.2, 4, 0.2))

FeaturePlot(ysexmc1, features = c('MPO', 'PRTN3','KDR','VIM','AFP','FOXA2'), ncol = 3)

ysexmc1$cluster[ysexmc1$SCT_snn_res.2 %in% c(2, 3, 7, 10)] <-'YS.EXMC1'
ysexmc1$cluster[ysexmc1$SCT_snn_res.2 %in% c(1, 4, 8, 12)] <-'YS.Endo1'
ysexmc1$cluster[ysexmc1$SCT_snn_res.2 %in% c(0, 5, 6, 9, 11)] <-'pMP'
ysexmc1$cluster[ysexmc1$SCT_snn_res.2 %in% c(14)] <-'Hypo'

ys$cluster[colnames(ys) %in% colnames(ysexmc1)] <-
  as.character(na.omit(ysexmc1$cluster[match(colnames(ys), colnames(ysexmc1))]))

ys2 <- ys[, !ys$cluster %in% c(15, 'MTB')] #filter cluster 15:low quality

ys2_cluster_order <- c('Hypo', "YS.Endo1",  "YS.Endo2",  "YS.Endo3",  "YS.EXMC1",  "YS.EXMC2",
                       "pMP", "pEry", "pMeg")



save(file = 'scs_250306.Rdata', scs, ys2, ys3)

######PAGA analysis ####
load('/data/gyd/project/CS6/review/scs_250306.Rdata') #SCS data with annotated cell clusters


ys3  <- scs[,c(colnames(scs)[scs$celuster %in% c('Anterior.pole', 'Epi')], 
               colnames(ys2))]

ys3$cluster <- ys3$celuster
ys3$cluster[colnames(ys3) %in% colnames(ys2)] <-
  as.character(na.omit(ys2$cluster[match(colnames(ys3), colnames(ys2))]))


ys3 <- NormalizeData(ys3)
ys3 <- FindVariableFeatures(ys3)
ys3 <- ScaleData(ys3)
ys3 <- RunPCA(ys3, npcs = 50)
ElbowPlot(ys3, ndims = 30)

ys3 <- RunUMAP(ys3, dims = 1:10, n.neighbors = 50,
               local.connectivity = 1, min.dist = 0.5, spread = 1)

col_ys3 <- c(col_ys2, col_scs_cluster['Anterior.pole'])

ys3$cluster <- factor(ys3$cluster, levels = c('Anterior.pole', 'Epi', ys2_cluster_order))

ys3$Cluster <- as.character(ys3$cluster)
ys3$Cluster[ys3$Cluster=='Anterior.pole'] <- 'Epi'
ys3$Cluster <- factor(ys3$Cluster, levels = rev(c('Epi', ys2_cluster_order)))


DimPlot(ys3, group.by = 'Cluster', label = F, cols = col_ys3, label.size = 4, repel = T)

DotPlot(ys3, group.by = 'Cluster', features = c('SOX2','SOX17', 'CST1', 'AFP', 'FGA', 'KDR', 'ANXA1',
                                                'MPO', 'PRTN3', 'GYPA', 'HBE1', 'PF4', 'GP9'), scale.by = 'size')+
  scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6','orange', 'red3','brown'))(100))+theme(axis.text.x = element_text(angle = 60, hjust = 1))


ys3_paga <- RunPAGA(ys3, group_by = 'cluster', nonlinear_reduction = 'umap')


##### expression of ligand and receptor genes ####
a <- DotPlot(ys2[, ys2$cluster %in% c('Hypo','YS.Endo1','YS.EXMC1','pMP')], group.by = 'cluster', features = rev(c('SMO', 'PTCH1', 'IHH')), scale.by = 'size')+coord_flip()+
  scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6','orange', 'red3','brown'))(100))+theme(axis.text.x = element_text(angle = 60, hjust = 1))


b <- DotPlot(ys2[, ys2$cluster %in% c('YS.Endo2','YS.Endo3','YS.EXMC2','pEry','pMeg')], group.by = 'cluster', features = rev(c('EPO', 'EPOR')), scale.by = 'size')+coord_flip()+
  scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6','orange', 'red3','brown'))(100))+theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggarrange(a, b, ncol = 2)
