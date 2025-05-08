#code for Fig.S4
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

