#code for Fig.S10
#software version: Seuratwrapper:V0.30; Seurat 4.3.0.1; clusterProfiler 4.2.2; pheatmap 1.0.12; ggplot2 3.5.1;
######cell origin of Pery2 and Pmega2#######
#oad Ery and Mk data from adult bone marrow and embryonic yolk sac
load('bm_ys_ery_mk.Rdata')

#calculate DEGs for ery/mega between embryoic and adult stages 
ery_mk_bm_ys$site[ery_mk_bm_ys$site=='BM'] <- 'Adult'
ery_mk_bm_ys$site[ery_mk_bm_ys$site=='YS'] <- 'Embryonic'

all_ery <- ery_mk_bm_ys[, ery_mk_bm_ys$celltype=='Ery']
all_ery <- SetIdent(all_ery, value = 'site')
degs_all_ery <- FindAllMarkers(all_ery, logfc.threshold = log(1.25), only.pos = T)
degs_all_ery <- degs_all_ery[degs_all_ery$p_val_adj<0.05, ]

degs_all_ery <- degs_all_ery[grep('^IG', degs_all_ery$gene, invert = T), ]#filter immunoglobulin genes
degs_all_ery <- degs_all_ery[order(degs_all_ery$avg_log2FC, decreasing = T), ]
degs_all_ery <- degs_all_ery[degs_all_ery$pct.1 > 0.2,  ]#obtain genes with high expressed proportion

a <- DotPlot(all_ery, features = rev(sc_tl_topgene(degs_all_ery, 20)), group.by = 'site', scale = T)+
  coord_flip()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+labs(x="", y="")+
  scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))

all_meg <- ery_mk_bm_ys[, ery_mk_bm_ys$celltype=='Mega']
all_meg <- SetIdent(all_meg, value = 'site')
degs_all_meg <- FindAllMarkers(all_meg, logfc.threshold = log(1.25), only.pos = T)
degs_all_meg <- degs_all_meg[degs_all_meg$p_val_adj<0.05, ]
degs_all_meg <- degs_all_meg[grep('^IG', degs_all_meg$gene, invert = T), ]
degs_all_meg <- degs_all_meg[order(degs_all_meg$avg_log2FC, decreasing = T),]
degs_all_meg <- degs_all_meg[degs_all_meg$pct.1>0.2, ]

b <- DotPlot(all_meg, features = rev(sc_tl_topgene(degs_all_meg, 20)), group.by = 'site', scale = T)+
  coord_flip()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+labs(x="", y="")+
  scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))

#check expression in ery and mega cluster in CS6 
cs6_hema2$cluster_final <- plyr::mapvalues(as.character(cs6_hema2$CellType_V2), c("Hypo.2", "YS.Endo_1", "YS.Endo_2", "YS.Endo_3",
                                                                                  "YS.EXMC_1", "YS.EXMC_2", "Myeloid Progenitor",
                                                                                  "Primitive_Ery1", "Primitive_Ery2", "Primitive_Mk1",
                                                                                  "Primitive_Mk2", "Connecting Stalk"), c('Hypo1', 'YS.Endo1', 'YS.Endo2', 'YS.Endo3',
                                                                                                                          'YS.EXMC1', 'YS.EXMC2', 'pMP', 'pEry1', 'pEry2', 'pMeg1', 'pMeg2', 'CS'))


c <- DotPlot(cs6_hema2[, cs6_hema2$cluster_final %in% c('pEry1', 'pEry2')], group.by = 'cluster_final',
        features = rev(sc_tl_topgene(degs_all_ery,20)), scale = F)+coord_flip()+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))


d <- DotPlot(cs6_hema2[, cs6_hema2$cluster_final %in% c('pMeg1', 'pMeg2')], group.by = 'cluster_final',
        features = rev(sc_tl_topgene(degs_all_meg,20)), scale = F)+coord_flip()+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))



ggarrange(a, b, ncol = 2, nrow =1)
ggarrange(c, d, ncol = 2, nrow =1)

#### integratative analysis for ery and mega lineage #####
load('cs10_cs11_cs15_integrated.Rdata')

#### integratative analysis for ery lineage #####
ery <- ys.integrated[, as.character(ys.integrated$Cluster) %in% c('ErP', 'pEry')]
ery <- NormalizeData(ery)
ery <- ScaleData(ery)
ery <- FindVariableFeatures(ery)
hvg_ery <- VariableFeatures(ery)

noise_ery <- sc_tl_detect_noise_genes(ery)
ery <- RunPCA(ery, features = setdiff(VariableFeatures(ery), noise_ery))


library(SeuratWrappers)
ery@assays$RNA@scale.data <- as.matrix(0)

ery <- RunFastMNN(SplitObject(ery, split.by = 'stage'),
                  features = setdiff(hvg_ery, noise_ery))
ery <- RunUMAP(ery, dims = 1:30, reduction = 'mnn')
ery <- FindNeighbors(ery, dims = 1:30, reduction = 'mnn')
ery <- FindClusters(ery, resolution = seq(0.2, 2, 0.2))

#Define clusters
ery$celltype <- plyr::mapvalues(ery$RNA_snn_res.0.4, 0:4, c('Late_Ery', 'Late_Ery', 'Late_Ery', 'Early_Ery', 'Erythroblast'))
col_ery <- setNames(c('#F1919C', '#Cf90EC', '#32a777'),  c('Erythroblast', 'Early_Ery', 'Late_Ery'))
DimPlot(ery, group.by = 'celltype', label = T, cols = col_ery)
ery$celltype <- factor(ery$celltype, levels = c('Erythroblast', 'Early_Ery', 'Late_Ery'))

#DEGs
ery <- SetIdent(ery, value = 'celltype')
degs_ery <- FindAllMarkers(ery, logfc.threshold = log(1.25), only.pos = T)
degs_ery <- degs_ery[degs_ery$p_val_adj<0.05, ]

average_ery <- AverageExpression(ery)
degs_ery2 <- sc_tl_filter_deg(degs_ery, average_data = average_ery$RNA)

DotPlot(ery, features = sc_tl_topgene(degs_ery, 10), group.by = 'celltype')+
  coord_flip()+theme(axis.text.x =element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = BlueAndRed(100))

DotPlot(ery, features = sc_tl_topgene(degs_ery2, 10), group.by = 'celltype')+
  coord_flip()+theme(axis.text.x =element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = BlueAndRed(100))

#GO analysis
sc_tl_goterms <- function(gene, org, title=NULL, col, length, on="BP", graph=FALSE){
  
  library(clusterProfiler)
  library(ggplot2)
  
  if(org=="mouse"){
    library(org.Mm.eg.db)
    entrez_gene<-bitr(gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
    enrichGOA <- enrichGO(gene = entrez_gene$ENTREZID, OrgDb = org.Mm.eg.db,keyType ="ENTREZID", ont = on,  pAdjustMethod = "BH",pvalueCutoff  = 1,qvalueCutoff  = 1,readable= TRUE)
    godata<-enrichGOA@result
    godata<-godata[order(godata$pvalue,decreasing = F),]
    
    godata_filter<-godata[1:length,]
    godata_filter<-godata_filter[order(godata_filter$pvalue,decreasing = T),]
    godata_filter$Description <- sapply(godata_filter$Description,function(x){
      
      num <- length(stringr::str_split(x," ")[[1]])
      if(num >5){
        x<- paste(paste(stringr::str_split(x," ")[[1]][1:ceiling(num/2)], collapse = " "),
                  paste(stringr::str_split(x," ")[[1]][(ceiling(num/2)+1):num], collapse = " "), sep='\n')
      }
      return(x)
    })
    
  }else if(org=='Rat'){
    library(org.Rn.eg.db)
    
    entrez_gene<-bitr(gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Rn.eg.db)
    enrichGOA <- enrichGO(gene = entrez_gene$ENTREZID, OrgDb = org.Rn.eg.db,keyType ="ENTREZID", ont = on,  pAdjustMethod = "BH",pvalueCutoff  = 1,qvalueCutoff  = 1,readable= TRUE)
    godata<-enrichGOA@result
    godata<-godata[order(godata$pvalue,decreasing = F),]
    
    godata_filter<-godata[1:length,]
    godata_filter<-godata_filter[order(godata_filter$pvalue,decreasing = T),]
    godata_filter$Description <- sapply(godata_filter$Description,function(x){
      
      num <- length(stringr::str_split(x," ")[[1]])
      if(num >5){
        x<- paste(paste(stringr::str_split(x," ")[[1]][1:ceiling(num/2)], collapse = " "),
                  paste(stringr::str_split(x," ")[[1]][(ceiling(num/2)+1):num], collapse = " "), sep='\n')
      }
      return(x)
    })
    
    
  }else if(org=='monkey'){
    library(org.Mmu.eg.db)
    
    entrez_gene<-bitr(gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mmu.eg.db)
    enrichGOA <- enrichGO(gene = entrez_gene$ENTREZID, OrgDb = org.Mmu.eg.db, keyType ="ENTREZID", ont = on,  pAdjustMethod = "BH",pvalueCutoff  = 1,qvalueCutoff  = 1,readable= TRUE)
    godata<-enrichGOA@result
    godata<-godata[order(godata$pvalue,decreasing = F),]
    
    godata_filter<-godata[1:length,]
    godata_filter<-godata_filter[order(godata_filter$pvalue,decreasing = T),]
    godata_filter$Description <- sapply(godata_filter$Description,function(x){
      
      num <- length(stringr::str_split(x," ")[[1]])
      if(num >5){
        x<- paste(paste(stringr::str_split(x," ")[[1]][1:ceiling(num/2)], collapse = " "),
                  paste(stringr::str_split(x," ")[[1]][(ceiling(num/2)+1):num], collapse = " "), sep='\n')
      }
      return(x)
    })
    
    
  }else if(org=='human'){
    library(org.Hs.eg.db)
    
    entrez_gene<-bitr(gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
    enrichGOA <- enrichGO(gene = entrez_gene$ENTREZID, OrgDb = org.Hs.eg.db,keyType ="ENTREZID", ont = on,  pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
    godata<-enrichGOA@result
    godata<-godata[order(godata$pvalue,decreasing = F),]
    
    godata_filter<-godata[1:length,]
    godata_filter<-godata_filter[order(godata_filter$pvalue,decreasing = T),]
    godata_filter$Description <- sapply(godata_filter$Description,function(x){
      
      num <- length(stringr::str_split(x," ")[[1]])
      if(num >5){
        x<- paste(paste(stringr::str_split(x," ")[[1]][1:ceiling(num/2)], collapse = " "),
                  paste(stringr::str_split(x," ")[[1]][(ceiling(num/2)+1):num], collapse = " "), sep='\n')
      }
      return(x)
    })
    
    
  } else{
    
    print('Error org was provided !')
  }
  
  if(graph==TRUE){
    
    print(goplot(enrichGOA))
  }
  
  godata_filter$Description<-factor(godata_filter$Description,levels = godata_filter$Description)
  a<-ggplot(data=godata_filter,aes(y=-log10(pvalue),x=Description))+geom_bar(stat = "identity",fill=col, width = 0.5)+coord_flip()+labs(title=title)+theme(axis.ticks = element_blank())
  print(a)
  result<-list(gene,enrichGOA,a)
  return(result)
  
}

#
go_early <- sc_tl_goterms(gene = degs_ery2$gene[degs_ery2$cluster=='Early_Ery'], org = 'human',
                          title = 'GO terms of Early_Ery', length = 10, col = col_ery[2], on = 'BP', graph = F)

go_early[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))

#
go_erythroblast <- sc_tl_goterms(gene = degs_ery2$gene[degs_ery2$cluster=='Erythroblast'], org = 'human',
                          title = 'GO terms of Erythroblast', length = 10, col = col_ery[1], on = 'BP', graph = F)

go_erythroblast[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))

#
go_late <- sc_tl_goterms(gene = degs_ery2$gene[degs_ery2$cluster=='Late_Ery'], org = 'human', title = 'GO terms of Late_Ery', length = 10, col = col_ery[3], on = 'BP', graph = F)
go_late[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))

#feature genes
library(ggpubr)
ggarrange(plotlist =  purrr::map(c('MYB', 'GATA2', 'GATA1', 'KLF1','NFE2', 
                                   'RPS18', 'RPS9', 'HMGA1', 'PKM', 'MT1H', 'ALAS2', 'UROS','BPGM', 'FECH', 'GYPC'), function(gene){
                                     FeaturePlot(ery, features = gene, order = T)+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6','grey', 'red'))(100))
             }), ncol = 5, nrow = 3)

#cell cycle analysis
genelist <- list(g1sGene = c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"),
                 g2mGene = c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"))

ery <- CellCycleScoring(ery, s.features = toupper(genelist$g1sGene),
                        g2m.features = toupper(genelist$g2mGene))

ery$Phase <- factor(ery$Phase, levels = c('G1', 'S', 'G2M'))
col_cycle <- c('G1' = 'blue', 'S' ='orange3', 'G2M' ='red3')

ery@meta.data %>% dplyr::select(celltype, Phase) %>% table() %>% data.frame() %>% 
  ggplot(aes(celltype, Freq, fill = Phase))+geom_bar(stat = 'identity', position = 'fill')+theme_bw()+
  theme(axis.text = element_text(size = 12, colour = 'black'))+labs(y='Cell Proportion')+scale_fill_manual(values = col_cycle)



####DEGs between pEry1 and  pEry2 in CS6 ####
#Ery
cs6_ery <- cs6_hema2[, cs6_hema2$CellType_V2 %in% c('Primitive_Ery1', 'Primitive_Ery2')]
cs6_ery <- SetIdent(cs6_ery, value = 'CellType_V2')
degs_cs6_ery <- FindAllMarkers(cs6_ery, logfc.threshold = log(1.25), only.pos = T)
degs_cs6_ery <- degs_cs6_ery[degs_cs6_ery$p_val_adj<0.05, ]

cs6_ery_data <- cs6_ery@assays$RNA@data
cs6_ery_data <- cs6_ery_data[rowMeans(cs6_ery_data)>0, ]
cs6_ery_scale <- t(scale(t(cs6_ery_data)))

average_cs6_ery <- sc_tl_average(cs6_ery_data, cs6_ery$CellType_V2)

pheatmap::pheatmap(average_cs6_ery[as.character(na.omit(sc_tl_topgene(degs_cs6_ery, 20))),],
                   cluster_rows = F, cluster_cols = F, color = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))



#### integratative analysis for mega lineage #####
mk <- ys.integrated[, as.character(ys.integrated$Cluster) %in% c('MkP', 'Mk')]
mk <- NormalizeData(mk)
mk <- ScaleData(mk)
mk <- FindVariableFeatures(mk)
hvg_mk <- VariableFeatures(mk)

noise_mk <- sc_tl_detect_noise_genes(mk)
mk <- RunPCA(mk, features = setdiff(VariableFeatures(mk), noise_mk))

mk@assays$RNA@scale.data <- as.matrix(0)
mk <- RunFastMNN(SplitObject(mk, split.by = 'stage'),
                  features = setdiff(hvg_mk, noise_mk))

mk <- RunUMAP(mk, dims = 1:30, reduction = 'mnn')
mk <- FindNeighbors(mk, dims = 1:30, reduction = 'mnn')
mk <- FindClusters(mk, resolution = seq(0.2, 2, 0.2))

#define clusters
mk$celltype <- plyr::mapvalues(mk$RNA_snn_res.2, 0:13,
                               c('MkP', 'MkP', 'thrombopoiesis_Mk', 'thrombopoiesis_Mk', 'thrombopoiesis_Mk', 'thrombopoiesis_Mk', 'MkP', 'MkP', 'thrombopoiesis_Mk',
                                 'thrombopoiesis_Mk', 'thrombopoiesis_Mk', 'thrombopoiesis_Mk', 'Immune_Mk', 'Niche_Mk'))

library(ggsci)
col_mk <- setNames(pal_d3(palette = 'category10')(4), c('MkP', 'thrombopoiesis_Mk', 'Immune_Mk', 'Niche_Mk'))

#feature genes
ggarrange(plotlist =  purrr::map(c('GATA2', 'PLEK', 'NFE2', 'GFI1B', 'PF4', 'ITGA2B', 'ITGB3',
                                   "GP1BA", 'GP5', 'GP9','CD47', "CD9", 'THBS1','MAX', 'MYLK'), function(gene){
                                     FeaturePlot(mk, features = gene, order = T)+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'red'))(100))
                                   }), ncol = 5, nrow = 3)

ggarrange(plotlist =  purrr::map(c('COL1A1', 'COL3A1', 'PTN','COLEC11','VIM',
                                   'PPBP', 'PTCRA','CCL5','HLA-A', 'BMP6'), function(gene){
                                     FeaturePlot(mk, features = gene, order = T)+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'red'))(100))
                                   }), ncol = 5, nrow = 2)

#DEGs
mk <- SetIdent(mk, value = 'celltype')
degs_mk <- FindAllMarkers(mk, logfc.threshold = log(1.25), only.pos = T)
degs_mk <- degs_mk[degs_mk$p_val_adj<0.05, ]

average_mk <- AverageExpression(mk)
degs_mk2 <- sc_tl_filter_deg(degs_mk, average_data = average_mk$RNA)

DotPlot(mk, features = rev(sc_tl_topgene(degs_mk, 10)), group.by = 'celltype')+
  coord_flip()+theme(axis.text.x =element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = BlueAndRed(100))

DotPlot(mk, features = rev(sc_tl_topgene(degs_mk2, 10)), group.by = 'celltype')+
  coord_flip()+theme(axis.text.x =element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = BlueAndRed(100))

#GO analysis
go_thro_mk <- sc_tl_goterms(gene = degs_mk2$gene[degs_mk2$cluster=='thrombopoiesis_Mk'], org = 'human', 
                            title = 'GO terms of thrombopoiesis_Mk', col = col_mk[2], length = 10, on = 'BP', graph = F)

go_thro_mk[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))

#
go_Immune_Mk <- sc_tl_goterms(gene = degs_mk2$gene[degs_mk2$cluster=='Immune_Mk'], org = 'human', 
                            title = 'GO terms of Immune_Mk', col = col_mk[3], length = 10, on = 'BP', graph = F)

go_Immune_Mk[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))

#
go_Niche_Mk <- sc_tl_goterms(gene = degs_mk2$gene[degs_mk2$cluster=='Niche_Mk'], org = 'human', 
                            title = 'GO terms of Niche_Mk', col = col_mk[4], length = 10, on = 'BP', graph = F)

go_Niche_Mk[[3]]+theme_bw()+theme(axis.text = element_text(colour = 'black', size = 12))



#key genes
DotPlot(cs6_hema2[, as.character(cs6_hema2$CellType_V2) %in% c('Primitive_Mk1', 'Primitive_Mk2')], features = c('ITGA2B','THBS1','GP9',
                     'PLEK','PF4','COL1A1', 'COL3A1','VIM','ACTN1','ACTC1','IGF2'), group.by = 'CellType_V2', scale = F, scale.by = 'size')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+scale_color_gradientn(colours = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))


#DEGs between pMeg1 and pMeg2 in CS6
cs6_Mk <- cs6_hema2[, cs6_hema2$CellType_V2 %in% c('Primitive_Mk1', 'Primitive_Mk2')]
cs6_Mk <- SetIdent(cs6_Mk, value = 'CellType_V2')
degs_cs6_Mk <- FindAllMarkers(cs6_Mk, logfc.threshold = log(1.25), only.pos = T)
degs_cs6_Mk <- degs_cs6_Mk[degs_cs6_Mk$p_val_adj<0.05, ]

cs6_Mk_data <- cs6_Mk@assays$RNA@data
cs6_Mk_data <- cs6_Mk_data[rowMeans(cs6_Mk_data)>0, ]
cs6_Mk_scale <- t(scale(t(cs6_Mk_data)))

average_cs6_Mk <- sc_tl_average(cs6_Mk_data, cs6_Mk$CellType_V2)

pheatmap::pheatmap(average_cs6_Mk[as.character(na.omit(sc_tl_topgene(degs_cs6_Mk, 20))),],
                   cluster_rows = F, cluster_cols = F, color = colorRampPalette(c('#E6E6E6', 'orange', 'red', 'brown'))(100))


