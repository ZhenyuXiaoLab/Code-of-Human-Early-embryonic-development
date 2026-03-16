## code for Fig.S8
#            Seurat       5.0.3
#        SeuratDisk       0.0.9020
#           ggplot2       3.4.4
#            dplyr        1.1.3
#            tibble       3.2.1
#        ggchicklet       0.4.1
#         tidyverse       2.0.0
#            ggsci        3.0.0
#         destiny        3.10.0
#  SingleCellExperiment 1.22.0
#         biomaRt        2.62.0
#         msigdbr        7.5.1
#           fgsea        1.32.2
### GSEA-GO

############ CS6_ave_hypo_GSEA_go
CS6_deg_ave <- read.csv('E:/Tina/BIT/CS6/pseudo/ave_hypo/cs6_ave_deg_norm.csv')
CS6_deg_ave_sub <- CS6_deg_ave[CS6_deg_ave$cluster %in% c('Hypo.1','Hypo.2','ave'),]
head(CS6_deg_ave_sub)

celltype_use <- c('Hypo.1','Hypo.2','ave')
deg <- CS6_deg_ave_sub[which(CS6_deg_ave_sub$cluster %in% celltype_use[1]),]

genelist=deg$avg_log2FC
names(genelist)=toupper(deg$gene)
genelist=sort(genelist,decreasing = T)
head(genelist)


msigdbr_show_species() 
m_df <- msigdbr(species = "Homo sapiens",category = "C5")
head(m_df)

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#summary(fgsea_sets)

fgseaRes<-fgsea(fgsea_sets, stats = genelist, nperm = 1000)

fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% head()

fgseaResTidy <- as.data.frame(fgseaResTidy)
fgseaResTidy$leadingEdge <- as.character(fgseaResTidy$leadingEdge)
write.csv(fgseaResTidy,'E:/Tina/BIT/CS6/pseudo/ave_hypo/cs6_all_deg_gsea_Hypo1.csv')


#######plot
ave_hypo_gsea_go_selected <- read.csv('E:/Tina/BIT/CS6/pseudo/ave_hypo/cs6_all_deg_gsea_select.csv')

  library(ggchicklet) 
  library(tidyverse)
  library(ggsci)
  
  dat <- ave_hypo_gsea_go_selected
  dat$celltype <- factor(dat$celltype, levels = unique(dat$celltype))
  fill_col = c('#7373E1','#CA8EC1','#EB7921')
  wid = 0.8
  y_largest = max(abs(dat$NES))

  p1 = ggplot(data = dat)+
    geom_chicklet(data = dat, color='lightgrey',width = wid,
                  aes(x = pathway_name, y = NES, fill = celltype))+
    scale_fill_manual(values = fill_col)+
    
    theme_classic(base_size = 14)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())+
    geom_text(data = dat,hjust=1,angle=0,
              aes(x = pathway_name, y=-0.05,label = pathway_name))+
    labs(x = NULL,y='NES')+
    theme(aspect.ratio = 1)+
    coord_flip()
```

### DC MAP

###############  Diffusion map ############### 
library(Seurat)
library(destiny)
library(SingleCellExperiment)
library(ggplot2)

cs7_tyser_endoderm <- readRDS('/sdc/hlzhu/cs7_tyser/cs7_tyser_endo.rds')
cs7_tyser_endo_color <- c('DE(NP)_DE2_AVE'='#B67AB5', 'DE(P)_DE1'='#4B102D',
                                           'Hypoblast'='#F8C9B6','YS Endoderm'='#BC806F')

seu_obj <- cs7_tyser_endoderm


emb = Embeddings(seu_obj,"umap")
dmm <- DiffusionMap(emb)

cc <- as.data.frame(dmm$DC1)
cc[,2] <- dmm$DC4
cc$endo_celltype <- seu_obj$endo_celltype
colnames(cc) <- c('DC_1','DC_2','endo_celltype')

diff_map <- ggplot(cc, aes(x=DC_1, y=DC_2, color=endo_celltype))+
  geom_point(size=6)+
  scale_color_manual(values = cs7_tyser_endo_color)+
  theme_classic()
diff_map



### dotplot

rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)

#### Data Loading
cs6_raw <- readRDS("D:/文件/bioinfo/DATA/人CS6data＆mk＆ms/CS6_velocity_2.rds")
# Convert cs6's assay to 'RNA'
cs6 = CreateSeuratObject(counts = cs6_raw@assays$Spatial@counts, meta.data = cs6_raw@meta.data)

load("D:/文件/bioinfo/DATA/人CS6data＆mk＆ms/mouse_ve.rda")
load("D:/文件/bioinfo/DATA/人CS6data＆mk＆ms/monkey_ve_pva.rda")
load("D:/文件/bioinfo/DATA/humanCS7&CS8/Nature_AVE_CS7_Tyser.rda")

cs7_ncb <- readRDS("D:/文件/bioinfo/DATA/humanCS7&CS8/CS7_AVE.rds")
cs8 <- readRDS("D:/文件/bioinfo/DATA/humanCS7&CS8/CS8_AVE.rds")

cs6$source <- "cs6"
cs7$source <- "cs7"
cs8$source <- "cs8"
AVE_Nature$source <- "cs7_nature"
mk_ve$source <- "mk"
ms_ve$source <- "ms"

# Rename required populations and set order
cs6 <- subset(cs6, subset = celltype_new %in% c('Hypo.1','Hypo.2','ave','Anterior pole'))
cs6$ig_celltype <- cs6$celltype_new
cs6@meta.data$ig_celltype[cs6@meta.data$ig_celltype == c('Hypo.1')] <- "Hypo1_cs6b"
cs6@meta.data$ig_celltype[cs6@meta.data$ig_celltype == c('Hypo.2')] <- "Hypo2_cs6b"
cs6@meta.data$ig_celltype[cs6@meta.data$ig_celltype == c('ave')] <- "AVE-like Hypo_cs6b"
cs6@meta.data$ig_celltype[cs6@meta.data$ig_celltype == c('Anterior pole')] <- "Anterior pole_cs6b"
cs6$ig_celltype <- factor(cs6$ig_celltype, levels = c("Hypo1_cs6b","Hypo2_cs6b",'AVE-like Hypo_cs6b','Anterior pole_cs6b'))

cs7_ncb$ig_celltype <- cs7_ncb$celltype_pseu

cs8$ig_celltype <- cs8$meso_celltype
cs8@meta.data$ig_celltype[cs8@meta.data$ig_celltype == c('YS.Endo')] <- "YS.Endo_cs8"
cs8@meta.data$ig_celltype[cs8@meta.data$ig_celltype == c('Visceral.Endo')] <- "Visceral.Endo_cs8"
cs8@meta.data$ig_celltype[cs8@meta.data$ig_celltype == c('AVE')] <- "AVE_cs8"
cs8$ig_celltype <- factor(cs8$ig_celltype, levels = c("YS.Endo_cs8","Visceral.Endo_cs8",'AVE_cs8'))

AVE_Nature$ig_celltype <- AVE_Nature$sub_cluster
AVE_Nature@meta.data$ig_celltype[AVE_Nature@meta.data$ig_celltype == c('YS Endoderm')] <- "YS Endoderm_cs7"
AVE_Nature@meta.data$ig_celltype[AVE_Nature@meta.data$ig_celltype == c('DE(P)')] <- "DE(P)_cs7"
AVE_Nature@meta.data$ig_celltype[AVE_Nature@meta.data$ig_celltype == c('Hypoblast')] <- "Hypoblast_cs7"
AVE_Nature@meta.data$ig_celltype[AVE_Nature@meta.data$ig_celltype == c('DE(NP)')] <- "DE(NP)_cs7"
AVE_Nature$ig_celltype <- factor(AVE_Nature$ig_celltype, levels = c("YS Endoderm_cs7", 'DE(P)_cs7', "Hypoblast_cs7", 'DE(NP)_cs7'))

mk_ve$ig_celltype <- mk_ve$celltype
mk_ve@meta.data$ig_celltype[mk_ve@meta.data$ig_celltype == c('VE')] <- "VE_mk"
mk_ve@meta.data$ig_celltype[mk_ve@meta.data$ig_celltype == c('AVE')] <- "AVE_mk"
mk_ve$ig_celltype <- factor(mk_ve$ig_celltype, levels = c("VE_mk", 'AVE_mk'))

###### Convert mouse gene names
library(biomaRt)

useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", mirror = "www")

# Connect to Ensembl database
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

# Get human-mouse gene comparison table
genes_conversion <- getLDS(attributes = c("mgi_symbol"),  # Mouse genes
                           filters = "mgi_symbol",
                           values = rownames(ms_ve),
                           mart = mouse_mart,
                           attributesL = c("hgnc_symbol"),  # Human genes
                           martL = human_mart)
gene_map <- setNames(genes_conversion$HGNC.symbol, genes_conversion$MGI.symbol)

# Replace gene names in mouse data

expression_matrix <- GetAssayData(ms_ve, slot = "counts")
rownames(expression_matrix)
old_genes <- rownames(expression_matrix)
new_genes <- ifelse(old_genes %in% names(gene_map), gene_map[old_genes], old_genes)
rownames(expression_matrix) <- new_genes

# Check for duplicate genes
duplicated_genes <- rownames(expression_matrix)[duplicated(rownames(expression_matrix))]
if(length(duplicated_genes) > 0) {
  cat("Duplicate gene names: ", duplicated_genes, "\n")
  
  # Merge expression values for duplicate genes, can choose sum, mean, etc.
  expression_matrix <- aggregate(expression_matrix, 
                                 by = list(rownames(expression_matrix)), 
                                 FUN = sum)  # Sum expression values for duplicate genes, or you can choose mean etc.
  
  # Set merged gene names
  rownames(expression_matrix) <- expression_matrix$Group.1
  expression_matrix <- expression_matrix[, -1]  # Delete Group.1 column
}

ms_ve_trans <- CreateSeuratObject(counts = expression_matrix)
# Add original Seurat object metadata to new object
ms_ve$ig_celltype <- ms_ve$celltype
ms_ve@meta.data$ig_celltype[ms_ve@meta.data$ig_celltype == c('EmVE')] <- "EmVE_ms"
ms_ve@meta.data$ig_celltype[ms_ve@meta.data$ig_celltype == c('ExVE')] <- "ExVE_ms"
ms_ve@meta.data$ig_celltype[ms_ve@meta.data$ig_celltype == c('AVE')] <- "AVE_ms"
ms_ve$ig_celltype <- factor(ms_ve$ig_celltype, levels = c("EmVE_ms","ExVE_ms",'AVE_ms'))
metadata <- ms_ve@meta.data
ms_ve_trans@meta.data <- metadata

##### Integrate data
dtlist <- list(cs6,cs7_ncb,cs8,AVE_Nature,mk_ve,ms_ve_trans)
features <- SelectIntegrationFeatures(object.list = dtlist)
# Add genes of interest
my_genes <- c('GATA4','GATA6','AFP','TTR','APOC1','APOA1',
              'SOX17','DKK1','HHEX','APOC3','CER1','LEFTY1',
              'CHRD','NODAL','OTX2','FGF17','GSC')
features <- union(features, my_genes)

dtlist <- lapply(X = dtlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features, npcs = 30)
})
# Perform integration
anchors <- FindIntegrationAnchors(object.list = dtlist, dims = 1:30, reduction = "rpca", k.anchor = 20)
dt_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

############ Dimension reduction and visualization
dt_integrated <- FindVariableFeatures(dt_integrated, selection.method = "vst", nfeatures = 2000)
dt_integrated <- NormalizeData(dt_integrated, normalization.method = 'LogNormalize', scale.factor = 10000)
# Scale data
dt_integrated <- ScaleData(dt_integrated, features = head(VariableFeatures(dt_integrated),1000))
# PCA dimensionality reduction
dt_integrated <- RunPCA(dt_integrated, features = VariableFeatures(dt_integrated))

DimPlot(dt_integrated, reduction = 'pca', group.by = 'source', pt.size = .4, label = T, repel = T) + 
  scale_color_manual(values = vibrant_colors)
DimPlot(dt_integrated, reduction = 'pca', group.by = 'ig_celltype', pt.size = .4, label = T, repel = T, ncol = 3)

dt_integrated <- FindNeighbors(dt_integrated, dims = 1:30, reduction = 'pca')
dt_integrated <- FindClusters(dt_integrated, resolution = 0.5)
dt_integrated <- RunUMAP(dt_integrated, dims = 1:30, reduction = 'pca', resolution = 0.8)

DimPlot(dt_integrated, reduction = 'umap', group.by = 'source', pt.size = .4, label = T, repel = T) + 
  scale_color_manual(values = vibrant_colors)

DimPlot(dt_integrated, reduction = 'umap', group.by = 'ig_celltype', split.by = 'source', pt.size = .4, label = T, repel = T, ncol = 3)

dt_integrated$ig_celltype <- factor(dt_integrated$ig_celltype, levels = c("EmVE_ms","ExVE_ms",'AVE_ms',
                                                                       "VE_mk", 'AVE_mk',
                                                                       "YS.Endo_cs8","Visceral.Endo_cs8",'AVE_cs8',
                                                                       "YS Endoderm_cs7", 'DE(P)_cs7',"Hypoblast_cs7",'DE(NP)_cs7',
                                                                       'Definitive.Endo','Hypoblast','Anterior hypoblast',
                                                                       "Hypo1_cs6b","Hypo2_cs6b",'AVE-like Hypo_cs6b','Anterior pole_cs6b'
))

# Gene lists
list_genes = list('Hypo' = c('PDGFRA','GATA6','GATA4','SOX17'),
                  'ave' = c('LEFTY2','HHEX','CER1', 
                            'PKDCC', 
                            'SHISA2','SFRP1','CDH2','SAMD3', 'OTX2','DKK1',
                            'SAT1','SERPINF1','SPOCK3','CCKBR','BEX1'),
                  'Anterior Pole'=c('POU5F1','CHRD','GSC', 'DDIT4','CRIP1', 'FGF17',
                                    'FOXA2','SHH','EOMES','NODAL','NOTO',  'FGF4', 
                                    'SNRPN','MT1X','SIX3','CYP26A1','LEFTY1'))

# Draw dotplot
p <- DotPlot(dt_integrated, features=rev(list_genes), group.by = 'ig_celltype', scale = TRUE) +
  theme(panel.grid.major = element_line(color = "#E4EEED", size = 0.5),  # Major grid lines
        panel.grid.minor = element_line(color = "#E4EEED", size = 0.2),  # Minor grid lines
        panel.background = element_rect(fill = "white", color = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.border = element_rect(color = "black"),
         panel.spacing = unit(1, "mm"),
         axis.title = element_blank(),
  )+
  scale_color_gradientn(values = seq(0, 1, length.out = 7), 
                        # colours = c('#f6f6f6','#f5f2f3','#f0ebeb','#F1E6E6','#f8e7eb','#dc8a9a','#be3b56','#d15c6c','#8d192b')) +
                        colours = c('#ACC7DC','#FAF8F4','#ECCDD2', '#DEA2AF')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_size_continuous(range = c(1,4)) 

p

# Split by developmental stage
p$data <- p$data %>% filter(id %in% c("Hypo1_cs6b", "Hypo2_cs6b", "AVE-like Hypo_cs6b", 'Anterior pole_cs6b'))
ggsave("2ave_cs6_dot.pdf", plot = p, width = 10, height = 2.2)

p$data <- p$data %>% filter(id %in% c("Definitive.Endo", "Hypoblast", "Anterior hypoblast"))
ggsave("2ave_cs7_ncb_dot.pdf", plot = p, width = 10, height = 2.1)

p$data <- p$data %>% filter(id %in% c("YS Endoderm_cs7", "DE(P)_cs7", "Hypoblast_cs7", "DE(NP)_cs7"))
ggsave("2ave_cs7_dot.pdf", plot = p, width = 10, height = 2.2)

p$data <- p$data %>% filter(id %in% c("YS.Endo_cs8", "Visceral.Endo_cs8", "AVE_cs8"))
ggsave("2ave_cs8_dot.pdf", plot = p, width = 10, height = 2.1)

p$data <- p$data %>% filter(id %in% c("VE_mk", "AVE_mk"))
ggsave("2ave_mk_dot.pdf", plot = p, width = 10, height = 2)

p$data <- p$data %>% filter(id %in% c("AVE_ms", "ExVE_ms", 'EmVE_ms'))
ggsave("2ave_ms_dot.pdf", plot = p, width = 10, height = 2.1)