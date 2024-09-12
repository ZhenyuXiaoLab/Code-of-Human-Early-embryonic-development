# FIG3: Differentiation of PrCP

``` R
seurat_g1 <- readRDS('./data/seurat_g1_slingshot_629.rds')
seurat_g1 <- RunSlingshot(srt = seurat_g1, group.by = "order_617", reduction = "umap", start = 'Inter.Epi')
p <- CellDimPlot(seurat_g1, group.by = "order_617", reduction = "umap", lineages_span = 0.5, lineages = paste0("Lineage", 1:3), palcolor = mycols[levels(seurat_my$order_617)], lineages_palcolor = c("#fe9929","#54278f", "#1c9099" )) 
ggsave('./plot/cs6_epi_slingshot_1.pdf', p, width = 7, height = 4.75)
p <- FeatureDimPlot(seurat_g1, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")
ggsave('./plot/cs6_epi_slingshot_2.pdf', p, width = 17, height = 5)

seurat_g1 <- RunDynamicFeatures(srt = seurat_g1, lineages = c("Lineage1", "Lineage2", "Lineage3"), n_candidates = 50, BPPARAM = BPPARAM)

# Function to extract and sort the top 50 genes
get_top_50_genes <- function(lineage_data) {
  # Sort by padjust value
  sorted_genes <- lineage_data[order(lineage_data$padjust), ]
  # Extract the top 50 gene names
  top_50_genes <- rownames(sorted_genes)[1:103]
  return(top_50_genes)
}

seurat_g1 <- RunDynamicFeatures(srt = seurat_g1, lineages = c("Lineage1", "Lineage2", "Lineage3"), features = unique(amnion_markers))
seurat_g1 <- RunDynamicFeatures(srt = seurat_g1, lineages = c("Lineage3"), features = unique(lineage3_genes), BPPARAM = BPPARAM)

# Extract the top 50 dynamic feature genes for each lineage
lineage1_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage1$DynamicFeatures)
lineage2_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage2$DynamicFeatures)
lineage3_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage3$DynamicFeatures)

p <- DynamicPlot(
    srt = seurat_g1, 
    lineages = c("Lineage3"), 
    group.by = "order_617",
    features = emt_genes,
    compare_lineages = TRUE, compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#1c9099",  # Set line color
    point_palcolor = mycols[levels(seurat_g1$order_617)][c('PrCP','Anterior.Epi','Inter.Epi','Posterior.Epi','AM.Ecto','AM')],  # Set point colors
    add_line = TRUE,
    add_interval = TRUE,
    line.size = 1,
    add_point = TRUE,
    add_rug = TRUE,
    flip = FALSE,
    reverse = FALSE,
    x_order = "value",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scp",
    ncol = 8
)
ggsave('./plot/lineage3_genes_emt_test.pdf', p, width = 24, height = 14)
genes <- unique(lineage2_genes, amnion_markers)
seurat_g1_s1 <- RunDynamicFeatures(srt = seurat_g1_s1, lineages = c("Lineage2"), features = unique(genes), BPPARAM = BPPARAM)
p <- DynamicPlot(
    srt = seurat_g1_s1, 
    lineages = c("Lineage2"), 
    group.by = "celltype",
    features = genes,
    compare_lineages = TRUE, compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#54278f", 
    point_palcolor = mycols[levels(seurat_g1$order_617)][c('PrCP','Anterior.Epi','Inter.Epi','Posterior.Epi','AM.Ecto','AM')],  # Set point colors
    add_line = TRUE,
    add_interval = TRUE,
    line.size = 1,
    add_point = TRUE,
    add_rug = TRUE,
    flip = FALSE,
    reverse = FALSE,
    x_order = "value",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scp",
    ncol = 8,
    heatmap_palcolor = heatcolor
)
ggsave('./plot/lineage2_genes_2.pdf', p, width = 24, height = 17)
Idents(seurat_g1) <- seurat_g1$order_617

###calculate RSS score
library(optparse)
op_list <- list(
make_option(c("-i", "--inrds"), type = "character", default = NULL, action = "store", help = "The input of Seurat RDS",metavar="rds"),
make_option(c("-d", "--ident"), type = "character", default = NULL, action = "store", help = "The sample Ident of Seurat object",metavar="idents"),
make_option(c("-s", "--size"),  type = "integer", default = NULL, action = "store", help = "The sample size of Seurat object",metavar="size"),
make_option(c("-l", "--label"), type = "character", default = "out", action = "store", help = "The label of output file",metavar="label"),
make_option(c("-a", "--assay"), type = "character", default = "Spatial", action = "store", help = "The assay of input file",metavar="assay")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

assay <- opt$assay

library(Seurat)
obj <- readRDS(opt$inrds)
if (!is.null(opt$ident)) {
Idents(obj) <-  opt$ident
size=opt$size
if (!is.null(size)) {
obj <- subset(x = obj, downsample = opt$size)
}
saveRDS(obj,"subset.rds")
}
if (is.null(opt$label)) {
label1 <- 'out'
}else{
label1 <- opt$label
}

library(SCopeLoomR)
outloom <- paste0(label1,".loom")
build_loom(file.name = outloom,dgem = obj@assays[[assay]]@counts)
write.table(obj@meta.data,'metadata_subset.xls',sep='\t',quote=F)



library(optparse)
op_list <- list(
make_option(c("-l", "--input_loom"), type = "character", default = NULL, action = "store", help = "The input of aucell loom file",metavar="rds"),
make_option(c("-m", "--input_meta"), type = "character", default = NULL, action = "store", help = "The metadata of Seurat object",metavar="idents"),
make_option(c("-a", "--assay"), type = "character", default = 'Spatial', action = "store", help = "The assay of Seurat object",metavar="assay"),
make_option(c("-c", "--celltype"), type = "character", default = NULL, action = "store", help = "The colname of metadata to calculate RSS",metavar="label")
)
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
celltype <- opt$celltype
message(paste0('细胞类型列为：',celltype))
assay <- opt$assay
loom <- open_loom(opt$input_loom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

meta <- read.table(opt$input_meta,sep='\t',header=T,stringsAsFactor=F)
cellinfo <- meta[,c(opt$celltype,paste0("nFeature_",assay),paste0("nCount_",assay))]
colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"

sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
try({
rssPlot <- plotRSS(rss)
save(regulonAUC,rssPlot,regulons,file='regulon_RSS.Rdata')
})

saveRDS(rss,paste0(celltype,"_rss.rds"))

source('/home/xlyang/python_work/pyscenic_pipeline/00.scripts/function_pyscenic_visualize.R')
plot_pyscenic(inloom='aucell.loom',incolor=incolor,inrss=paste0(celltype,"_rss.rds"),inrds='subset.rds',infun='median', ct.col=celltype,inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(pheatmap)
library(cowplot)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(stringr)



colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
incolor<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))

plot_pyscenic <- function(inloom='aucell.loom',incolor=incolor,inrss='seurat_annotations_rss.rds',inrds='subset.rds',infun='median', ct.col='seurat_annotations',inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50){
###load data
loom <- open_loom(inloom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

embeddings <- get_embeddings(loom)
close_loom(loom)

rss <- readRDS(inrss)
sce <- readRDS(inrds)

##calculate  RSS fc
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
                dat= data.frame(
                regulon  = rownames(rss),
                cluster =  colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, get(infun))
               )
             }))

df$fc = df$sd.1 - df$sd.2

#select top regulon
ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)

ntopgene <- unique(ntopg$regulon)
write.table(ntopgene,'sd_regulon_RSS.list',sep='\t',quote=F,row.names=F,col.names=F)
#plot rss by cluster

#using plotRSS
rssPlot <- plotRSS(rss)
regulonsToPlot <- rssPlot$rowOrder
rp_df <- rssPlot$df

write.table(regulonsToPlot,'rss_regulon.list',sep='\t',quote=F,row.names=F,col.names=F)
write.table(rp_df,'rssPlot_data.xls',sep='\t',quote=F)
nlen <- length(regulonsToPlot)
hei <- ceiling(nlen)*0.4
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
lgroup <- levels(rssPlot$df$cellType)

nlen2 <- length(lgroup)
wei <- nlen2*2
pdf(paste0('regulons_RSS_',ct.col,'_in_dotplot.pdf'),wei,hei)
print(rssPlot$plot)
dev.off()

# sd top gene
anrow = data.frame( group = ntopg$cluster)
lcolor <- incolor[1:length(unique(ntopg$cluster))]
names(lcolor) <- unique(anrow$group)
annotation_colors <- list(group=lcolor)

pn1 = rss[ntopg$regulon,]
pn2 = rss[unique(ntopg$regulon),]
rownames(pn1) <-  make.unique(rownames(pn1))
rownames(anrow) <- rownames(pn1)
scale='row'
hei <- ceiling(length(ntopg$regulon)*0.4)
pdf(paste0('regulon_RSS_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
)
print(
pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons')
)
dev.off()

#plotRSS gene

pn2 = rss[unique(rp_df$Topic),]
scale='row'
hei <- ceiling(length(unique(rp_df$Topic))*0.4)
pdf(paste0('regulon_RSS_in_plotRSS_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn2,scale=scale,show_rownames = T, main='plotRSS unique regulons')
)
dev.off()

#all regulons

hei <- ceiling(length(rownames(rss))*0.2)
pdf(paste0('all_regulons_RSS_in_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(rss,scale=scale,show_rownames = T,main='all regulons RSS')
)
dev.off()
#plot rss by all cells
if (is.null(inregulons)){
inregulons <- regulonsToPlot
}else{
inregulons <- intersect(inregulons,rownames(rss))
regulonsToPlot <- inregulons

}
pn3=as.matrix(regulonAUC@assays@data$AUC)
regulon <- rownames(pn3)
#regulon <- inregulons
pn3 <- pn3[regulon,]
#pn3 <- pn3[,sample(1:dim(pn3)[2],500)]

sce$group1=sce@meta.data[,ct.col]

meta <- sce@meta.data
meta <- meta[order(meta$group1),]
#meta <- meta[colnames(pn3),]
ancol = data.frame(meta[,c('group1')])
colnames(ancol) <- c('group1')
rownames(ancol) <- rownames(meta)
lcolor <- incolor[1:length(unique(ntopg$cluster))]
names(lcolor) <- unique(ntopg$cluster)
annotation_colors <- list(group1 =lcolor)

df1 <- ancol
df1$cell <- rownames(df1)
df1 <- df1[order(df1$group1),]
pn3 <- pn3[,rownames(df1)]
torange=c(-2,2)
pn3 <- scales::rescale(pn3,to=torange)
pn3 <- pn3[,rownames(ancol)]

scale='none'
hei <- ceiling(length(unique(regulon))*0.2)
pdf(paste0('all_regulon_activity_in_allcells.pdf'),10,hei)
print(
pheatmap(pn3,annotation_col = ancol,scale=scale,annotation_colors=annotation_colors,show_rownames = T,show_colnames = F,cluster_cols=F)
)
#pheatmap(pn3,scale=scale,show_rownames = T, show_colnames = F,cluster_cols=F)
dev.off()

#plot in seurat
regulonsToPlot = inregulons
sce$sub_celltype <- sce@meta.data[,ct.col]
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]

cellClusters <- data.frame(row.names = colnames(sce),
                           seurat_clusters = as.character(sce$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(sce),
                        celltype = sce$sub_celltype)

sce@meta.data = cbind(sce@meta.data ,t(sub_regulonAUC@assays@data@listData$AUC[regulonsToPlot,]))
Idents(sce) <- sce$sub_celltype

nlen <- length(regulonsToPlot)
hei <- ceiling(nlen)*0.4
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
nlen2 <- length(unique(sce$sub_celltype))
wei <- nlen2*2
pdf('regulons_activity_in_dotplot.pdf',wei,hei)
print(DotPlot(sce, features = unique(regulonsToPlot)) + coord_flip()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
      scale_color_gradientn(colours = blu)
      )
dev.off()

hei=ceiling(nlen/4)*4
pdf('regulons_activity_in_umap.pdf',16,hei)
print(RidgePlot(sce, features = regulonsToPlot , ncol = 4))
print(VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ))
print(FeaturePlot(sce, features = regulonsToPlot))
dev.off()

grn <- read.table(ingrn,sep='\t',header=T,stringsAsFactors=F)
inregulons1=gsub('[(+)]','',inregulons)

c1 <- which(grn$TF %in% inregulons1)
grn <- grn[c1,]
#edge1 <- data.frame()
#node1 <- data.frame()
pdf(paste0(ntop2,'_regulon_netplot.pdf'),10,10)
for (tf in unique(grn$TF)) {
tmp <- subset(grn,TF==tf)
if (dim(tmp)[1] > ntop2) {
tmp <- tmp[order(tmp$importance,decreasing=T),]
tmp <- tmp[1:ntop2,]
}
node2 <- data.frame(tmp$target)
node2$node.size=1.5
node2$node.colour <- 'black'
colnames(node2) <- c('node','node.size','node.colour')
df1 <- data.frame(node=tf,node.size=2,node.colour='#FFDA00')
node2 <- rbind(df1,node2)


edge2 <- tmp
colnames(edge2) <- c('from','to','edge.width')
edge2$edge.colour <- "#1B9E77"
torange=c(0.1,1)
edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)

graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
  scale_edge_width_continuous(range = c(1,0.2)) +geom_node_point(aes(colour = node.colour, size = node.size))+ theme_void() +
      geom_node_label(aes(label = node,colour = node.colour),size = 3.5, repel = TRUE)
p1 <- p1 + scale_color_manual(values=c('#FFDA00','black'))+scale_edge_color_manual(values=c("#1B9E77"))
print(p1)
}
dev.off()
#plot activity heatmap
meta <- sce@meta.data
celltype <- ct.col
cellsPerGroup <- split(rownames(meta),meta[,celltype])
sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
scale='row'
rss <- regulonActivity_byGroup
hei <- ceiling(length(regulonsToPlot)*0.4)
pn1 <- rss[regulonsToPlot,]
pdf(paste0('regulon_activity_in_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn1,scale=scale,show_rownames = T, main='regulons activity')
)
dev.off()
saveRDS(rss,'all_regulons_activity_in_',ct.col,'.rds')
hei <- ceiling(length(rownames(rss))*0.2)
pdf(paste0('all_regulons_activity_in_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(rss,scale=scale,show_rownames = T,main='all regulons activity')
)
dev.off()

##calculate fc
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
                dat= data.frame(
                regulon  = rownames(rss),
                cluster =  colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, get(infun))
               )
             }))

df$fc = df$sd.1 - df$sd.2

#select top regulon
ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)

ntopgene <- unique(ntopg$regulon)
write.table(ntopgene,'sd_regulon_activity.list',sep='\t',quote=F,row.names=F,col.names=F)

anrow = data.frame( group = ntopg$cluster)
lcolor <- incolor[1:length(unique(ntopg$cluster))]
names(lcolor) <- unique(anrow$group)
annotation_colors <- list(group=lcolor)
pn1 = rss[ntopg$regulon,]
pn2 = rss[unique(ntopg$regulon),]
rownames(pn1) <-  make.unique(rownames(pn1))
rownames(anrow) <- rownames(pn1)
scale='row'
hei <- ceiling(length(ntopg$regulon)*0.4)
pdf(paste0('regulon_activity_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
)
print(
pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons ')
)
dev.off()

}

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(pheatmap)
library(cowplot)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(stringr)



colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
len <- 100
incolor<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))

plot_pyscenic <- function(inloom='aucell.loom',incolor=incolor,inrss='seurat_annotations_rss.rds',inrds='subset.rds',infun='median', ct.col='seurat_annotations',inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50){
###load data
loom <- open_loom(inloom)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)

embeddings <- get_embeddings(loom)
close_loom(loom)

rss <- readRDS(inrss)
sce <- readRDS(inrds)

##calculate  RSS fc
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
                dat= data.frame(
                regulon  = rownames(rss),
                cluster =  colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, get(infun))
               )
             }))

df$fc = df$sd.1 - df$sd.2

#select top regulon
ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)

ntopgene <- unique(ntopg$regulon)
write.table(ntopgene,'sd_regulon_RSS.list',sep='\t',quote=F,row.names=F,col.names=F)
#plot rss by cluster

#using plotRSS
rssPlot <- plotRSS(rss)
regulonsToPlot <- rssPlot$rowOrder
rp_df <- rssPlot$df

write.table(regulonsToPlot,'rss_regulon.list',sep='\t',quote=F,row.names=F,col.names=F)
write.table(rp_df,'rssPlot_data.xls',sep='\t',quote=F)
nlen <- length(regulonsToPlot)
hei <- ceiling(nlen)*0.4
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
lgroup <- levels(rssPlot$df$cellType)

nlen2 <- length(lgroup)
wei <- nlen2*2
pdf(paste0('regulons_RSS_',ct.col,'_in_dotplot.pdf'),wei,hei)
print(rssPlot$plot)
dev.off()

# sd top gene
anrow = data.frame( group = ntopg$cluster)
lcolor <- incolor[1:length(unique(ntopg$cluster))]
names(lcolor) <- unique(anrow$group)
annotation_colors <- list(group=lcolor)

pn1 = rss[ntopg$regulon,]
pn2 = rss[unique(ntopg$regulon),]
rownames(pn1) <-  make.unique(rownames(pn1))
rownames(anrow) <- rownames(pn1)
scale='row'
hei <- ceiling(length(ntopg$regulon)*0.4)
pdf(paste0('regulon_RSS_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
)
print(
pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons')
)
dev.off()

#plotRSS gene

pn2 = rss[unique(rp_df$Topic),]
scale='row'
hei <- ceiling(length(unique(rp_df$Topic))*0.4)
pdf(paste0('regulon_RSS_in_plotRSS_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn2,scale=scale,show_rownames = T, main='plotRSS unique regulons')
)
dev.off()

#all regulons

hei <- ceiling(length(rownames(rss))*0.2)
pdf(paste0('all_regulons_RSS_in_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(rss,scale=scale,show_rownames = T,main='all regulons RSS')
)
dev.off()
#plot rss by all cells
if (is.null(inregulons)){
inregulons <- regulonsToPlot
}else{
inregulons <- intersect(inregulons,rownames(rss))
regulonsToPlot <- inregulons

}
pn3=as.matrix(regulonAUC@assays@data$AUC)
regulon <- rownames(pn3)
#regulon <- inregulons
pn3 <- pn3[regulon,]
#pn3 <- pn3[,sample(1:dim(pn3)[2],500)]

sce$group1=sce@meta.data[,ct.col]

meta <- sce@meta.data
meta <- meta[order(meta$group1),]
#meta <- meta[colnames(pn3),]
ancol = data.frame(meta[,c('group1')])
colnames(ancol) <- c('group1')
rownames(ancol) <- rownames(meta)
lcolor <- incolor[1:length(unique(ntopg$cluster))]
names(lcolor) <- unique(ntopg$cluster)
annotation_colors <- list(group1 =lcolor)

df1 <- ancol
df1$cell <- rownames(df1)
df1 <- df1[order(df1$group1),]
pn3 <- pn3[,rownames(df1)]
torange=c(-2,2)
pn3 <- scales::rescale(pn3,to=torange)
pn3 <- pn3[,rownames(ancol)]

scale='none'
hei <- ceiling(length(unique(regulon))*0.2)
pdf(paste0('all_regulon_activity_in_allcells.pdf'),10,hei)
print(
pheatmap(pn3,annotation_col = ancol,scale=scale,annotation_colors=annotation_colors,show_rownames = T,show_colnames = F,cluster_cols=F)
)
#pheatmap(pn3,scale=scale,show_rownames = T, show_colnames = F,cluster_cols=F)
dev.off()

#plot in seurat
regulonsToPlot = inregulons
sce$sub_celltype <- sce@meta.data[,ct.col]
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]

cellClusters <- data.frame(row.names = colnames(sce),
                           seurat_clusters = as.character(sce$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(sce),
                        celltype = sce$sub_celltype)

sce@meta.data = cbind(sce@meta.data ,t(sub_regulonAUC@assays@data@listData$AUC[regulonsToPlot,]))
Idents(sce) <- sce$sub_celltype

nlen <- length(regulonsToPlot)
hei <- ceiling(nlen)*0.4
blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
nlen2 <- length(unique(sce$sub_celltype))
wei <- nlen2*2
pdf('regulons_activity_in_dotplot.pdf',wei,hei)
print(DotPlot(sce, features = unique(regulonsToPlot)) + coord_flip()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
      scale_color_gradientn(colours = blu)
      )
dev.off()

hei=ceiling(nlen/4)*4
pdf('regulons_activity_in_umap.pdf',16,hei)
print(RidgePlot(sce, features = regulonsToPlot , ncol = 4))
print(VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ))
print(FeaturePlot(sce, features = regulonsToPlot))
dev.off()

grn <- read.table(ingrn,sep='\t',header=T,stringsAsFactors=F)
inregulons1=gsub('[(+)]','',inregulons)

c1 <- which(grn$TF %in% inregulons1)
grn <- grn[c1,]
#edge1 <- data.frame()
#node1 <- data.frame()
pdf(paste0(ntop2,'_regulon_netplot.pdf'),10,10)
for (tf in unique(grn$TF)) {
tmp <- subset(grn,TF==tf)
if (dim(tmp)[1] > ntop2) {
tmp <- tmp[order(tmp$importance,decreasing=T),]
tmp <- tmp[1:ntop2,]
}
node2 <- data.frame(tmp$target)
node2$node.size=1.5
node2$node.colour <- 'black'
colnames(node2) <- c('node','node.size','node.colour')
df1 <- data.frame(node=tf,node.size=2,node.colour='#FFDA00')
node2 <- rbind(df1,node2)


edge2 <- tmp
colnames(edge2) <- c('from','to','edge.width')
edge2$edge.colour <- "#1B9E77"
torange=c(0.1,1)
edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)

graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
  scale_edge_width_continuous(range = c(1,0.2)) +geom_node_point(aes(colour = node.colour, size = node.size))+ theme_void() +
      geom_node_label(aes(label = node,colour = node.colour),size = 3.5, repel = TRUE)
p1 <- p1 + scale_color_manual(values=c('#FFDA00','black'))+scale_edge_color_manual(values=c("#1B9E77"))
print(p1)
}
dev.off()
#plot activity heatmap
meta <- sce@meta.data
celltype <- ct.col
cellsPerGroup <- split(rownames(meta),meta[,celltype])
sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
scale='row'
rss <- regulonActivity_byGroup
hei <- ceiling(length(regulonsToPlot)*0.4)
pn1 <- rss[regulonsToPlot,]
pdf(paste0('regulon_activity_in_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn1,scale=scale,show_rownames = T, main='regulons activity')
)
dev.off()
saveRDS(rss,'all_regulons_activity_in_',ct.col,'.rds')
hei <- ceiling(length(rownames(rss))*0.2)
pdf(paste0('all_regulons_activity_in_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(rss,scale=scale,show_rownames = T,main='all regulons activity')
)
dev.off()

##calculate fc
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
                dat= data.frame(
                regulon  = rownames(rss),
                cluster =  colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, get(infun))
               )
             }))

df$fc = df$sd.1 - df$sd.2

#select top regulon
ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)

ntopgene <- unique(ntopg$regulon)
write.table(ntopgene,'sd_regulon_activity.list',sep='\t',quote=F,row.names=F,col.names=F)

anrow = data.frame( group = ntopg$cluster)
lcolor <- incolor[1:length(unique(ntopg$cluster))]
names(lcolor) <- unique(anrow$group)
annotation_colors <- list(group=lcolor)
pn1 = rss[ntopg$regulon,]
pn2 = rss[unique(ntopg$regulon),]
rownames(pn1) <-  make.unique(rownames(pn1))
rownames(anrow) <- rownames(pn1)
scale='row'
hei <- ceiling(length(ntopg$regulon)*0.4)
pdf(paste0('regulon_activity_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
print(
pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
)
print(
pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons ')
)
dev.off()

}




```

