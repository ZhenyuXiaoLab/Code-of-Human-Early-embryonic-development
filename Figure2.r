# Fig2: Cellular communication of PrCP

``` R
library(CellChat)
library(patchwork)
library(writexl)
options(stringsAsFactors = FALSE)
library(Seurat)
library(SeuratDisk)
mycols <- readRDS('./data/cols_726.rds')
seurat_obj <- readRDS('./data/seurat_emb_629.rds')
DimPlot(seurat_obj)
seurat_obj$order_617 <- droplevels(seurat_obj$order_617)
cellchat <- createCellChat(object = seurat_obj, meta = seurat_obj@meta.data, group.by = "order_617")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Use the Secreted Signaling set
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- subsetDB(CellChatDB)
# Set the database to be used in the cellchat object
cellchat@DB <- CellChatDB.use

# Preprocessing - filter out signaling-related genes
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the cellchat communication network, note that there are many important parameters to adjust here
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# Extract the cell network framework
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save network results as a table
df.net <- subsetCommunication(cellchat)
save.name <- file.path('./table', "LR_emb.csv")
write.csv(df.net, save.name, row.names = FALSE)
save.name <- file.path('./table', "LR_emb.xlsx")
write_xlsx(df.net, path = save.name)

# Visualize communication between selected categories
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", sources.use = c('PrCP','Connecting.Stalk'))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", sources.use = c('PrCP','Connecting.Stalk'))

pdf("./plot/xinhao_circle.pdf", width=9, height=7.5)
print(netVisual_chord_gene(cellchat, sources.use = c('PrCP','Connecting.Stalk'),  lab.cex = 0.5, legend.pos.y = 30, directional = 1, color.use = mycols[levels(seurat_obj$order_617)], transparency = 0.1))
dev.off()
pdf("./plot/xinhao_stenth_circle.pdf", width=5, height=5)
print(netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", sources.use = c('PrCP','Connecting.Stalk'), color.use = mycols, alpha.edge = 0.9))
dev.off()

dplt <- netVisual_bubble(cellchat, sources.use = c('PrCP','Connecting.Stalk'), remove.isolate = TRUE)
ggsave('./dplot.pdf', dplt, width = 7, height = 6)

df.net <- subsetCommunication(cellchat, sources.use = c('PrCP','Connecting.Stalk'))

# Custom bubble plot
df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")

# Process data, e.g., set p-value categories and transform communication probability
df.net$pval[df.net$pval > 0.05] = 1
df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
df.net$pval[df.net$pval <= 0.01] = 3
df.net$prob[df.net$prob == 0] <- NA
df.net$prob.original <- df.net$prob
df.net$prob <- -1/log(df.net$prob)
df.net <- subset(df.net, source.target%!in%c('PrCP -> PrCP', 'Connecting.Stalk -> Connecting.Stalk'))
df <- df.net

angle.x <- 45
line.on <- TRUE
line.size <- 0.2
color.text.use <- TRUE
color.text <- NULL
angle = c(0, 45, 90)
hjust = c(0, 1, 1)
vjust = c(0, 1, 0.5)
vjust.x = vjust[angle == angle.x]
hjust.x = hjust[angle == angle.x]

# Set color mapping based on user selection
n.colors <- 10
color.heatmap = c("#3288bd", "#66c2a5", "#fcbba1", "#fc9272", "#fb6a4a", '#de2d26', '#a50f15')
direction <- 1
if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
        RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
        (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
} else {
    color.use <- color.heatmap
}
if (direction == -1) {
    color.use <- rev(color.use)
}
font.size <- 10
font.size.title <- 10
show.legend <- TRUE
color.grid <- "grey90"
grid.on <- TRUE
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
title.name <- NULL

# Create bubble plot using ggplot
g <- ggplot(df, aes(x = interaction_name_2, y = source.target, 
        color = prob, size = pval)) + geom_point(pch = 16) + 
        theme_linedraw() + theme(panel.grid.major = element_blank()) + 
        theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
            vjust =  vjust.x), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")

g <- g + scale_radius(range = c(2*min(df$pval), 2*max(df$pval)), 
        breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
            sort(unique(df$pval))], name = "p-value")

# Set bubble color gradient
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),  values = c(0, 0.05, 0.1, 0.2, 0.4,0.5,0.75,1),
        na.value = "white", limits = c(quantile(df$prob, 
            0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
        breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
            1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob."))
} else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
        na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
            title = "Commun. Prob."))
}

# Adjust font size and legend of the chart
g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))

# Add grid lines 
# if (grid.on) {
#     if (length(unique(df$interaction_name_2)) > 1) {
#         g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
#             0.5, 1), lwd = 0.1, colour = color.grid)
#     }
#     if (length(unique(df$source.target)) > 1) {
#         g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$source.target)) - 
#             0.5, 1), lwd = 0.1, colour = color.grid)
#     }
# }

# add the title
if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
}
if (!show.legend) {
    g <- g + theme(legend.position = "none")
}

return(g)

ggsave('./plot/cellchat_prcp_cs.pdf', g, width = 8, height = 5)

# Filter specific signaling pathways
# Extract WNT, SHH, FGF, NODAL pathways
specific_pathways <- c("WNT")
filtered_interactions <- CellChatDB$interaction[grep(paste(specific_pathways, collapse="|"), CellChatDB$interaction$pathway_name, ignore.case=TRUE), ]
filtered_interactions <- CellChatDB$interaction[grep(paste0("\\b(", paste(specific_pathways, collapse="|"), ")\\b"), CellChatDB$interaction$pathway_name, ignore.case=TRUE), ]

library(SCP)
library(ComplexHeatmap)
common_genes <- intersect(c(Nodal_ligand, Nodal_receptor, Nodal_cofactors, Nodal_antagonists), rownames(seurat_obj))
length(common_genes)
common_genes <- intersect(c(HH_ligand, HH_receptor), rownames(seurat_obj))
length(common_genes)
common_genes <- intersect(c(RA_ligand, RA_receptor), rownames(seurat_obj))
length(common_genes)
common_genes <- intersect(unique(c(FGF_ligand, FGF_receptor, FGF_all)), rownames(seurat_obj))
common_genes <- intersect(unique(PCP_all), rownames(seurat_obj))

if (!"AP_distance" %in% colnames(seurat_obj@meta.data)) {
  stop("Seurat object does not contain AP_distance column.")
}

seurat_obj <- RunDynamicFeatures(srt = seurat_obj, lineages = c("AP_distance"), n_candidates = 200, features = common_genes, BPPARAM = BiocParallel::SerialParam())

ht <- DynamicHeatmap(
  srt = seurat_obj, lineages = c("AP_distance"), use_fitted = TRUE, r.sq = 0, padjust = 1, dev.expl = 0, min_expcells = 2,
  heatmap_palette = "viridis", cell_annotation = "celltype",
  pseudotime_label = 25, pseudotime_label_color = "red", pseudotime_label_linetype = 0,
  pseudotime_label_linewidth = 0, height = 9, width = 10,
  cell_annotation_palcolor = list(celltype = c(mycols[levels(seurat_obj$celltype)])),
  show_row_names = TRUE
)

pdf('./plot/heat_PCP.pdf', width = 15, height = 17)
ht$plot
dev.off()

```

