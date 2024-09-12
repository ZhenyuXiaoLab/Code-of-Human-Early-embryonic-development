########integrated data from cs6 and cs7 ############
com_cs6_cs7 <- merge(cs7_exmc[co_genes_cs6_cs7, ],
                     cs6_exmc[co_genes_cs6_cs7,
                              sample(colnames(cs6_exmc), ncol(cs7_exmc))]) #subsample cs6 data

com_cs6_cs7v3 <- SCTransform(com_cs6_cs7v3, vars.to.regress = c('batch', 'nCount_RNA'))
com_cs6_cs7v3 <- RunFastMNN(SplitObject(com_cs6_cs7v3, split.by = 'stage')[c('CS7', 'CS6')], k=50)

com_cs6_cs7v3 <- RunUMAP(com_cs6_cs7v3, reduction = 'mnn', dims = 1:30, local.connectivity = 20, min.dist = 0.5, spread = 1)
DimPlot(com_cs6_cs7v3, group.by = 'cluster_final', cols = all_colors$col_cs6_cs7)+
DimPlot(com_cs6_cs7v3, group.by = 'stage', cols = all_colors$group)

##############monocle3 analysis
library(monocle3)
com_cs6_cs7v3@meta.data <- data.frame(com_cs6_cs7v3@meta.data, t(com_cs6_cs7v3@assays$RNA@data), check.names = F)
com_cs6_cs7_monocle3 <- as.cell_data_set(com_cs6_cs7v3[, com_cs6_cs7v3$cluster_final %in% c("pMP_CS6",
                                                                                       "pMP_CS7",
                                                                                       "YSMP_CS7",
                                                                                       "Mac_CS7",
                                                                                       "EC_CS7")])


com_cs6_cs7_monocle3 <- cluster_cells(com_cs6_cs7_monocle3, reduction_method = 'UMAP',
                                      cluster_method = 'louvain', resolution = 1)

com_cs6_cs7_monocle3 <- learn_graph(com_cs6_cs7_monocle3, use_partition = FALSE)

com_cs6_cs7_monocle3 <- order_cells(com_cs6_cs7_monocle3, reduction_method = 'UMAP',
                                    root_cells = colnames(com_cs6_cs7v3)[com_cs6_cs7v3$cluster_final %in% c("EC_CS7","pMP_CS6")])

plot_cells(com_cs6_cs7_monocle3, color_cells_by = 'cluster_final',
           show_trajectory_graph = T, label_branch_points = F,
           label_roots = F, label_leaves = F, 
           label_cell_groups = F, label_groups_by_cluster = F,
           labels_per_group = F, rasterize = T, cell_size = 1.2)+
  scale_color_manual('Cluster',values = all_colors$col_cs6_cs7)

plot_cells(com_cs6_cs7_monocle3, color_cells_by = 'cluster_final',
           show_trajectory_graph = T, label_branch_points = F,
           label_roots = F, label_leaves = F, 
           label_cell_groups = F, label_groups_by_cluster = F,
           labels_per_group = F, rasterize = T, cell_size = 1.2)+
  scale_color_manual('Cluster',values = all_colors$col_cs6_cs7)+


plot_cells(com_cs6_cs7_monocle3, color_cells_by = 'stage',
             show_trajectory_graph = T, label_branch_points = F,
             label_roots = F, label_leaves = F, 
             label_cell_groups = F, label_groups_by_cluster = F,
             labels_per_group = F, rasterize = T, cell_size = 1.2)+
  scale_color_manual('group',values = as.character(col_group))+
  
plot_cells(com_cs6_cs7_monocle3, color_cells_by = 'pseudotime',
             show_trajectory_graph = T, label_branch_points = F,
             label_roots = F, label_leaves = F, 
             label_cell_groups = F, label_groups_by_cluster = F,
             labels_per_group = F, rasterize = T, cell_size = 1.2)+
  scale_color_gradientn("Pseudotime",colours = colors.use$gradient)

###############
FeaturePlot(com_cs6_cs7v3[, com_cs6_cs7v3$cluster_final %in% c("pMP_CS6",
                                                               "pMP_CS7",
                                                               "YSMP_CS7",
                                                               "Mac_CS7",
                                                               "EC_CS7")], 
            features = c('KDR',  'S100P', 'MPO', 'AZU1', 'PTPRC', 'CSF1R', 'C1QA', 'CD163'),
            cols = c('grey', 'red'), ncol = 4, order = T)

###############
w1 <- com_cs6_cs7v3[, as.character(com_cs6_cs7v3$cluster_final) %in% c('pMP_CS6',
                                                                       'pMP_CS7',
                                                                    'Mac_CS7')]

w1$pseudotime <- com_cs6_cs7_monocle3@principal_graph_aux$UMAP$pseudotime[colnames(w1)]
w1$cluster_monocle <- com_cs6_cs7_monocle3@clusters$UMAP$clusters[colnames(w1)]


####
dynamic_genes <- sc_tl_find_gene_along_paseudotime(pseudotime_data = w1@meta.data[, c('pseudotime', 'cluster_final')],
                                                   data = w1@assays$RNA@data,
                                                   col = list(cluster_final = all_colors$col_cs6_cs7, pseudotime = colors.use$gradient),
                                                   k = 5, pval = 1e-2)

w1_data <- w1@assays$RNA@data
phmat <- t(scale(t(w1_data)))
phmat <- phmat[, colnames(w1)[order(w1$pseudotime)]]
w1_smooth_data <- sc_tl_smooth_data(phmat, k=10)

dynamic_genes <- dynamic_genes[[2]]
w1_anno <- w1@meta.data

ph_w1 <- pheatmap(w1_smooth_data[intersect(rownames(dynamic_genes)[order(dynamic_genes$gam.pval, decreasing = F)][1:1000], rownames(phmat)),
                        rownames(w1_anno)[order(w1_anno$pseudotime, decreasing = F)]],
                  cluster_rows = T, cluster_cols = F,
                  cutree_rows = 4, clustering_method = 'ward.D2',
                  color = colors.use$BluewhiteRed, annotation_col = w1_anno[, c('cluster_final', 'pseudotime')],
                  annotation_colors = list(cluster_final = all_colors$col_cs6_cs7, pseudotime = colors.use$gradient),
                  show_rownames = F, show_colnames = F)

w1_gene_anno <- data.frame(cutree(ph_w1$tree_row, k = 4))
names(w1_gene_anno) <- 'Gene_pattern'
table(w1_gene_anno$Gene_pattern)

# w1_gene_anno$Gene_pattern <- plyr::mapvalues(w1_gene_anno$Gene_pattern, c(2, 1, 3, 4),
#                                              paste('Pattern', c(1, 2, 3, 4), sep=""))

w1_gene_anno$Gene_pattern <- plyr::mapvalues(w1_gene_anno$Gene_pattern, c(2, 3, 1, 4),
                                             paste('Pattern', c(1, 2, 3, 4), sep=""))

all_colors$gene_pattern <- c('Pattern1' = 'blue3',
                             'Pattern2' = 'darkgreen',
                             'Pattern3' = 'orange3',
                             'Pattern4' = 'brown')

pheatmap(w1_smooth_data[rownames(w1_gene_anno)[order(w1_gene_anno$Gene_pattern, decreasing = F)],
                        rownames(w1_anno)[order(w1_anno$pseudotime, decreasing = F)]],
         cluster_rows = F, cluster_cols = F,
         cutree_rows = 4, clustering_method = 'ward.D2',
         color = colors.use$BluewhiteRed, annotation_col = w1_anno[, c('cluster_final', 'pseudotime')],
         annotation_colors = list(cluster_final = all_colors$col_cs6_cs7[c(3, 9, 11)],
                                  pseudotime = colors.use$gradient,
                                  Gene_pattern = all_colors$gene_pattern),
         show_rownames = F, show_colnames = F, annotation_row = w1_gene_anno, annotation_names_row = F)


