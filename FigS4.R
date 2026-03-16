#code for Fig.S4
# key Software versions:
# Seurat v5.0.2
# ggplot2 v3.4.4
## code for Fig.S6
# Seurat	5.1.0
# tidyr	1.3.1
# dplyr	1.1.4
# tibble	3.2.1.9016
# ggplot2	3.4.4
	### vlnplot
big_vln_plot <- function (seurat_obj, 
                          genes = c(), 
                          assay = "RNA", 
                          slot = "data",
                          reverse = TRUE, 
                          cols = NULL, 
                          add_line = TRUE, 
                          group.by = NULL,
                          cluster_col = FALSE, 
                          cluster_row = FALSE, 
                          group.level = NULL,
                          hide_triangle = FALSE) {
  
  require(Seurat)
  require(tidyr)
  require(dplyr)
  require(tibble)
  require(ggplot2)
  
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) {
    stop("No valid genes provided or found in the Seurat object.")
  }
  
  if (slot == "scale.data") {
    seurat_obj <- ScaleData(seurat_obj, features = genes, assay = assay)
  }
  
  vars_to_fetch <- c(genes, if (is.null(group.by)) "ident" else group.by)
  dt <- FetchData(seurat_obj, vars = vars_to_fetch, slot = slot, assays = assay)
  
  dt <- dt %>% 
    rownames_to_column("cell") %>% 
    pivot_longer(cols = all_of(genes), names_to = "genes", values_to = "expressions") %>% 
    rename(cluster = if (is.null(group.by)) "ident" else group.by)
  
  dt$genes <- factor(dt$genes, levels = genes)
  
  if (!is.null(group.level)) {
    dt$cluster <- factor(dt$cluster, levels = group.level)
  }
  
  if (is.null(cols)) {
    cols <- (scales::hue_pal())(length(unique(dt$cluster)))
  }
  
  if (cluster_col || cluster_row) {
    require(ggdendro)
    
    if (reverse) {
      if (cluster_col) {
        gene_dist <- dist(t(dt %>% pivot_wider(names_from = genes, values_from = expressions) %>% select(-cell, -cluster)))
        gene_hclust <- hclust(gene_dist)
        gene_order <- gene_hclust$labels[gene_hclust$order]
        dt$genes <- factor(dt$genes, levels = gene_order)
      }
    } else {
      if (cluster_row) {
        cluster_data <- dt %>% 
          group_by(cluster, genes) %>% 
          summarise(mean_expression = mean(expressions, na.rm = TRUE), .groups = "drop") %>% 
          pivot_wider(names_from = genes, values_from = mean_expression) %>% 
          column_to_rownames("cluster")
        
        if (any(is.na(cluster_data))) {
          warning("NAs found in cluster data. Removing NAs for clustering.")
          cluster_data <- na.omit(cluster_data)
        }
        
        if (nrow(cluster_data) > 1) {
          cluster_dist <- dist(cluster_data)
          cluster_hclust <- hclust(cluster_dist)
          cluster_order <- rownames(cluster_data)[cluster_hclust$order]
          dt$cluster <- factor(dt$cluster, levels = cluster_order)
        } else {
          warning("Not enough data for clustering after removing NAs.")
        }
      }
    }
  }
  
  max_expr_clusters <- dt %>% 
    group_by(genes, cluster) %>% 
    summarise(mean_expression = mean(expressions, na.rm = TRUE), .groups = "drop") %>% 
    group_by(genes) %>% 
    summarise(max_cluster = cluster[which.max(mean_expression)], max_expression = max(mean_expression))
  
  p <- ggplot(data = dt, aes(x = if (reverse) expressions else genes, 
                             y = if (reverse) cluster else expressions, 
                             fill = cluster)) +
    geom_violin(scale = "width", draw_quantiles = if (add_line) c(0.25, 0.5, 0.75) else NULL, 
                color = "black", size = 0.45, alpha = 0.8) +
    scale_fill_manual(values = cols, limits = rev(levels(dt$cluster))) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(size = 14, angle = if (reverse) 60 else 45, hjust = if (reverse) 1 else 1), 
          axis.text.y = element_text(size = 16), 
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16), 
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 16, angle = 60), 
          legend.title = element_text(size = 16), 
          legend.text = element_text(size = 15))
  
  custom_breaks <- function(x) { c(0, max(x)) }
  custom_labels <- function(x) { c("0", format(round(max(x), 1), nsmall = 1)) }
  
  if (reverse) {
    p <- p + 
      facet_grid(cols = vars(genes), scales = "free_x") +
      labs(x = "Log Normalized Expression", y = "") +
      scale_x_continuous(breaks = custom_breaks, labels = custom_labels)
  } else {
    p <- p + 
      facet_grid(rows = vars(cluster), scales = "free_y") +
      labs(x = "", y = "Log Normalized Expression") +
      scale_y_continuous(breaks = custom_breaks, labels = custom_labels)
  }
  
  if (!add_line && !hide_triangle) {
    star_data <- dt %>% 
      group_by(genes, cluster) %>% 
      summarise(mean_expression = mean(expressions), .groups = "drop") %>% 
      left_join(max_expr_clusters, by = "genes") %>% 
      filter(cluster == max_cluster)
    
    if (reverse) {
      p <- p + 
        geom_point(data = star_data, aes(x = mean_expression, y = cluster), 
                   shape = 24, size = 2.5, color = "black", fill = "black")
    } else {
      p <- p + 
        geom_point(data = star_data, aes(x = genes, y = mean_expression), 
                   shape = 24, size = 2.5, color = "black", fill = "black")
    }
  }
  
  options(repr.plot.width = 12, repr.plot.height = if (reverse) 4 else 12)
  return(p)
}

p <- big_vln_plot(seurat_obj, 
                  genes = unique(c("GATA4","GATA6","MYL7","LUM", "ISM2","GATA3","GJA5", "KRT8","CGA", 
                                   "KRT18", "INSL4", "PRG2", "MMP12", "PPBP", "DKK1", "CXCL8", "PLAC8", 
                                   "CD82", "KRT7", "ERVFRD-1", "EPCAM", "SOX9","LGR5","PECAM1","ACKR1", 
                                   "RGS5", "KDR", "PTPRC", "PRAP1", "IGF1","HAND2","PAEP", 'ACTA2',
                                   'MYH11','VIM','DCN',"PCNA","MKI67")), 
                  slot = 'scale.data', 
                  reverse = TRUE, 
                  add_line = FALSE, 
                  cols = rev(mycols[unique(seurat_obj$celltype)]), 
                  group.by = 'celltype', 
                  cluster_row = FALSE, 
                  cluster_col = FALSE, 
                  group.level = unique(seurat_obj$celltype), 
                  assay = 'Spatial')
```

### GO
#same as the code provided in EXT.Fig2

### project 
#Before executing the following projection code, we first converted monkey and mouse gene names to their human orthologous gene names.
# Mapping prediction
anchors <- FindTransferAnchors(reference = seurat_obj_ref, query = seurat_obj_my, reference.reduction = "pca",k.filter = NA,mapping.score.k = 50,dims = 1:10,k.anchor = 10,k.score = 10) #,k.filter = NA
seurat_obj_my <- MapQuery(anchorset = anchors, reference = seurat_obj_ref, query = seurat_obj_my,
                          refdata = "celltype", reference.reduction = "pca", 
                          reduction.model = ifelse(use_raw_umap, "umapnew", "umap"))
    
meta_my <- cbind(seurat_obj_my@meta.data, as.data.frame(seurat_obj_my@reductions$ref.umap@cell.embeddings))
meta_ref <- cbind(seurat_obj_ref@meta.data, as.data.frame(seurat_obj_ref@reductions[[ifelse(use_raw_umap, "umapnew", "umap")]]@cell.embeddings))
    
options(repr.plot.width=7, repr.plot.height=6)
p1 <- Dimplot_sort(seurat_obj_ref, group_name = c('celltype'), order_list = 'auto', save_name='', color_palette = mycols, reduction = "umap")
    
library(ggplot2)
p2 <- ggplot(meta_ref, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "grey", size = 0.01, alpha = 0.5) +
  geom_point(data = meta_my, aes(x = refUMAP_1, y = refUMAP_2, color = celltype), size = 0.75, alpha = 1, shape = 17) +
  scale_color_manual(values = mycols) +
  theme_minimal() +
  labs(x = "", y = "", color = "Cell Type") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

