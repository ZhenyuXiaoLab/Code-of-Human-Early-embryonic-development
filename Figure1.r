# Fig1: Cluster annotation

``` R
library(Seurat)
#library(SeuratData)
#library(writexl) #excel
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)
color_feature <- mako(9999)

####plot umap and heatmap
Dimplot_sort(seurat_obj,group_name = c('celltype'), order_list =  list(celltype = levels(seurat_obj@meta.data$celltype)), color_palette = col_all,save_name = './plot/fig1_umap_711.pdf')
#pheatmap
library(Seurat)
library(pheatmap)
  # 计算平均基因表达
  mean_gene_exp <- as.matrix(
    data.frame(
      Seurat::AverageExpression(seurat_obj,
                                features = c('GSC',"SIX3", "NODAL", "FGF4", "FOXA2", "CHRD", "NOTO", "CER1", "POU5F1", "UCHL1", "SOX2", "OTX2", "TBXT", 'MIXL1', 'MESP1','MSX1',"BMP4", "ID2", "GABRP", "ISL1", "SLN", "MYL7", "DLK1", "LUM", "VIM", "KDR", "GATA6", "SOX17", "FOXA3", "GATA4", "GATA3", "GATA2", "TFAP2A", "PAGE4", "ERVW-1", "CGA", "MFSD2A", "SDC1", "HLA-G", "LAIR2", "APOE", "TTR", "AFP", "TF", "HBZ", "PF4", "PPBP", "GYPC"),
                                group.by = 'celltype',
                                assays = 'Spatial',
                                layer = 'data'
      )
    )
  )
  
  htdf <- scale(t(mean_gene_exp), scale = TRUE, center = TRUE)
Cell_Type <- gsub("^Spatial\\.", "", rownames(htdf))
Cell_Type[8] <- 'Connecting Stalk'
     rownames(htdf) <- Cell_Type
annotation_row_df <- as.data.frame(Cell_Type)
rownames(annotation_row_df) = annotation_row_df$Cell_Type
colors <- mycols[levels(seurat_obj$celltype)]

#names(colors) <- gsub("[- ]", ".", names(mycols))
#colors<-  colors[order1]
colors_list <- list('Cell_Type' = colors)
  p <- pheatmap(htdf, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           #annotation_col = annotation_col,
           annotation_row = annotation_row_df,
           annotation_colors = colors_list,
           #annotation_legend_side = "right",
           show_rownames = TRUE,
           #show_colnames = !reverse,
           color = colorRampPalette(c(blues))(100),
           main = "Marker Genes")
pdf('./plot/fig_marker_heatmap1_t.pdf',width = 12, height = 5)
    print(p)
dev.off()

####get spatial plot of each slice 
library(Seurat)
library(ggplot2)
library(patchwork)

dimplot_spatial_lable <- function(seurat_obj, group_by, split_by, width, height, savename = "merged_spatial_lable.pdf", 
                                  mycols, lable = TRUE, split = FALSE, slice_num = c()) {
    dynamic_cols <- function(number_of_legends) {
        if (number_of_legends <= 10) {
            return(1)
        } else if (number_of_legends <= 20) {
            return(2)
        } else {
            return(3)
        }
    }
    
    split_objects <- SplitObject(seurat_obj, split.by = split_by)
    split_objects <- split_objects[order(as.numeric(names(split_objects)))]
    
    if (length(slice_num) > 0) {
        split_objects <- split_objects[names(split_objects) %in% slice_num]
    }
    
    plots_list <- list()
    
    for (i in names(split_objects)) {
        if (nrow(split_objects[[i]]@meta.data) > 0) {
            message(paste("正在绘制切片:", i))
            p <- DimPlot(split_objects[[i]], reduction = "spatial", 
                         pt.size = 1e-04, raster = FALSE, alpha = 1, group.by = group_by, 
                         cols = mycols, label = lable, label.size = 1.25) + 
                theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                      axis.text.x = element_blank(), axis.text.y = element_blank(), 
                      axis.ticks = element_blank(), axis.line = element_blank()) + 
                ggtitle(i)
            
            if (lable) {
                p <- p + guides(color = guide_legend(ncol = dynamic_cols(length(unique(split_objects[[i]]@meta.data[[group_by]]))), 
                                  override.aes = list(size = 3)))
            } else {
                p <- p + theme(legend.position = "none")
            }
            
            plots_list[[i]] <- p
            
            if (split) {
                if (!dir.exists(savename)) {
                    dir.create(savename, recursive = TRUE)
                }
                ggsave(file.path(savename, paste0(i, ".pdf")), plot = p, 
                       width = width, height = height, device = "pdf", 
                       limitsize = FALSE)
                message(paste("已保存图像:", file.path(savename, paste0(i, ".pdf"))))
            }
        }
    }
    
    if (!split) {
        plots_list <- Filter(function(x) !is.null(x), plots_list)
        n_col <- 7
        n_row <- ceiling(length(plots_list) / n_col)
        reordered_plots <- rev(plots_list)
        reordered_plots <- Filter(function(x) !is.null(x), reordered_plots)
        combined_plot <- wrap_plots(reordered_plots, ncol = n_col)
        ggsave(savename, plot = combined_plot, width = width, 
               height = height, device = "pdf", limitsize = FALSE)
        message(paste("已保存合并图像:", savename))
    }
}
dimplot_spatial_lable(seurat_obj, group_by = 'order_617', split_by = 'slice_num', width = 35, height = 25, savename = './plot/FIG1_labled_spatial.pdf',mycols = mycols,lable = TRUE)
                                  
####GO Enrichment
library(SCP)
p <- FeatureHeatmap(
  srt = seurat_obj,
    group.by = 'celltype',
    group_palcolor = list(mycols[levels(seurat_obj$celltype)]),
  features = markers$gene,
  feature_split =markers$cluster,
    feature_split_palcol = list(mycols[levels(seurat_obj$celltype)]),
  #cell_annotation = c("Phase"),  
    cell_annotation_palcolor = list(c( color_feature[65],color_feature[30], color_feature[45] )),
  #exp_method = "raw",
  #heatmap_palette = "cividis",
    heatmap_palcol = hcols ,
    slot = "counts",
    height = 20, width = 11,border = FALSE
)
                                
library(ggplot2)
library(dplyr)
library(stringr)
plot_combined_go_enrichment <- function(file_paths, colors, use_p_adjust = TRUE) {
  combined_data <- data.frame()
  
  for (file_path in file_paths) {
    df <- read.delim(file_path, header=TRUE)
    cell_type <- str_extract(basename(file_path), "^[^_]+")
    df$cell_type <- cell_type
    combined_data <- rbind(combined_data, df)
  }
  
  if (use_p_adjust) {
    combined_data$neg_log10_p <- -log10(combined_data$p.adjust)
  } else {
    combined_data$neg_log10_p <- -log10(combined_data$pvalue)
  }
  
  color_map <- color_palette
  
  combined_data <- combined_data %>%
    group_by(cell_type) %>%
    arrange(desc(neg_log10_p), .by_group = TRUE) %>%
    ungroup()
  
  p <- ggplot(combined_data, aes(x=reorder_within(Description, neg_log10_p, cell_type),
                            y=neg_log10_p, 
                            fill=cell_type)) +
    geom_bar(stat="identity") +
    geom_text(aes(x=reorder_within(Description, neg_log10_p, cell_type), y=0, label=Description), 
              hjust=0, color="black", size=3.5, position=position_dodge(width=0.9)) +
    scale_fill_manual(values = color_map) +
    coord_flip() +
    scale_x_reordered() +
    theme_minimal() +
    labs(x="GO Term", y="-log10(p-value)") +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 1, strip.position = "left")
  
  return(plot)
}

file_paths <- list.files(path = "./table/sl_go_fig1/", 
                         pattern = "_GO\\.txt$", 
                         full.names = TRUE)

colors <- color_palette[levels(seurat_obj$celltype)]

combined_plot <- plot_combined_go_enrichment(file_paths, colors, use_p_adjust = FALSE)

ggsave('./plot/combined_go_enrichment.pdf', combined_plot, width = 12, height = 20)

library(ggplot2)

plot_bar_sl <- function(path, color, terms = NULL) {
  df <- read.delim(path, header=TRUE)
  
  if (is.null(terms)) {
    message("terms parameter is NULL, automatically selecting the first three GO terms")
    df <- df[order(df$p.adjust), ]
    terms <- df$ID[1:min(3, nrow(df))]
  }
  
  df <- df[df$ID %in% terms, ]
  df$neg_log10_p_adj <- -log10(df$p.adjust)
  
  p <- ggplot(df, aes(x=reorder(Description, neg_log10_p_adj), y=neg_log10_p_adj, label=Description)) +
    geom_bar(stat="identity", fill=color) +
    geom_text(aes(x=Description, y=0, label=Description), hjust=0, color="black", size=3.5, position=position_dodge(width=0.9)) +
    coord_flip() +
    theme_minimal() +
    labs(x=NULL, y=NULL) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

  return(plot)
}

```

