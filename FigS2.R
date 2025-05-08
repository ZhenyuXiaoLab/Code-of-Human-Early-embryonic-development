#code for Fig.S2
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
            message(paste("Plotting slice:", i))
            p <- DimPlot(split_objects[[i]], reduction = "spatial", 
                         pt.size = 1e-04, raster = FALSE, alpha = 1, group.by = group_by, 
                         cols = mycols, label = lable, label.size = 1.25) + 
                theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                      axis.text.x = element_blank(), axis.text.y = element_blank(), 
                      axis.ticks = element_blank(), axis.line = element_blank()) + 
                ggtitle(i) + guides(color = guide_legend(ncol = dynamic_cols(length(unique(split_objects[[i]]@meta.data[[group_by]]))), 
                                    override.aes = list(size = 3)))
            plots_list[[i]] <- p
            
            if (split) {
                # Create directory if it doesn't exist
                if (!dir.exists(savename)) {
                    dir.create(savename, recursive = TRUE)
                }
                # Save each plot to the specified folder
                ggsave(file.path(savename, paste0(i, ".pdf")), plot = p, 
                       width = width, height = height, device = "pdf", 
                       limitsize = FALSE)
                message(paste("Saved plot:", file.path(savename, paste0(i, ".pdf"))))
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
        message(paste("Saved combined plot:", savename))
    }
}

dimplot_spatial_lable(seurat_obj, group_by = 'order_612', split_by = 'slice_num', width = 5, height = 4, savename = './plot/each_slices', mycols = mycols, lable = FALSE, split = TRUE, slice_num = c(1:49))