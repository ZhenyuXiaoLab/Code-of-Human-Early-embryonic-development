#code for Fig.S3
# Software versions:
# clusterProfiler v4.1.2
# ggplot2 v3.4.4
# igraph v2.0.3
###The function for plotting feature plots utilizes optimized functions in scCustomize. 
#For details, please refer to the tutorial provided by the author.

### GO
``` R
library(dplyr)
library(readr)
library(fs)

extract_DEGs <- function(dataFile, resultdir = '') {
  # Read data file
  data <- read_csv(dataFile)
  
  # Check if output directory exists, create if not
  if (!dir_exists(resultdir)) {
    dir_create(resultdir)
  }
  
  # Check if 'cluster' column exists in data
  if ("cluster" %in% names(data)) {
    # Group data
    data_list <- data %>%
      group_by(cluster) %>%
      group_split()
    for (data_subset in data_list) {
      cluster_name <- unique(data_subset$cluster)
      # Clean cluster names
      cluster_cleaned <- gsub('"', '', cluster_name)
      cluster_cleaned <- gsub('/', '_', cluster_cleaned)
      cluster_cleaned <- gsub(' ', '_', cluster_cleaned)
      # Create filename
      file_name <- paste0(resultdir, "/", cluster_cleaned, "sig_genes.txt")
      # Extract data and save to file
      data_subset %>%
        dplyr::select(gene, avg_log2FC) %>%
        write_tsv(file_name, col_names = TRUE)
    }
  } else {
    # Add 'change' column based on avg_log2FC
    data <- data %>%
      mutate(change = ifelse(avg_log2FC > 0, 'up', 'down'))
    
    # Filter for 'up' and 'down' separately and save to respective files
    data %>% 
      filter(change == 'up') %>%
      select(gene, avg_log2FC) %>%
      write_tsv(paste0(resultdir, "/up_sig_genes.txt"), col_names = TRUE)
    
    data %>% 
      filter(change == 'down') %>%
      select(gene, avg_log2FC) %>%
      write_tsv(paste0(resultdir, "/down_sig_genes.txt"), col_names = TRUE)
  }
}
#extract_DEGs('./cs8_markers.csv','./')
library('Seurat')
## Execute GO_sub.sh before running this
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(writexl)
go_gene_search <- function(path = './', result.dir = 'GO_table') {
  # Ensure the result directory exists
  if (!dir.exists(result.dir)) {
    dir.create(result.dir)
  }

  # Find files ending with sig_genes.txt in specified path
  files <- list.files(path, pattern = "sig_genes.txt$", full.names = TRUE)

  # Process each file
  for (file in files) {
    # Read file contents
    rt <- read.table(file, sep="\t", check.names=FALSE, header=TRUE)
    genes <- as.vector(rt[,1])

    # Remove NA values before Entrez ID lookup
    genes_no_na <- as.vector(genes[!is.na(genes) & genes != ""])
    print(file)
    # Find corresponding Entrez IDs
    entrezIDs <- mget(genes_no_na, org.Hs.egSYMBOL2EG, ifnotfound=NA)
    entrezIDs <- as.character(entrezIDs)

    # Combine original data with Entrez IDs
    out <- cbind(rt, entrezID=entrezIDs)

    # Construct output filename
    filename <- basename(file)
    file_label <- gsub("sig_genes.txt", "", filename)
    output_file <- paste0(file_label, "go.txt")
    output_filepath <- file.path(result.dir, output_file)

    # Write combined output to disk
    write.table(out, file=output_filepath, sep="\t", quote=FALSE, row.names=FALSE)
  }

  cat("Transform completed, result saved in:", result.dir, "\n")
}

# Example
go_gene_search(path = './GO_table/', result.dir= './GO_table/GO_transfer')
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(writexl)
go_enrich <- function(path = './', input.pattern = 'go.txt', result.dir = 'GO_enrichment_results') {
  # Ensure result directory exists
  if (!dir.exists(result.dir)) {
    dir.create(result.dir)
  }

  # Find files matching pattern in specified path
  input_files <- list.files(path, pattern = input.pattern, full.names = TRUE)

  # Process each input file
  for(input_file in input_files) {
    # Read input file contents
    rt <- read.table(input_file, sep="\t", header=TRUE, check.names=FALSE)

    # Remove rows with NA entrezID
    rt <- rt[!is.na(rt$entrezID),]
    genes <- rt$entrezID

    # Perform GO enrichment analysis
    kk <- enrichGO(gene          = genes,
                   OrgDb         = org.Hs.eg.db,
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   ont           = "all",
                   readable      = TRUE)

    # Save enrichment results
    output_file <- gsub(pattern = input.pattern, replacement = "_GO_enrichment.txt", x = basename(input_file))
    ii <- as.data.frame(kk)
    write.table(kk, file = file.path(result.dir, output_file), sep = "\t", quote = FALSE, row.names = FALSE)
    write_xlsx(ii, paste0(result.dir,'/', output_file, ".xlsx"))

    # Create barplot
    barplot_file <- gsub(pattern = ".txt", replacement = "_barplot.pdf", x = output_file)
    print(barplot_file)
    bp <- barplot(kk, drop = TRUE, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free')
    show(bp)
    ggsave(bp,filename = file.path(result.dir, barplot_file), width = 10,height = 12)

    # Create bubble chart
    bubble_file <- gsub(pattern = ".txt", replacement = "_bubble.pdf", x = output_file)
    dp <- dotplot(kk, showCategory = 10, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free')
    show(dp)
    ggsave(dp,filename = file.path(result.dir, bubble_file), width = 10,height = 12)
  }

  cat("GO enrichment analysis completed, results saved in", result.dir, "\n")
}

# Example usage - note input.pattern should match files from previous function
# go_enrich(path = './GO_1/GO_result1/', input.pattern = 'go.txt', result.dir= './GO_1/GO_enrichment_results_2')

data <- read.table("./GO_1/pmega_GO_enrichment.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
net_plot <- function(data, level1, level2, seurat_obj, ident = 'Erythroblast', save.path = '', loop.time = 10) {
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(fields) # For drawing color bars
  #library(extrafont) # Optional font configuration
   # font_import(prompt = FALSE)
   #fonts()
  #loadfonts(device = "pdf")
  # Ensure save directory exists
  if (save.path != '' && !dir.exists(save.path)) {
    dir.create(save.path, recursive = TRUE)
  }

  for (seed_val in 1:loop.time) {
    set.seed(seed_val)

    # Process data by given column names
    tidy_data <- data %>%
      dplyr::select(!!sym(level1), !!sym(level2)) %>%
      tidyr::separate_rows(!!sym(level2), sep = "/") %>%
      dplyr::rename(GO = !!sym(level1), Gene = !!sym(level2))

    # Create graph object
    g <- graph_from_data_frame(tidy_data, directed = TRUE)

    # Customize vertex appearance
    V(g)$size <- ifelse(V(g)$name %in% tidy_data$GO, 20, 10)
    V(g)$color <- "black" # Default color before expression level coloring
    V(g)$frame.color <- NA
    V(g)$label <- V(g)$name
    V(g)$label.color <- "black"
    V(g)$label.cex <- ifelse(V(g)$name %in% tidy_data$GO, 4.5, 3.8) # Adjust font size
    #V(g)$label.cex <- ifelse(V(g)$name %in% tidy_data$GO, 1, 0.8) # Adjust font size
    V(g)$label.font <- ifelse(V(g)$name %in% tidy_data$GO, 2, 2) # Bold for GO terms

    # Reverse all edge directions (optional)
    g <- reverse_edges(g)

    # Customize edge appearance
    E(g)$curved <- 0.3 # Set curvature
    E(g)$arrow.size <- 0.8 # Set arrow size
    E(g)$color <- "darkgray"
    E(g)$width <- 2.8
    #par(family="arial") # Specific font

    # Define layout, multiple options available
    cat("Calculating layout, seed =", seed_val, "\n")
    #layout <- layout_with_drl(g, options = list(maxiter = 2000))
    layout  <- layout_nicely(g, dim = 2) # My preferred option
    #layout <- layout_with_fr(g) 
    # layout <- layout_as_tree(g) 
    # layout <- layout_with_fr(g)  * 2  # Scale layout to make edges thinner, clearer for many nodes

    # Find markers and merge expression data # This code may be redundant, consider moving outside loop
    cat("Calculating DEGs for color mapping, ident =", ident, "\n")
    biaoda <- FindMarkers(seurat_obj, ident.1 = ident, features = unique(tidy_data$Gene))
    biaoda$Gene <- row.names(biaoda)
    merged_data <- tidy_data %>%
      left_join(biaoda, by = "Gene")

    # Create color mapping function for expression data
    min_log2fc <- min(merged_data$avg_log2FC, na.rm = TRUE)
    max_log2fc <- max(merged_data$avg_log2FC, na.rm = TRUE)
     #my_palette <- colorRampPalette(c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8","#253494"))(30) # Yellow-blue
      my_palette <- colorRampPalette(c("#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf", "#1c9099","#016c59"))(30)
     # my_palette <- colorRampPalette(c("#ffffd4", "#fee391", "#fec44f", "#fe9929", "#d95f0e","#993404"))(30) # Brown
    get_color <- function(log2fc, min_val, max_val) {
      sapply(log2fc, function(val) {
        if (is.na(val)) return("grey")
        index <- findInterval(val, seq(min_val, max_val, length.out = 30))
        my_palette[index]
      })
    }

    # Apply color mapping to vertices
    V(g)$color <- sapply(V(g)$name, function(name) {
      gene_log2fc <- unique(merged_data$avg_log2FC[merged_data$Gene == name])
      if (length(gene_log2fc) == 1 && !is.na(gene_log2fc)) {
        get_color(as.numeric(gene_log2fc), min_log2fc, max_log2fc)
      } else {
         "#7bccc4"
        #"#FFAD5C" # Default color when mapping fails (genes should all map successfully)
      }
    })

    # Set output PDF filename
    file_name <- if (save.path == '') {
      paste0("net_plot_seed_", seed_val, ".pdf")
    } else {
      paste0(save.path, "/net_plot_seed_", seed_val, ".pdf")
    }

    # Output to PDF file
    pdf(file_name, width = 30, height = 30) #,family = 'ArialMT'

    # Draw graph
    plot(g, layout = layout, main = paste("GO_net_plot_seed", seed_val),
         edge.width = E(g)$width,
         edge.color = E(g)$color,
         vertex.frame.color = V(g)$frame.color,
         vertex.label = V(g)$label, vertex.label.dist = 0,
         vertex.size = V(g)$size) #, vertex.label.cex = V(g)$label.cex,vertex.label.family="ArialMT"

# Draw color bar
# Save old graphics parameters
old_par <- par(no.readonly = TRUE)
# Adjust margins to leave space at bottom (numbers can be adjusted as needed)
par(mar = c(6, 4, 4, 2) + 0.1)
# Draw horizontal color bar at bottom of current plot
image.plot(z = matrix(0:2, nrow = 1), col = my_palette, legend.only = TRUE,
           horizontal = TRUE, axis.args = list(at = seq(0, 2, by = 0.5), labels = seq(0, 2, by = 0.5), las = 1))
# Add color bar title
mtext("Log2FC", side = 1, line = 5)
# Restore old graphics parameters
par(old_par)
    
    dev.off()
    # Draw graph
    plot(g, layout = layout, main = paste("GO_net_plot_seed", seed_val),
         edge.width = E(g)$width,
         edge.color = E(g)$color,
         vertex.frame.color = V(g)$frame.color,
         vertex.label = NA, vertex.label.dist = 0,
         vertex.size = V(g)$size, vertex.label.cex = V(g)$label.cex)

# Draw color bar
# Save old graphics parameters
old_par <- par(no.readonly = TRUE)
# Adjust margins to leave space at bottom (numbers can be adjusted as needed)
par(mar = c(6, 4, 4, 2) + 0.1)
# Draw horizontal color bar at bottom of current plot
image.plot(z = matrix(0:2, nrow = 1), col = my_palette, legend.only = TRUE,
           horizontal = TRUE, axis.args = list(at = seq(0, 2, by = 0.5), labels = seq(0, 2, by = 0.5), las = 1))
# Add color bar title
mtext("Log2FC", side = 1, line = 5)
# Restore old graphics parameters
par(old_par)
    
  }
}

# Example function call
# net_plot(data = data, level1 = 'Description', level2 = 'geneID', seurat_obj = seurat_obj_GO, ident = 'MP',save.path = 'testnet')
```



``` R
library(SCP)
p <- FeatureHeatmap(
  srt = seurat_obj,
    group.by = 'celltype',
    group_palcolor = list(mycols[levels(seurat_obj$celltype)]),
  features = markers$gene,
  feature_split =markers$cluster,
    feature_split_palcol = list(mycols[levels(seurat_obj$celltype)]),
  cell_annotation = c("Phase"),  # The cell annotation you want to add
    cell_annotation_palcolor = list(c( color_feature[65],color_feature[30], color_feature[45] )),
  #exp_method = "raw",
  #heatmap_palette = "cividis",
    heatmap_palcol = hcols ,
    slot = "scale.data",
    height = 20, width = 11,border = FALSE,
      ht_params = list(
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
       gap = unit(0, "mm"),
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm")
  )
)
ggsave('./plot/feature_plot_711.pdf',p$plot,width = 16,height = 25)

######Bar plot
## GO bar plot visualization
mycols <- readRDS('./data/cols_726.rds')
library(ggplot2)
library(dplyr)
library(stringr)

plot_combined_go_enrichment <- function(file_paths, colors, use_p_adjust = TRUE) {
  # Initialize an empty dataframe to store all data
  all_data <- data.frame()
  
  for (file_path in file_paths) {
    df <- read.delim(file_path, header=TRUE)
    # Extract filename as cell type
    cell_type <- str_extract(basename(file_path), "^[^_]+")
    df$cell_type <- cell_type
    all_data <- rbind(all_data, df)
  }
  
  # Calculate -log10(p.adjust) or -log10(pvalue)
  if (use_p_adjust) {
    all_data$neg_log10_p <- -log10(all_data$p.adjust)
  } else {
    all_data$neg_log10_p <- -log10(all_data$pvalue)
  }
  
  # Assign colors to each cell_type
  color_map <- colors
  
  # Sort dataframe
  all_data <- all_data %>%
    group_by(cell_type) %>%
    arrange(desc(neg_log10_p), .by_group = TRUE) %>%
    ungroup()
  
  # Generate combined bar plot
  p <- ggplot(all_data, aes(x=reorder_within(Description, neg_log10_p, cell_type), 
                            y=neg_log10_p, 
                            fill=cell_type)) +
    geom_bar(stat="identity") +
    geom_text(aes(x=reorder_within(Description, neg_log10_p, cell_type), y=0, label=Description), 
              hjust=0, color="black", size=3.5, position=position_dodge(width=0.9)) +
    scale_fill_manual(values = color_map) +
    coord_flip() +
    scale_x_reordered() + # Reorder x-axis
    theme_minimal() +
    labs(x="GO Term", y="-log10(p-value)") + # Display x-axis title
    theme(
      legend.position = "none", # Remove legend
      panel.grid.major = element_blank(), # Remove grid lines
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(), # Remove y-axis text
      axis.ticks.y = element_blank() # Remove y-axis ticks
    ) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 1, strip.position = "left")
  
  return(p)
}

# Get all GO enrichment result files in directory
file_paths <- list.files(path = "./table/sl_go_fig1/", 
                         pattern = "_GO\\.txt$", 
                         full.names = TRUE)

# Predefined set of colors (assuming you already have this color vector)
colors <- mycols[levels(seurat_obj$celltype)]

# Generate combined bar plot
combined_plot <- plot_combined_go_enrichment(file_paths, colors, use_p_adjust = FALSE)

# Save plot
ggsave('./plot/combined_go_enrichment.pdf', combined_plot, width = 12, height = 20)
```



### cell cycle

``` R
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

#' Analyze and visualize cell cycle phases using Seurat metadata
#'
#' @param meta_data DataFrame containing cell metadata
#' @param cluster_col Column name containing cluster/cell type information
#' @param color_palette Color palette for cell types
#' @param threshold Score threshold for phase assignment (default=3)
#' @return DataFrame with cell cycle phase assignments
#' @examples
#' cycle_df <- cell_cycle_analyze_seurat(seurat_obj@meta.data, "celltype", mycols, threshold=0.05)
cell_cycle_analyze_seurat <- function(meta_data, cluster_col, color_palette, threshold = 3) {
  
  # Extract S.Score and G2M.Score from metadata
  g1sScore <- meta_data$S.Score
  g2mScore <- meta_data$G2M.Score
  cell_type <- meta_data[[cluster_col]]
  
  # Create dataframe for visualization
  cell_cycle_data <- data.frame(G1S = g1sScore, G2M = g2mScore, cell_type = cell_type)
  
  # Assign cell cycle phases based on scores:
  # Q: Both scores below threshold (Quiescent)
  # G1: Only G1/S score above threshold
  # S: Both scores above threshold with G1/S > G2/M
  # G2/M: Both scores above threshold with G2/M > G1/S
  cell_cycle_data$phase <- "G1"
  cell_cycle_data$phase[g1sScore < threshold & g2mScore < threshold] <- "Q"
  cell_cycle_data$phase[g1sScore < g2mScore] <- "G2/M"
  cell_cycle_data$phase[g1sScore > g2mScore & g2mScore > threshold] <- "S"
  cell_cycle_data$phase <- factor(cell_cycle_data$phase, levels = c("Q", "G1", "S", "G2/M"))
  
  # Create scatter plot showing cell cycle scores
  a <- ggplot(cell_cycle_data) +
    geom_point(aes(x = G1S, y = G2M, color = cell_type, shape = phase), size = 2) +
    scale_colour_manual(values = color_palette) +
    xlab("G1/S phase score") + ylab("G2/M phase score") +
    ggtitle("Cell Cycle Analysis") +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_linerange(aes(x = threshold, ymin = 0, ymax = threshold)) +
    geom_path(data = data.frame(x = c(0, threshold, max(g1sScore)), y = c(threshold, threshold, max(g1sScore))),
              aes(x = x, y = y)) +
    geom_segment(aes(x = threshold, y = threshold, xend = max(g1sScore), yend = threshold))
  
  # Create stacked bar plot showing phase distribution
  b <- ggplot(cell_cycle_data %>%
                count(cell_type, phase) %>%
                group_by(cell_type) %>%
                arrange(desc(phase)) %>%
                mutate(pct = n / sum(n), ypos = cumsum(pct) - 0.5 * pct)) +
    geom_bar(aes(x = factor(cell_type), y = pct, fill = phase), stat = "identity", width = 0.8) +
    geom_text(aes(x = cell_type, y = ypos, label = paste0(n, "\n", sprintf("%1.1f", pct * 100), "%")), size = 2.5) +
    scale_fill_manual(values = c("darkseagreen", "plum2", "dodgerblue", "khaki2")) +
    ggtitle("Cell cycle Distribution of Cells") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
  
  # Save combined plot
  ggsave('./plot/cellcycle_tre.pdf', plot_grid(a, b), width = 12, height = 7)
  
  return(cell_cycle_data)
}

#' Advanced cell cycle analysis using expression matrix
#'
#' @param data_exprMat Expression matrix (genes x cells)
#' @param anno Annotation dataframe containing cell type info
#' @param col Color palette for cell types
#' @param thre Score threshold (NULL for auto)
#' @param org Organism ("hsa" for human, else mouse)
#' @param label Whether to show labels
#' @return DataFrame with cell cycle phase assignments
cell_cycle_analyze <- function(data_exprMat, anno, col, thre=NULL, org, label){
  
  # Subset expression matrix to match annotation
  data_exprMat <- data_exprMat[,rownames(anno)]
  
  # Define phase-specific gene markers
  if(org=="hsa"){
    g1sGene <- intersect(toupper(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")), rownames(data_exprMat))
    g2mGene <- intersect(toupper(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")), rownames(data_exprMat))
  } else {
    g1sGene <- intersect(c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"), rownames(data_exprMat))
    g2mGene <- intersect(c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"), rownames(data_exprMat))
  }
  
  # Calculate phase scores
  g1sScore <- colMeans(data_exprMat[g1sGene,], na.rm = T)
  g2mScore <- colMeans(data_exprMat[g2mGene,], na.rm = T)
  
  # Create dataframe for visualization
  cell_cycle_data <- data.frame(G1S = g1sScore, G2M = g2mScore, 
                               row.names = colnames(data_exprMat), 
                               cell_type = factor(anno$cluster, levels = label))
  
  # Set default threshold if not provided
  if(is.null(thre)) thre <- 3
  
  # Assign phases based on threshold
  cell_cycle_data$phase[cell_cycle_data$G1S < thre & cell_cycle_data$G2M < thre] <- "Q"
  cell_cycle_data$phase[!(cell_cycle_data$G1S < thre & cell_cycle_data$G2M < thre) & 
                        (cell_cycle_data$G1S < cell_cycle_data$G2M)] <- "G2/M"
  cell_cycle_data$phase[!(cell_cycle_data$G1S < thre & cell_cycle_data$G2M < thre) & 
                        (cell_cycle_data$G1S > cell_cycle_data$G2M) & 
                        cell_cycle_data$G2M < thre] <- "G1"
  cell_cycle_data$phase[!(cell_cycle_data$G1S < thre & cell_cycle_data$G2M < thre) & 
                        (cell_cycle_data$G1S > cell_cycle_data$G2M) & 
                        cell_cycle_data$G2M > thre] <- "S"
  
  cell_cycle_data$phase <- factor(cell_cycle_data$phase, levels = c("Q", "G1", "S", "G2/M"))
  
  # Calculate plot limits
  ccth <- thre
  ccx <- ceiling(max(cell_cycle_data$G1S))
  ccy <- ceiling(max(cell_cycle_data$G2M))
  
  # Create scatter plot
  a <- ggplot(cell_cycle_data) + 
    geom_point(mapping = aes(g1sScore, g2mScore, color = cell_type, shape=phase), size = 2) + 
    scale_colour_manual(values = col) +
    xlab("G1/S phase score") + ylab("G2/M phase score") + 
    ggtitle("Cell Cycle Analysis") +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, ccx)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, ccy)) +
    geom_linerange(mapping = aes(x = ccth, ymin = 0, ymax = ccth)) + 
    geom_path(data=data.frame(x = c(0, thre, ccx), y = c(thre, thre, ccx)), aes(x=x, y=y)) +
    geom_segment(mapping = aes(x = ccth, y = ccth, xend = ccx, yend = ccth))
  
  # Create stacked bar plot
  b <- ggplot(cell_cycle_data %>% 
                count(cell_type, phase) %>% 
                group_by(cell_type) %>% 
                arrange(desc(phase)) %>%
                mutate(pct = n/sum(n), ypos = cumsum(pct) - 0.5*pct)) + 
    geom_bar(mapping = aes(x = factor(cell_type), y = pct, fill = phase), 
             stat = "identity", width = 0.8) + 
    geom_text(mapping = aes(x = cell_type, y = ypos, 
                            label = paste0(n, "\n", sprintf("%1.1f", pct*100), "%")),
              size = 2.5) +
    scale_fill_manual(values = c("darkseagreen", "plum2", "dodgerblue", "khaki2")) +
    ggtitle("Cell cycle Distribution of Cells") + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  ggsave('./plot/cycle.pdf', plot_grid(a, b), width = 14, height = 6.5)
  return(cell_cycle_data)
}

#' Create stacked bar plot showing proportions
#'
#' @param meta_data DataFrame containing categorical variables to plot
#' @param group1 Primary grouping variable (x-axis)
#' @param group2 Secondary grouping variable (fill)
#' @param height Plot height
#' @param width Plot width
#' @param save_name Output file path
#' @param mycols Color palette
#' @param level1 Order of levels for group1
#' @param level2 Order of levels for group2
#' @param show_prob Whether to show percentage labels
plot_prop_barplot <- function(meta_data, group1, group2, height = NULL, width = NULL, 
                             save_name = './plot/prop_bar.pdf', mycols = NULL, 
                             level1 = NULL, level2 = NULL, show_prob = TRUE) {
  
  # Convert grouping variables to factors
  meta_data[[group1]] <- as.factor(meta_data[[group1]])
  meta_data[[group2]] <- as.factor(meta_data[[group2]])
  
  # Set level orders if specified
  if (!is.null(level1)) {
    meta_data[[group1]] <- factor(meta_data[[group1]], levels = level1)
  }
  if (!is.null(level2)) {
    meta_data[[group2]] <- factor(meta_data[[group2]], levels = level2)
  }
  
  # Calculate proportions
  cell_prop <- as.data.frame(prop.table(table(meta_data[[group1]], meta_data[[group2]])))
  colnames(cell_prop) <- c('cluster', group2, 'proportion')
  
  # Calculate total counts for labels
  total_counts <- cell_prop %>%
    group_by(!!sym(group2)) %>%
    summarise(total_number = sum(proportion)) %>%
    ungroup() %>%
    mutate(ratio = scales::percent(total_number))
  
  # Create base plot
  p <- ggplot(cell_prop, aes(x = !!sym(group2), y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    ggtitle("") +
    theme_bw() +
    theme(axis.ticks.length = unit(0.2, 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    guides(fill = guide_legend(title = NULL))
  
  # Add percentage labels if requested
  if (show_prob) {
    p <- p + geom_text(data = total_counts, 
                       aes(x = !!sym(group2), y = 1, label = ratio),
                       inherit.aes = FALSE, vjust = -0.2) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  }
  
  # Apply custom colors
  if (!is.null(mycols)) {
    mycols <- setNames(mycols, levels(meta_data[[group1]]))
    p <- p + scale_fill_manual(values = mycols)
  }
  
  # Save plot
  ggsave(filename = save_name, plot = p, width = width, height = height)
}

# Example usage with Seurat object
# seurat_obj <- readRDS('./data/fig1_seurat5_711.rds')
data_exprMat <- as.matrix(GetAssayData(seurat_obj, slot = 'counts'))
anno <- seurat_obj@meta.data[, c('celltype', 'barcode')]
anno$cluster <- as.character(anno$celltype)

# Run cell cycle analysis
cycle_df <- cell_cycle_analyze(data_exprMat, anno = anno, col = mycols, 
                              thre = 0.1, org = 'hsa', label = TRUE)

# Create proportion bar plot
plot_prop_barplot(meta_data = cycle_df, 
                 group1 = "phase", 
                 group2 = "cell_type", 
                 height = 6, 
                 width = 13, 
                 save_name = "./prop_bar.pdf", 
                 mycols = c("#99322e", "#f2a83b", "#cccd42", '#0000dc'), 
                 level1 = c('G2/M', 'S', 'G1', 'Q'), 
                 level2 = levels(cycle_df$cell_type), 
                 show_prob = FALSE)
```