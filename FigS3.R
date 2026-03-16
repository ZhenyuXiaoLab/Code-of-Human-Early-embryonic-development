#code for Fig.S3
# key Software versions:
# clusterProfiler 4.1.2
# ggplot2 3.4.4
# igraph 2.0.3
# SCP 0.5.3
### GO netplot by igraph|SCP

####---- Load Packages ----####
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(cowplot)
library(SCP)
library(dplyr)
library(BiocParallel)

####---- Load Data ----####
GO_result <- read_delim(file = "./table/blood_go/GO_result/blood_sl_up.csv", delim = ",", col_names = T)

GO_result_edge <- GO_result %>%
  dplyr::select(2, 3, 11) %>%
  dplyr::mutate(ID2 = str_c(ONTOLOGY, ID, sep = "_")) %>%
  dplyr::mutate(Type = "GO") %>%
  dplyr::select(Type, everything())

# Node file
node_file1 <- GO_result_edge %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(Type_sum = sum(Count)) %>%
  purrr::set_names(c("node", "node_size")) %>%
  dplyr::mutate(node_level = "Type") %>%
  dplyr::mutate(type = node)

node_file2 <- GO_result_edge %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::summarise(ONTOLOGY_sum = sum(Count)) %>%
  purrr::set_names(c("node", "node_size")) %>%
  dplyr::mutate(node_level = "ONTOLOGY") %>%
  dplyr::mutate(type = node)

node_file3 <- GO_result_edge %>%
  dplyr::group_by(ID2) %>%
  dplyr::summarise(ID2_sum = Count) %>%
  purrr::set_names(c("node", "node_size")) %>%
  dplyr::mutate(node_level = "ID") %>%
  dplyr::mutate(type = str_remove(string = node, pattern = "_.*")) %>%
  dplyr::mutate(node = str_remove(string = node, pattern = "^.*_"))

node_file <- rbind(node_file1, node_file2, node_file3)

# Edge file
edge_file1 <- GO_result_edge %>%
  dplyr::group_by(Type, ONTOLOGY) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(1, 2) %>%
  purrr::set_names(c("from", "to"))

edge_file2 <- GO_result_edge %>%
  dplyr::group_by(ONTOLOGY, ID) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(1, 2) %>%
  purrr::set_names(c("from", "to"))

edge_file <- rbind(edge_file1, edge_file2)

graph_country <- tbl_graph(node_file, edge_file)

####---- Plot GO Graph ----####
ggraph(graph_country, layout = 'dendrogram', circular = TRUE) +
  geom_edge_diagonal(aes(color = node1.node), alpha = 1/3) +
  geom_node_point(aes(size = node_size, color = type), alpha = 1/3) +
  coord_fixed() +
  scale_size(range = c(3, 15)) +
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = node,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = type
    ),
    size = 2, hjust = 'outward'
  ) +
  geom_node_text(
    aes(label = node,
        filter = !leaf,
        color = type),
    fontface = "bold",
    size = 3,
    family = "sans"
  ) +
  theme_nothing() +
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

ggsave(filename = "GO_graph.pdf", height = 8, width = 8)

####---- SCP Enrichment Analysis ----####
blood_markers <- read.csv('./table/blood_go/BLOOD_deg.csv')

blood_up <- blood_markers %>%
  filter(avg_log2FC > 0) %>%
  pull(gene)

blood_down <- blood_markers %>%
  filter(avg_log2FC < 0) %>%
  pull(gene)

all_genes <- c(blood_up, blood_down)
gene_groups <- factor(c(rep("UP", length(blood_up)),
                        rep("DOWN", length(blood_down))))

# Perform enrichment analysis
enrichment_results <- RunEnrichment(
  srt = NULL,
  geneID = blood_up,
  db = c("GO_BP"),
  minGSSize = 10,
  maxGSSize = 500,
  GO_simplify = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)

# Load seurat object for grouped analysis
pancreas_sub <- RunEnrichment(
  srt = seurat_obj, group_by = "celltype", db = c("GO", "KEGG"),
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05",
  BPPARAM = BiocParallel::SerialParam()
)

# Plot enrichment map
p <- EnrichmentPlot(
  srt = pancreas_sub, group_by = "celltype", group_use = "Blood",
  plot_type = "enrichmap", db = c("GO"), topTerm = 160
)

pdf('./plot/Blood_go_clustered.pdf', width = 11, height = 11)
print(p)
dev.off()
```

####The clustering and bubble plot generation code uses the same implementation as shown in FIG1.R.
#and is not repeated here.
#3D modeling was accomplished using the modeling software MAYA.