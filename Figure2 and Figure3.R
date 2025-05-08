# code for Figure2 and Figure3
## slingshot

library(Seurat)
library(SCP)
library(ggplot2)

# Load Seurat object with slingshot trajectory
seurat_g1 <- readRDS('./data/seurat_g1_slingshot_629.rds')

# Run Slingshot for trajectory analysis
seurat_g1 <- RunSlingshot(srt = seurat_g1, group.by = "order_617", reduction = "umap", start = 'Inter.Epi')

# Plot trajectories on UMAP
p <- CellDimPlot(
  seurat_g1, 
  group.by = "order_617", 
  reduction = "umap", 
  lineages_span = 0.5, 
  lineages = paste0("Lineage", 1:3),
  palcolor = mycols[levels(seurat_my$order_617)],
  lineages_palcolor = c("#fe9929","#54278f", "#1c9099")
) 
ggsave('./plot/cs6_epi_slingshot_1.pdf', p, width = 7, height = 4.75)

# Plot lineage features on UMAP
p <- FeatureDimPlot(
  seurat_g1, 
  features = paste0("Lineage", 1:3), 
  reduction = "UMAP", 
  theme_use = "theme_blank"
)
ggsave('./plot/cs6_epi_slingshot_2.pdf', p, width = 17, height = 5)

# Identify dynamic features along trajectories
seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage1", "Lineage2", "Lineage3"), 
  n_candidates = 50,
  BPPARAM = BPPARAM
)

# Function to extract and sort top 50 genes
get_top_50_genes <- function(lineage_data) {
  # Sort by padjust value
  sorted_genes <- lineage_data[order(lineage_data$padjust), ]
  # Extract top 50 gene names
  top_50_genes <- rownames(sorted_genes)[1:50]
  return(top_50_genes)
}

# Run DynamicFeatures with custom gene sets
seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage1", "Lineage2", "Lineage3"),
  features = unique(amnion_markers)
)

seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage3"),
  features = unique(lineage3_genes),
  BPPARAM = BPPARAM
)

# Extract top dynamic feature genes for each lineage
lineage1_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage1$DynamicFeatures)
lineage2_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage2$DynamicFeatures)
lineage3_genes <- get_top_50_genes(seurat_g1@tools$DynamicFeatures_Lineage3$DynamicFeatures)
lineage3_genes <- c(lineage3_genes, 'CGB', 'OTX2', 'UCHL1', 'TBXT', 'POU5F1')

# Define gene sets for analysis
amnion_markers <- c("ISL1", "TFAP2B", "TFAP2A", "WNT6", "GABRP", "HEY1", "BAMBI", "DLX5", "SOX4", "GRHL1", 
                    "MEIS1", "BMP4", "PRTG", "DSP", "TGFBI", "MEST", "MSX2", "MITF", "VTCN1", "IGFBP3", 
                    "PRKD1", "KCNMA1", "STC1", "TCF4", "HAND1", "WNT11", "TBX3", "CDX2", "FURIN", "GATA6", 
                    "SALL4", "KRT19", "AQP3", "AMOTL1", "DAB2", "EMP2", "CA12", "CITED2", "TEAD1", "FOLR1", 
                    "PTGES")

lienage2_genes <- c(
  "CDX2", "HAND1", "PTN", "PGF", "CD55", "TBX3", "GATA6", "SALL4", "CITED2", "WNT11", "FURIN", 
  "KRT19", "AQP3", "AMOTL1", "DAB2", "EMP2", "CA12", "TEAD1", "FOLR1", "PTGES", "GABRP", "VTCN1", 
  "TFAP2A", "TFAP2B", "IGFBP3", "PRKD1", "KCNMA1", "STC1", "TCF4", "ISL1", "DLX5", "PRTG", "TGFB1", 
  "WNT6", "BMP4", "MSX2", "DSP", "MEIS1", "SOX4", "BAMBI", "HEY1", "MEST", "GRHL1", "MITF"
)

EMT_genes <- c("AXL", "BMI1", "CDH1", "CDKN2A", "E2-2", "E47", "EHMT2", "EPAS1", "ETS1", "EZH2", 
               "FOXC2", "G9A", "GSC", "HDAC1", "HDAC2", "HDAC3", "HIF1A", "KLF8", "LOX", "LOXL2", 
               "LSD1", "MCT1", "PIK3CA", "PRRX1", "PTEN", "RB1", "SIP1", "SIRT1", "SIX1", "SLUG", 
               "SMAD2", "SNAI1", "SNAI2", "SNAIL", "SNAIL1", "SNAIL2", "SUV39H1", "SUZ12", "T", 
               "TCF4", "TP53", "TWIST", "TWIST1", "ZEB1", "ZEB2")

# PRCP genes analysis
prcp_genes <- c("GDF3", "WNT3A", "FGF17", "GSC", "SAT1", "LHX1", "CHRD", "POU5F1", "NOTO", "HHEX", 
                "NODAL", "FOXA2", "MIXL1", "SP5", "TBXT", "SHISA2", "CDH2", "FZD8", 
                "MT1H", "SIX3", "OTX2", "SOX17", "SFRP1", "CER1", "DKK1", "MESP1", "FOXJ1", "DKK4")

prcp_genes <- intersect(prcp_genes, rownames(seurat_g1))

library(BiocParallel)
# Set BPPARAM for single-thread processing
BPPARAM <- SerialParam()

seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage1"),
  features = unique(prcp_genes),
  BPPARAM = BPPARAM
) 

seurat_l1 <- subset(seurat_g1, Lineage1 >= 0)

# Plot dynamic features for PRCP lineage
p <- DynamicPlot(
    srt = seurat_l1, 
    lineages = c("Lineage1"), 
    group.by = "celltype", # Don't use 'order' as it may cause errors
    features = prcp_genes,
    compare_lineages = TRUE, 
    compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#fe9929",  # Line color: #1c9099 for Lineage3, #54278f for Lineage2, #fe9929 for Lineage1
    point_palcolor = list(order_617 = mycols[levels(seurat_l1$order_617)][c('Inter.Epi', 'Anterior.Epi', 'PrCP')]),
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
ggsave('./plot/lingage1_genes_prcp_1.pdf', p, width = 24, height = 12)

# EMT genes analysis
seurat_g1 <- RunDynamicFeatures(
  srt = seurat_g1, 
  lineages = c("Lineage3"),
  features = unique(emt_genes),
  BPPARAM = BPPARAM
) 

p <- DynamicPlot(
    srt = seurat_g1, 
    lineages = c("Lineage3"), 
    group.by = "order_617",
    features = emt_genes,
    compare_lineages = TRUE, 
    compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#1c9099",  # Line color: #1c9099 for Lineage3
    point_palcolor = mycols[levels(seurat_g1$order_617)][c('PrCP', 'Anterior.Epi', 'Inter.Epi', 'Posterior.Epi', 'AM.Ecto', 'AM')],
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
ggsave('./plot/lingage3_genes_emt_test.pdf', p, width = 24, height = 14)

# Combine gene sets for Lineage2 analysis
genes <- unique(c(lienage2_genes, amnion_markers))

seurat_g1_s1 <- RunDynamicFeatures(
  srt = seurat_g1_s1, 
  lineages = c("Lineage2"),
  features = unique(genes),
  BPPARAM = BPPARAM
)

# Plot dynamic features for Lineage2
p <- DynamicPlot(
    srt = seurat_g1_s1, 
    lineages = c("Lineage2"), 
    group.by = "celltype",
    features = genes,
    compare_lineages = TRUE, 
    compare_features = FALSE,
    pt.size = 0.5,
    line_palcolor = "#54278f",  # Line color: #54278f for Lineage2
    point_palcolor = mycols[levels(seurat_g1$order_617)][c('PrCP', 'Anterior.Epi', 'Inter.Epi', 'Posterior.Epi', 'AM.Ecto', 'AM')],
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
ggsave('./plot/lingage2_genes_2.pdf', p, width = 24, height = 17)


## pyScenic

### create_loom_input.R

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



### pyscenic_from_loom.sh

``` sh
input_loom=out.loom
n_workers=20
#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-h|--help:\t Usage information"
exit 1
}
#get value
while getopts :i:n:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done
#database path
tfs=/home/xlyang/python_work/pyscenic_pipeline/01.database/human_base/hs_hgnc_tfs.txt
feather=/home/xlyang/python_work/pyscenic_pipeline/01.database/human_base/*.feather
tbl=/home/xlyang/python_work/pyscenic_pipeline/01.database/human_base/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
pyscenic=/sdc/xlyang/software/anaconda3/envs/pyscenic/bin/pyscenic

# grn
 $pyscenic grn \
 --num_workers $n_workers \
 --output grn.tsv \
 --method grnboost2 \
 $input_loom  $tfs

# cistarget
$pyscenic ctx \
grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
ctx.csv \
--output aucell.loom \
--num_workers $n_workers
```

### calcRSS_by_scenic.R

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

### Execute the above script to run pyscenic

``` R
inscp='../00.scripts/'
inrds='./seurat_g1_slingshot_629.rds'
ingroup='celltype'

# Step 1: Activate snapatacv1 environment and run the first R script
source /home/xlyang/software/anaconda3/bin/activate snapatacv1
Rscript ${inscp}/create_loom_input.R -i ${inrds} -d ${ingroup} -l out -a Spatial
sleep 30

# Step 2: Activate pyscenic environment and run the shell script
source /home/xlyang/software/anaconda3/bin/activate pyscenic
sh ${inscp}/pyscenic_from_loom.sh -i out.loom -n 10

# Step 3: Reactivate snapatacv1 environment and run the second R script
source /home/xlyang/software/anaconda3/bin/activate snapatacv1
Rscript ${inscp}/calcRSS_by_scenic.R -l aucell.loom -m metadata_subset.xls -c ${ingroup} -a Spatial

# Optionally, deactivate the environment at the end of the script
conda deactivate
```

### Visualize the RSS

``` R
library(Seurat)
library(pheatmap)

df <- readRDS('../celltype2_rss.rds')

# Define Min-Max normalization function
min_max_normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

# Perform Min-Max normalization for each row (each gene)
#htdf <- t(apply(df, 1, min_max_normalize))
htdf <- scale(df, scale = TRUE, center = TRUE)
# Calculate row standard deviation and filter
row_sds <- apply(htdf, 1, sd)
filtered_mat <- htdf[row_sds > 0.95, ]
htdf <- filtered_mat

# Create column annotation
Cell_Type <- colnames(df)
annotation_col_df <- data.frame(Cell_Type = Cell_Type)
rownames(annotation_col_df) <- Cell_Type

# Define column order
#custom_order <- c("CTB_pri", "CTB_ccc", "iEVT_naïve", "iEVT_mature", "GC")
htdf <- htdf[, custom_order]
annotation_col_df <- annotation_col_df[custom_order, , drop = FALSE]

# Define colors for cell types
#mycols <- c('CTB_pri'='#35648F', 'CTB_ccc'='#BF3F45', 'iEVT_naïve'='#E29F46', 'iEVT_mature'='#BAC76B', 'GC'='#8C29C9')
colors_list <- list('Cell_Type' = mycols)

# Create heatmap
p <- pheatmap(htdf, 
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              annotation_col = annotation_col_df,
              annotation_colors = colors_list,
              show_rownames = TRUE,
              show_colnames = TRUE,
              color = colorRampPalette(c('#253494', '#2c7fb8', '#c7e9b4', '#ffffcc'))(100),
              main = "RSS Score (Z-Score Normalized)")

print(p)
pdf('./circular_heatmap.pdf',width = 5,height = 7)
print(p)
dev.off()

library(dendextend) # Add this line
# Load necessary packages
library(circlize)
library(ComplexHeatmap)

# Define color mapping
mycol2 <- colorRamp2(c(-1.7, 0.3, 2.3), c("#57ab81", "white", "#ff9600"))

# Draw circular heatmap
circos.clear()
circos.par(gap.after = c(22))
circos.heatmap(hfdf, col = mycol2, dend.side = "inside", rownames.side = "outside", 
               track.height = 0.38, rownames.col = "black", rownames.cex = 0.9, 
               rownames.font = 1, cluster = TRUE, dend.track.height = 0.18,
               dend.callback = function(dend, m, si) {
                 color_branches(dend, k = length(mycols), col = mycols)
               })

# Add legend
lg <- Legend(title = "Exp", col_fun = mycol2, direction = "vertical")
grid.draw(lg)

# Add column names
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if (CELL_META$sector.numeric.index == 1) {
    cn <- colnames(hfdf)
    n <- length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.8, "mm"),
                7.8 + (1:n) * 1.1,
                cn, cex = 0.8, adj = c(0, 1), facing = "inside")
  }
}, bg.border = NA)

circos.clear()
```



## FASTMNN

``` R
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(SeuratDisk)

# Load Seurat object and metadata
seurat_obj <- LoadH5Seurat('./data/cs6_all.h5seurat')
meta_cs6 <- readRDS('./data/CS6_all_meta_250402.rds')
head(meta_cs6)

# Add metadata and subset cells
seurat_obj <- AddMetaData(seurat_obj, meta_cs6)
unique(seurat_obj@meta.data$celltype_new)
seurat_obj <- subset(seurat_obj, celltype_new %in% c('Posterior.Epi', 'Hypo.2', 'Hypo.1', 'ave', 'Anterior.Epi', 'Inter.Epi', 'Anterior pole'))

# Set RNA assay as default and remove Spatial assay
seurat_obj[['RNA']] = seurat_obj[['Spatial']]
DefaultAssay(seurat_obj) = 'RNA'
seurat_obj[['Spatial']] = NULL
seurat_obj$Day = 'CS6'
seurat_obj$celltype = seurat_obj$celltype_new

# Load and process LTQ data
seurat_ltq <- readRDS('/home/xlyang/download/refer_data/ltx_14d/ltq_14.rds')
unique(seurat_ltq@meta.data$Day)
unique(seurat_ltq@meta.data$Group)
seurat_ltq$celltype <- seurat_ltq@meta.data$Group
seurat_ltq <- subset(seurat_ltq, Day %in% c('D12', 'D14') & Group %in% c('EPI', 'PrE', 'PSA-EPI'))

# Define color palettes
ltq_col <- c('EPI' = '#813a90', 'PSA-EPI' = '#904410', 'PrE' = '#62509d')
CS6_color <- c('Anterior pole' = '#E86c79',  # "PrCP"
               "Anterior.Epi" = '#253494',
               "Inter.Epi" = '#006837', "Posterior.Epi" = '#1d91c0',
               "AM.Ecto" = '#4682B4', "AM" = '#bcbddc',
               "AM.EXMC" = '#9e9ac8', "Connecting Stalk" = '#fdd0a2',
               "Hypo.1" = '#6a51a3', "Hypo.2" = '#ae017e',
               "PL.EXMC" = '#ccebc5', "CTB" = '#2b8cbe',
               "CTB.Fusion" = '#7bccc4', "STB" = '#58BCE8',
               "MTB" = '#084081', "YS.Endo_1" = '#FFAA92',
               "YS.Endo_2" = '#8FB0FF', "YS.Endo_3" = '#00C2A0',
               "YS.EXMC_1" = '#6F0062', "YS.EXMC_2" = '#EEC3FF',
               "Myeloid Progenitor" = '#D16100', "Primitive_Ery1" = 'magenta',
               "Primitive_Ery2" = '#B79762', "Primitive_Mk1" = 'red3',
               "Primitive_Mk2" = 'purple', "YS.Endo" = '#ec7014',
               "YS.EXMC" = '#993404', "Blood" = '#e31a1c',
               "grey" = '#EBEBEB', "Epi" = '#006837', "Hypo" = '#6a51a3',
               'ave' = 'green')
mycols <- c(ltq_col, CS6_color)

# Prepare datasets for integration
seuratcs6 = seurat_obj
seuratcs6 <- subset(seuratcs6, features = rownames(seurat_ltq))
seurat_ltq <- subset(seurat_ltq, features = rownames(seuratcs6))
seurat_ltq$sample <- 'E12-E14'
seuratcs6$sample <- 'cs6'

# Merge datasets and prepare for integration
seurat_obj = merge(seuratcs6, seurat_ltq)
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")
seurat_obj <- JoinLayers(seurat_obj)

# Normalize data and select integration features
seurat_obj <- NormalizeData(seurat_obj)
features <- SelectIntegrationFeatures(object.list = SplitObject(seurat_obj, 'sample'))
VariableFeatures(seurat_obj) = features

# Split objects for integration
split_objects <- SplitObject(seurat_obj, split.by = "sample")
print(names(split_objects))

# Reorder objects to ensure cs6 is first
split_objects <- split_objects[c("cs6", "E12-E14")]

# Run FastMNN integration
seurat_obj <- RunFastMNN(
  object.list = split_objects,
  k = 8,
  merge.order = c("cs6", "E12-E14"),
  auto.merge = FALSE
)

# Run UMAP and visualize results
seurat_obj <- RunUMAP(seurat_obj, reduction = "mnn", dims = 1:10)

library(scplotter)
p1 = CellDimPlot(
  seurat_obj,
  group_by = "celltype",
  reduction = "umap",
  label = TRUE,
  palcolor = mycols,
  pt.size = 1.5,
  sizes.highlight = 2,
  theme = "theme_blank",
  legend.position = "right",
  raster = FALSE,
  highlight = 'sample == "E12-E14"'
)

# Save plot
library(ggplot2)
ggsave('./plot/e1214_cs6_fmnn_p1.pdf', p1, width = 6.5, height = 5.75)
```

## dotplot

``` R
rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)

#### Data reading
cs6_raw <- readRDS("D:/文件/bioinfo/DATA/人CS6data＆mk＆ms/CS6_velocity_2.rds")
# Convert cs6 assay to 'RNA'
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

# Rename required cell populations and set order
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
AVE_Nature$ig_celltype <- factor(AVE_Nature$ig_celltype, levels = c("YS Endoderm_cs7", 'DE(P)_cs7',"Hypoblast_cs7",'DE(NP)_cs7'))

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

# Get human-mouse gene mapping table
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
  cat("Duplicate gene names:", duplicated_genes, "\n")
  
  # Merge expression of duplicate genes (can choose sum, mean, etc.)
  expression_matrix <- aggregate(expression_matrix, 
                                 by = list(rownames(expression_matrix)), 
                                 FUN = sum)  # Here we sum expression of duplicate genes, or you can choose mean etc.
  
  # Set merged gene names
  rownames(expression_matrix) <- expression_matrix$Group.1
  expression_matrix <- expression_matrix[, -1]  # Remove Group.1 column
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

##### Data integration
dtlist <- list(cs6, cs7_ncb, cs8, AVE_Nature, mk_ve, ms_ve_trans)
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

############ Dimensionality reduction and visualization
dt_integrated <- FindVariableFeatures(dt_integrated, selection.method = "vst", nfeatures = 2000)
dt_integrated <- NormalizeData(dt_integrated, normalization.method = 'LogNormalize', scale.factor = 10000)
# Scale data
dt_integrated <- ScaleData(dt_integrated, features = head(VariableFeatures(dt_integrated), 1000))
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
                  'Anterior Pole' = c('POU5F1','CHRD','GSC', 'DDIT4','CRIP1', 'FGF17',
                                    'FOXA2','SHH','EOMES','NODAL','NOTO',  'FGF4', 
                                    'SNRPN','MT1X','SIX3','CYP26A1','LEFTY1'))

# Draw dotplot
p <- DotPlot(dt_integrated, features = rev(list_genes), group.by = 'ig_celltype', scale = TRUE) +
  theme(panel.grid.major = element_line(color = "#E4EEED", size = 0.5),  # Major grid lines
        panel.grid.minor = element_line(color = "#E4EEED", size = 0.2),  # Minor grid lines
        panel.background = element_rect(fill = "white", color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(color = "black"),
         panel.spacing = unit(1, "mm"),
         axis.title = element_blank(),
  ) +
  scale_color_gradientn(values = seq(0, 1, length.out = 7), 
                        # colours = c('#f6f6f6','#f5f2f3','#f0ebeb','#F1E6E6','#f8e7eb','#dc8a9a','#be3b56','#d15c6c','#8d192b')) +
                        colours = c('#ACC7DC','#FAF8F4','#ECCDD2', '#DEA2AF')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
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
```

