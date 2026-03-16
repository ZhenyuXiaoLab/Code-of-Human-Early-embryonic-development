## code for Fig.S9
#Seurat 5.0.3
#ggplot2 3.4.4

library(ggplot2)
library(Seurat)

# Load preprocessed single-cell RNA-seq data
CS6_epi_only <- readRDS("D:/data/bioinfo/DATA/humanCS6data＆mk＆ms/CS6_epi_only_cs.rds")

# Find variable features
CS6_epi_only <- FindVariableFeatures(CS6_epi_only, selection.method = "vst", nfeatures = 2000)
# Normalize data using LogNormalize method with scale factor 10000
CS6_epi_only <- NormalizeData(CS6_epi_only, normalization.method = 'LogNormalize', scale.factor = 10000)
# Scale the data
CS6_epi_only <- ScaleData(CS6_epi_only, features = VariableFeatures(CS6_epi_only))
# Run PCA dimensionality reduction
CS6_epi_only <- RunPCA(CS6_epi_only, features = VariableFeatures(CS6_epi_only))
# Set cell type factor levels in reverse order
CS6_epi_only$celltype_new <- factor(CS6_epi_only$celltype_new, levels = rev(c('Anterior pole',
                                                                             'Anterior.Epi',
                                                                             'Inter.Epi',
                                                                             'Posterior.Epi')))

# Define FGF pathway genes
FGF_gene <- list(
  ligand = c('FGF1','FGF17','FGF4','FGF5', 'FGF8', 'FGF2'), 
  receptor = c('FGFR4', 'FGFR3', 'FGFR1', 'FGFR2'),         # Corrected spelling
  transcription = c('ETS2', 'ETS1')                         # Keep original category name
)

# Generate dot plot for FGF genes
p1 <- DotPlot(CS6_epi_only, features = FGF_gene, group.by = 'celltype_new', scale = TRUE) +
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
                        colours = c('#CCFDFE', '#CDD6DD', '#CDACA1', '#E29693')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_size_continuous(range = c(2, 6))

p1

# Define WNT pathway genes
WNT_gene <- list(
  ligand = c('WNT3A','WNT9B','WNT2', 'WNT2B',   
             'WNT5A', 'WNT5B', 'WNT6',  
             'WNT11'),  # Note: Wnt3 in mouse corresponds to human WNT3
  receptor = c('FZD1', 'FZD2', 'FZD3', 'FZD4', 'FZD5',  
               'FZD7', 'FZD8', 'ROR2', 'ROR1', 'RYK', 'LRP6', 
               'PTK7','FZD6', 'LRP5'),  # Core receptors
  co_receptor = c('SDC4','GPC4','GPC2',  'SDC2','GPC1', 'SDC3','SDC1','GPC6',    # Syndecan family
                  'GPC3'),  # Glypican family
  inhibitor = c('CER1','DKK1', 'WIF1', 'SFRP2', 'SFRP1', 'SFRP5', 
                'WISE'),  # Sclerostin/Wise is SOST gene product
  agonist = c('RSPO4', 'RSPO2', 'RSPO1', 'RSPO3'  # R-spondin family
  ),  # Norrin disease protein
  co_protein = c('LGR5',   # LGR receptor family
                 'ZNRF3','LGR4')  # E3 ubiquitin ligase
)

# Generate dot plot for WNT genes
p2 <- DotPlot(CS6_epi_only, features = WNT_gene, group.by = 'celltype_new', scale = TRUE) +
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
                        colours = c('#CCFDFE', '#CDD6DD', '#CDACA1', '#E29693')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_size_continuous(range = c(2, 6))

p2

# Define NODAL pathway genes
NODAL_gene <- list(
  ligand = c('NODAL'),                        # TGF-β superfamily core ligand (OMIM: 601265)
  co_factor = c('TDGF1'),                     # Co-receptor (TDGF1 also known as CRIPTO)
  inhibitor = c('LEFTY1', 'LEFTY2', 'CER1'),  # Secreted antagonists
  receptor = c('ACVR1B','ACVR2A', 'ACVR2B' ), # Receptor complex (ACVR1B is ALK4)
  transcription_factor = c('SMAD3', 'SMAD2', 'SMAD4') # R-Smad/Co-Smad system
)

# Generate dot plot for NODAL genes
p3 <- DotPlot(CS6_epi_only, features = NODAL_gene, group.by = 'celltype_new', scale = TRUE) +
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
                        colours = c('#CCFDFE', '#CDD6DD', '#CDACA1', '#E29693')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_size_continuous(range = c(2, 6))

p3

# Define BMP pathway genes
BMP_gene <- list(
  ligand = c('BMP7','BMP4','BMP5','GDF6'), # TGF-β superfamily bone morphogenetic protein subfamily
  receptor = c(
    'BMPR2',    # BMPR-II (OMIM: 600799)
    'ACVR2A',   # ActR-II (OMIM: 102581)
    'ACVR2B',
    'BMPR1A',   # ALK-6 (OMIM: 603248)
    'ACVR1'     # ALK-2 (OMIM: 102576)
  ),
  co_factor = c('LMK1'), # Auxiliary regulatory factors (need to verify latest nomenclature)
  inhibitor = c(
    'CHRD','SKI','BTG2','DACH1',   # Receptor tyrosine kinases (possible cross-regulation)
    'BTG1',
    'NOG',      # B-cell translocation gene family
    'SKIL','BTG3' # Transcriptional repressor complex members
  ),
  transcription_factor = c(
    'SMAD5','CREBBP',   # BMP-specific R-Smads
    'SMAD4','KAT2A',    # Co-Smad (OMIM: 600993)
    'YY1','EP300','MSX1','ZEB2','SMAD8' # Co-activators/chromatin regulators
  ),
  target_genes = c(
    'ID3','TCF4','CREB3L1', 'PRRX2','TGFB1I1','TCF7','HEY1','ID1','ID2','ID4','GATA3',   # DNA-binding inhibitor factors
    'GATA2','SNAI1','SNAI2',
    'TBX2'
  )
)

# Generate dot plot for BMP genes
p4 <- DotPlot(CS6_epi_only, features = BMP_gene, group.by = 'celltype_new', scale = TRUE) +
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
                        colours = c('#CCFDFE', '#CDD6DD', '#CDACA1', '#E29693')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_size_continuous(range = c(2, 6))

p4

# Define TGF-β pathway genes
TGFb_gene <- list(
  ligand = c(
    'TGFB1',   # TGF-β1 (OMIM: 190180)
    'TGFB3' 
  ),
  receptor = c(
    'F2RL1',   # PAR-1 (thrombin receptor, OMIM: 187930)
    'TGFBR3',  # TGF-β type II receptor (OMIM: 190182)
    'TGFBR1',  # β-glycan (auxiliary receptor, OMIM: 600742)
    'F2R'      # PAR-2 (OMIM: 600934)
  ),
  transcription_factor = c(
    'GSC',     # Goosecoid (embryonic dorsal lip marker)
    'OTX2',
    'HHEX',
    'GRB2',
    'SMAD2',   # R-Smad (OMIM: 603109)
    'SMAD7'    # I-Smad (inhibitory, OMIM: 602932)
  )
)

# Generate dot plot for TGF-β genes
p5 <- DotPlot(CS6_epi_only, features = TGFb_gene, group.by = 'celltype_new', scale = TRUE) +
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
                        colours = c('#CCFDFE', '#CDD6DD', '#CDACA1', '#E29693')) +
  labs(x = NULL) +
  guides(size = guide_legend(order = 3),
         color = guide_legend(title = "avg.exp")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_size_continuous(range = c(2, 6))

p5
