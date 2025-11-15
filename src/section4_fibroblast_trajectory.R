library(reticulate)
use_condaenv('py311')
py_config()

#install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
#install.packages("BiocManager")
#BiocManager::install(c("GSVA", "GSEABase"))
# 
# if (!requireNamespace("zellkonverter", quietly = TRUE)) {
#   BiocManager::install("zellkonverter")
# }
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterExperiment")
# 

# Bioconductor packages
library(ComplexHeatmap)
library(tradeSeq)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)

# CRAN packages
library(pheatmap)
library(data.table)
library(slingshot)
library(tidyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(ggplot2)
library(circlize)
library(grid)
library(tibble)

library(SingleCellExperiment)
library(zellkonverter)


gene_groups <- list(
  VECs = c("IFI27", "SOX17", "SPRY1", "TSC22D1", "CTNNB1", "ICAM1", "VCAM1"),
  LECs = c("ECSCR", "FABP4", "LINC00636", "LYVE1", "MMRN1", "PROX1", "TFPI"),
  ECs = c("AQP1", "CDH5", "CLDN5", "GNG11", "IGFBP7", "MYCT1", "PLVAP", "RAMP2",
          "RNASE1", "SELE", "SPARCL1", "TFF3", "TM4SF1", "VWF"),
  Fbs = c("ADAM12", "APCDD1", "APOD", "ASPN", "CALD1", "COCH", "COL5A2", "COL6A1", "COL6A2",
          "COL6A3", "CXCL12", "FBLN1", "GEM", "GPR4", "HAPLN1", "HTRA1", "IGFBP5", "LEPR",
          "LSAMP", "LUM", "MFAP4", "MFAP5", "MGP", "MMP2", "MMP27", "MYL9", "NDUFA4L2",
          "PDGFRA", "PLAC9", "POSTN", "PTGDS", "SFRP2", "SLPI", "SOD3", "TAGLN", "THY1", "VIM"),
  Mph = c("AIF1", "BASP1", "C15orf48", "C1QA", "C5AR1", "CD68", "CD83", "CD93", "COTL1",
          "CTSZ", "FCER1G", "GPR183", "HMGB2", "HMGN2", "ID2", "IL1B", "INHBA", "INSIG1",
          "IRF4", "LYZ", "RGS1", "STMN1", "TSPAN33", "TYROBP"),
  T_cell = c("ACER1", "ALOX5AP", "ARHGDIB", "BHLHE41", "CD3D", "CD3E", "CD3G", "CD40LG", "CD52",
             "CD69", "CD8A", "CES1", "CXCR4", "DUSP2", "GATA2", "GIMAP7", "GZMK", "IFNG", "IL2RA",
             "IL32", "ITM2A", "KIT", "KLRB1", "LCK", "LENG8", "MKI67", "MT1F", "NLGN4Y", "NMB",
             "NR4A2", "PBXIP1", "PTPRCAP", "RASSF8", "RGS5", "RRM2", "SLC2A3", "UBE2C"),
  DCs = c("AXL", "BIRC3", "C1orf54", "CALCRL", "CBFA2T3", "CCL19", "CCL22", "CCR7", "CD1A",
          "CD1B", "CLDN1", "CLEC10A", "CLEC9A", "CPVL", "ENPP1", "FCER1A", "FSCN1", "GPR157",
          "IDO1", "IL3RA", "IL4I1", "IRF8", "LGALS2", "LST1", "LTB", "MARCKSL1", "MNDA",
          "PKIB", "PTGER3", "RGCC", "S100B", "SLC8A1", "SYNPO2", "TFPI2", "TMEM150C",
          "TMEM176A", "TUBB2B", "WDFY4")
)



# Define custom color mappings
stage_colors <- c(
  'nodular' = '#ed322f',
  'plaque' = '#ffdc5e',
  'patch' = '#51f512',
  'control' = 'blue'
)

broad_cell_types_colors <- c(
  'Lymphatic Endothelial Cells' = '#ffb695',
  'Macrophages' = '#ff40ff',
  'Vascular Endothelial Cells' = '#a4e000',
  'Fibroblasts' = '#c7d0c0',
  'T-cells' = '#941100',
  'Dendritic cells' = '#ff9300',
  'Macrophages_M1-like' = '#611861',
  'Macrophages_M2-like' = '#ff40ff'
)


infection_status_colors <- c(
  'uninfected' = '#d3d3d3',
  'infected' = 'red'
)

niche_colors <- c(
  "Basal Dermis"= "#b299e3",
  "Differentiated Epidermis"= "#ffd000",
  "Stroma"= "#646500",
  "TA VEC Stroma"= "#00e50c",
  "VEC Stroma"= "#cccc33",
  "Macrophage Immune Stroma"= "#00dbf4",
  "T-cell Immune Stroma"= "#0051f9",
  "Immune"= "#c100f9",
  "Tumor Core"= "#450000",
  "Tumor"= "#eb0000",
  "Tumor Boundary"= "#faa0aa"
)

leiden_colors <- c(
  "0" = "blue",
  "1" = "orange",
  "2" = "green",
  "3" = "red",
  "4" = 'purple'
)

group_colors <- c(
  VECs = "#a4e000",
  LECs = "#ffb695",
  ECs = "red",
  Fbs = "#c7d0c0",
  Mph = "#ff40ff",
  T_cell = "#941100",
  DCs = "#ff9300"
)


# Custom cell type colors
sub_cell_types_color_mapping <- c(
  'Keratinocytes' = '#181c82',
  'Differentiated Keratinocytes' = '#471fc7',
  'Spinous to Granular Cells' = '#034cff',
  'Pilosebaceous Cells' = '#bbbde2',
  'Melanocytes' = '#00bbbf',
  'Vascular Endothelial Cells' = '#a4e000',
  'Lymphatic Endothelial Cells' = '#ffb695',
  'Proliferating Lymphatic Endothelial Cells' = '#906855',
  'Pericytes' = '#9f7704',
  'Fibroblasts' = '#c7d0c0',
  'Pro-inflammatory Fibroblasts' = '#7cd28e',
  'Mesenchymal Fibroblasts' = '#3d8e27',
  'Myofibroblasts' = '#007c1d',
  'Secretory-papillary Fibroblasts' = '#076018',
  'Secretory-reticular Fibroblasts' = '#324708',
  'Macrophages' = '#ff40ff',
  'Dendritic cells' = '#ff9300',
  'B-cells' = '#f12d00',
  'T-cells' = '#941100',
  'Cd4' = '#600c09',
  'Cd4 Rgcc' = '#3c0c09',
  'Cd8 Exhausted' = '#240c09'
)

setwd("C:/Users/ard212.PITT/development/KS_xenium")

# Increase allowed memory size to 2GB (adjust if needed)
options(future.globals.maxSize = 2 * 1024^3)  # 2GB


# output_prefix = 'with_macrophages'
output_prefix = 'fibroblasts_and_LEC'


dir.create(paste0("results/pseudotime_analysis/", output_prefix, "/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("figures/pseudotime_analysis/", output_prefix, "/"), recursive = TRUE, showWarnings = FALSE)



# Construct file paths dynamically
expression_matrix <- as.matrix(read.csv(paste0("data/expression_matrix_subset_", output_prefix, ".csv"), row.names=1))
cell_metadata <- read.csv(paste0("data/cell_metadata_subset_", output_prefix, ".csv"), row.names=1)
gene_metadata <- read.csv(paste0("data/gene_metadata_subset_", output_prefix, ".csv"), row.names=1)
gene_names <- read.csv(paste0("data/gene_names_subset_", output_prefix, ".csv"), header=FALSE, stringsAsFactors=FALSE)[[1]]
phate_coords <- read.csv(paste0("data/phate_embedding_", output_prefix, ".csv"), row.names=1)
leiden_clusters <- read.csv(paste0("data/leiden_cluster_ids_", output_prefix, ".csv"), row.names=1)


expression_matrix <- t(expression_matrix)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = expression_matrix,
  assay = "RNA", 
  meta.data = cell_metadata
)

seurat_obj@meta.data = cell_metadata
seurat_obj

# Load PHATE coords

# Clean names
rownames(phate_coords) <- trimws(rownames(phate_coords))

colnames(seurat_obj) <- trimws(colnames(seurat_obj))

# Find common cells
common_cells <- intersect(colnames(seurat_obj), rownames(phate_coords))

# Subset and reorder Seurat object and PHATE coords
seurat_obj <- subset(seurat_obj, cells = common_cells)
phate_coords <- phate_coords[common_cells, , drop = FALSE]


seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Now make sure they're aligned
stopifnot(identical(rownames(phate_coords), colnames(seurat_obj)))  # This MUST pass

# Add PHATE
seurat_obj[["phate"]] <- CreateDimReducObject(
  embeddings = as.matrix(phate_coords),
  key = "PHATE_",
  assay = DefaultAssay(seurat_obj)
)

seurat_obj$leiden <- factor(leiden_clusters[colnames(seurat_obj), 1])


####



pdf(paste0("figures/pseudotime_analysis/", output_prefix, "/phate_umap_plots.pdf"),
    width = 18, height = 10)

# Individual UMAP plots
p1 <- DimPlot(seurat_obj, reduction = "phate", group.by = "Stage", label = FALSE) +
  ggtitle("Stage") +
  scale_color_manual(values = stage_colors)

p2 <- DimPlot(seurat_obj, reduction = "phate", group.by = "broad_cell_types", label = FALSE) +
  ggtitle("Broad Cell Types") +
  scale_color_manual(values = broad_cell_types_colors)

p3 <- DimPlot(seurat_obj, reduction = "phate", group.by = "infection_status", label = FALSE) +
  ggtitle("Infection Status") +
  scale_color_manual(values = infection_status_colors)

p4 <- DimPlot(seurat_obj, reduction = "phate", group.by = "niche_with_tumor_proximity", label = FALSE) +
  ggtitle("Niches")

# Arrange plots in 2x2 layout
combined_plot <- (p1 | p2) / (p3 | p4)

# Display
print(combined_plot)

dev.off()



###

DimPlot(seurat_obj, reduction = "phate", group.by = "leiden", label = FALSE) +
  ggtitle("Leiden")



# Run slingshot
start_cluster <- c('0') # Adjust this based on your data

sce <- as.SingleCellExperiment(seurat_obj)

# Add reduced dimensions from UMAP
reducedDims(sce)$PHATE <- seurat_obj@reductions$phate@cell.embeddings


message("Running Slingshot for cluster ", start_cluster, "...")
sce <- slingshot(sce, 
                 clusterLabels = 'leiden',
                 reducedDim = 'PHATE',
                 start.clus = start_cluster,
                 end.clus = c('1'))





#####################



# Plot Trajectories with Cell Types
curve_list_1 <- slingCurves(sce)
curve_df_1 <- do.call(rbind, lapply(seq_along(curve_list_1), function(i) {
  crv <- curve_list_1[[i]]
  if (is.null(crv)) return(NULL)
  data.frame(
    PHATE_1 = crv$s[,1],
    PHATE_2 = crv$s[,2],
    curve_id = paste0("lineage_", i),
    ord = seq_len(nrow(crv$s))
  )
}))

p1 <- DimPlot(seurat_obj, reduction = "phate", group.by = "broad_cell_types", label = FALSE) +
  ggtitle(paste("Broad Cell Types + Slingshot (Cluster", start_cluster, ")")) +
  scale_color_manual(values = broad_cell_types_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1)

ggsave(
  filename = paste0("figures/pseudotime_analysis/", output_prefix, "_phate_with_trajectories.pdf"),
  plot = p1, width = 8, height = 6
)

p1



# Set the output PDF file and size (adjust width/height as needed)
pdf(paste0("figures/pseudotime_analysis/", output_prefix, "/phate_plots_with_trajectories.pdf"),
    width = 12, height = 10)

# Your plotting code
p1 <- DimPlot(seurat_obj, reduction = "phate", group.by = "broad_cell_types", label = FALSE) +
  ggtitle("Broad Cell Types") +
  scale_color_manual(values = broad_cell_types_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1)

p2 <- DimPlot(seurat_obj, reduction = "phate", group.by = "Stage", label = FALSE) +
  ggtitle("Stage") +
  scale_color_manual(values = stage_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1)

p3 <- DimPlot(seurat_obj, reduction = "phate", group.by = "infection_status", label = FALSE) +
  ggtitle("Infection Status") +
  scale_color_manual(values = infection_status_colors) + 
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1)

p4 <- DimPlot(seurat_obj, reduction = "phate", group.by = "niche_with_tumor_proximity", label = FALSE) +
  ggtitle("Niches") +
  scale_color_manual(values = niche_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1)


#p4 <- DimPlot(seurat_obj, reduction = "phate", group.by = "leiden", label = TRUE) +
#  ggtitle("Leiden Clusters")

combined_plot <- (p1 | p2) / (p3 | p4)

print(combined_plot)

# Close the PDF device
dev.off()


# Set the output PDF file and size (adjust width/height as needed)
pdf(paste0("figures/pseudotime_analysis/", output_prefix, "/phate_plots_with_trajectories.pdf"),
    width = 12, height = 10)

# Plotting code with legends removed
p1 <- DimPlot(seurat_obj, reduction = "phate", group.by = "broad_cell_types", label = FALSE) +
  ggtitle("Broad Cell Types") +
  scale_color_manual(values = broad_cell_types_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1) +
  theme(legend.position = "none")

p2 <- DimPlot(seurat_obj, reduction = "phate", group.by = "Stage", label = FALSE) +
  ggtitle("Stage") +
  scale_color_manual(values = stage_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1) +
  theme(legend.position = "none")

p3 <- DimPlot(seurat_obj, reduction = "phate", group.by = "infection_status", label = FALSE) +
  ggtitle("Infection Status") +
  scale_color_manual(values = infection_status_colors) + 
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1) +
  theme(legend.position = "none")

p4 <- DimPlot(seurat_obj, reduction = "phate", group.by = "niche_with_tumor_proximity", label = FALSE) +
  ggtitle("Niches") +
  scale_color_manual(values = niche_colors) +
  geom_path(data = curve_df_1, aes(x = PHATE_1, y = PHATE_2, group = curve_id), color = "black", size = 1) +
  theme(legend.position = "none")

# Combine and print
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

# Close the PDF device
dev.off()


sds <- SlingshotDataSet(sce)
pt_matrix <- slingPseudotime(sce)
phate_coords <- reducedDims(sce)$PHATE

saveRDS(sce, file = paste0("results/pseudotime_analysis/", output_prefix, "/", output_prefix, "_sce.rds"))
saveRDS(sds, file = paste0("results/pseudotime_analysis/", output_prefix, "/", output_prefix, "_slingshot_dataset.rds"))
saveRDS(pt_matrix, file = paste0("results/pseudotime_analysis/", output_prefix, "/", output_prefix, "_pseudotime_matrix.rds"))




# Drop problematic entries from colData (like slingshot output)
to_drop <- c("slingshot", "tradeSeq", "crv")
existing <- colnames(colData(sce))
drop_cols <- intersect(to_drop, existing)
colData(sce) <- colData(sce)[, !colnames(colData(sce)) %in% drop_cols, drop = FALSE]

# Also clear metadata to be safe
metadata(sce) <- list()

cd <- colData(sce)
cd <- cd[, !duplicated(colnames(cd)), drop = FALSE]
colData(sce) <- cd

# Now write to h5ad
writeH5AD(sce, file = paste0("results/pseudotime_analysis/", output_prefix, "/", output_prefix, "_sce.h5ad"), compression = "gzip")

