library(reticulate)
use_condaenv('py311')
py_config()

#install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
#install.packages("BiocManager")
#BiocManager::install(c("GSVA", "GSEABase"))

if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  BiocManager::install("zellkonverter")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterExperiment")


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
  'Dendritic cells' = '#ff9300'
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

dir.create("results/pseudotime_analysis/with_macrophages", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/pseudotime_analysis/with_macrophages", recursive = TRUE, showWarnings = FALSE)





# Load the extracted data
expression_matrix <- as.matrix(read.csv("data/expression_matrix_subset_with_macrophages.csv", row.names=1))
cell_metadata <- read.csv("data/cell_metadata_subset_with_macrophages.csv", row.names=1)
gene_metadata <- read.csv("data/gene_metadata_subset_with_macrophages.csv", row.names=1)
gene_names <- read.csv("data/gene_names_subset_with_macrophages.csv", header=FALSE, stringsAsFactors=FALSE)[[1]]
phate_coords <- read.csv("data/phate_embedding_with_macrophages.csv", row.names=1)
leiden_clusters <- read.csv("data/leiden_cluster_ids_with_macrophages.csv", row.names = 1)

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

# Perform Leiden clustering with resolution 0.5
seurat_obj <- FindNeighbors(seurat_obj, reduction = "phate", dims = 1:2)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, algorithm = 4)  # algorithm 4 is for Leiden

# Assign the new clusters
seurat_obj$leiden <- Idents(seurat_obj)

# Convert to SingleCellExperiment for Slingshot
sce <- as.SingleCellExperiment(seurat_obj)

# Add PHATE coordinates as reduced dimensions
reducedDims(sce)$PHATE <- seurat_obj@reductions$phate@cell.embeddings

# Run Slingshot with PHATE coordinates
# Using cluster 1 as starting point (macrophages) and clusters 2 and 0 as endpoints
sce <- slingshot(sce, 
                 clusterLabels = 'leiden',
                 reducedDim = 'PHATE',
                 start.clus = '1',
                 end.clus = c('2', '0'))

# Get the curves and pseudotime
curves <- slingCurves(sce)
pseudotime <- slingPseudotime(sce, na = FALSE)
cellWeights <- slingCurveWeights(sce)

# Create a data frame for plotting the curves
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  crv <- curves[[i]]
  if (is.null(crv)) return(NULL)
  data.frame(
    PHATE_1 = crv$s[,1],
    PHATE_2 = crv$s[,2],
    curve_id = paste0("lineage_", i),
    ord = seq_len(nrow(crv$s))
  )
}))

# Plot trajectories with cell types
p1 <- DimPlot(seurat_obj, reduction = "phate", group.by = "broad_cell_types", label = FALSE) +
  ggtitle("Broad Cell Types + Slingshot Trajectories") +
  scale_color_manual(values = broad_cell_types_colors) +
  geom_path(data = curve_df, aes(x = PHATE_1, y = PHATE_2, group = curve_id), 
            color = "black", size = 1)

# Save the plot
ggsave(
  filename = paste0("figures/pseudotime_analysis/with_macrophages/", output_prefix, "_phate_with_trajectories.pdf"),
  plot = p1, 
  width = 10, 
  height = 8
)

# Plot individual lineage pseudotime plots
n_trajectories <- ncol(pseudotime)
n_cols <- ceiling(sqrt(n_trajectories))
n_rows <- ceiling(n_trajectories / n_cols)

# Calculate dimensions for PNG (in pixels)
png_width <- 800 * n_cols
png_height <- 800 * n_rows

# Save as PNG
png(paste0("figures/pseudotime_analysis/with_macrophages/", output_prefix, "_pseudotime_by_lineage.png"), 
    width = png_width, height = png_height, res = 100)
par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 2, 1))  # Reduced margins

# Color scale
col_vector <- colorRampPalette(c("red", "blue"))(100)

# Plot each lineage separately
for (i in seq_len(n_trajectories)) {
  pt <- pseudotime[, i]
  colors <- rep("lightgray", length(pt))
  colors[!is.na(pt)] <- col_vector[cut(pt[!is.na(pt)], breaks = 100)]
  
  plot(
    seurat_obj@reductions$phate@cell.embeddings,
    col = colors,
    pch = 16,
    cex = 0.3,
    asp = 1,
    main = paste("Pseudotime - Lineage", i),
    axes = FALSE  # Remove axes to save space
  )
  
  # Extract only the i-th curve
  curve <- curves[[i]]
  if (!is.null(curve)) {
    lines(curve$s, lwd = 2, col = "black")
  }
}
dev.off()

# Add all pseudotime values to Seurat object metadata
for(i in 1:n_trajectories) {
  seurat_obj[[paste0("pseudotime_traj_", i)]] <- pseudotime[,i]
}

# Save the updated Seurat object
saveRDS(seurat_obj, file = paste0("results/pseudotime_analysis/with_macrophages/", output_prefix, "_seurat_with_pseudotime.rds"))

# Prepare for differential expression analysis
# Select variable genes for analysis
sel_cells <- split(colnames(seurat_obj), seurat_obj$leiden)
sel_cells <- unlist(lapply(sel_cells, function(x) {
  set.seed(1)
  return(sample(x, min(20, length(x))))
}))

# Calculate gene variance
gv <- as.data.frame(na.omit(scran::modelGeneVar(seurat_obj@assays$RNA@data[, sel_cells])))
gv <- gv[order(gv$bio, decreasing = T), ]
sel_genes <- sort(rownames(gv)[1:500])

# Fit GAM model for trajectory analysis
BiocParallel::register(BiocParallel::SnowParam())
sceGAM <- fitGAM(
  counts = drop0(seurat_obj@assays$RNA@data[sel_genes, sel_cells]),
  pseudotime = pseudotime[sel_cells, ],
  cellWeights = cellWeights[sel_cells, ],
  nknots = 5,
  verbose = TRUE,
  parallel = TRUE,
  sce = TRUE
)

# Find genes that change along pseudotime
res_time <- na.omit(associationTest(sceGAM, contrastType = "consecutive"))
res_time <- res_time[res_time$pvalue < 1e-3, ]
res_time <- res_time[res_time$waldStat > mean(res_time$waldStat), ]
res_time <- res_time[order(res_time$waldStat, decreasing = T), ]

# Find genes that differ between lineages
res_diff <- na.omit(diffEndTest(sceGAM))
res_diff <- res_diff[res_diff$pvalue < 1e-3, ]
res_diff <- res_diff[res_diff$waldStat > mean(res_diff$waldStat), ]
res_diff <- res_diff[order(res_diff$waldStat, decreasing = T), ]

# Save results
write.csv(res_time, "results/pseudotime_analysis/with_macrophages/genes_changing_with_time.csv")
write.csv(res_diff, "results/pseudotime_analysis/with_macrophages/genes_different_between_lineages.csv")

# Plot top changing genes
top_genes <- rownames(res_time)[1:16]
p3 <- FeaturePlot(seurat_obj, 
                 features = top_genes,
                 reduction = "phate",
                 ncol = 4,
                 order = TRUE) +
  plot_annotation(title = "Top Genes Changing with Pseudotime")

ggsave(
  filename = "figures/pseudotime_analysis/with_macrophages/top_changing_genes.pdf",
  plot = p3,
  width = 16,
  height = 16
)

# Save all results
saveRDS(sce, file = "results/pseudotime_analysis/with_macrophages/slingshot_results.rds")
saveRDS(sceGAM, file = "results/pseudotime_analysis/with_macrophages/gam_results.rds")
saveRDS(pseudotime, file = "results/pseudotime_analysis/with_macrophages/pseudotime_matrix.rds")

