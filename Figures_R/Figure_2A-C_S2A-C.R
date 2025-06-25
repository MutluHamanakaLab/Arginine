# Load Libraries
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)

#Load Data
level3.seurat <- readRDS("./PATH/TO/DIRECTORY/fibroblast-lean-seurat-20240308.RDS") # SCT transformed preprocessed fibroblast population

# Assigning new cluster identities
level3.seurat_L <- FindClusters(level3.seurat, resolution = 0.3)

# Assigning new cluster identities
new.cluster.ids <- c("Myofibroblasts", "PLIN2 High Fibroblasts", "Fibroblasts", "Fibroblasts", "Myofibroblasts", "HAS1 High Fibroblasts", "PLIN2 High Fibroblasts",
                     "Myofibroblasts", "Myofibroblasts")
names(new.cluster.ids) <- levels(level3.seurat_L)
level3.seurat_L <- RenameIdents(level3.seurat_L, new.cluster.ids)

# Define the new order of clusters
new_order <- c("HAS1 High Fibroblasts", "PLIN2 High Fibroblasts", "Myofibroblasts", "Fibroblasts")

# Reorder the levels of the Idents
level3.seurat_L@active.ident <- factor(level3.seurat_L@active.ident, levels = new_order)

# Simplify the label
DefaultAssay(level3.seurat_L) <- "RNA" # rds file is at SCT assay
new.cluster.ids <- c(
  "Fibroblasts" = "Fibroblasts",
  "Myofibroblasts" = "Myofibroblasts",
  "PLIN2 High Fibroblasts" = "PLIN2+",
  "HAS1 High Fibroblasts" = "HAS1 High"
)

level3.seurat_L <- RenameIdents(level3.seurat_L, new.cluster.ids)

# Normalize RNA assay
level3.seurat_L <- NormalizeData(level3.seurat_L)
level3.seurat_L <- ScaleData(level3.seurat_L)
level3.seurat_L <- FindVariableFeatures(level3.seurat_L)

# Figure 2A
scDimplot <- DimPlot(level3.seurat_L)+ xlim(-8,8) + ylim(-5,5)
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_Legend.png",scDimplot, width = 25, height = 6, units = "cm")

scDimplot <- DimPlot(level3.seurat_L) + NoLegend() + xlim(-8,8) + ylim(-5,5)
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_noLegend.png",scDimplot, width = 20, height = 6, units = "cm")

# Figure S2A
scDimplot <- DimPlot(level3.seurat_L, split.by = "condition")+ xlim(-8,8) + ylim(-5,5)
ggsave("./PATH/TO/DIRECTORY/scDimplot_level3_split_re.png",scDimplot, width = 25, height = 6, units = "cm")

# Organize Idents for the plot later (Default)
Idents(level3.seurat_L) <- factor(Idents(level3.seurat_L), levels = c("HAS1 High", "PLIN2+", "Myofibroblasts", "Fibroblasts"))

# ---
# Dot Plot
# Define Gene List
arg_metabolism <- c("ARG2", "ASS1", "ASL", "NOS3")

# Figure 2B
scDotplot <- DotPlot(level3.seurat_L, features = arg_metabolism) + RotatedAxis()
W <- length(arg_metabolism) * 1.3
H <- length(unique(level3.seurat_L@active.ident)) *1.1
ggsave("./PATH/TO/DIRECTORY/scDotplot_level3_arg_HIGH_RNA_Normalized.png",scDotplot, height = H, width = W)
H <- length(unique(level3.seurat_L@active.ident)) *0.5
ggsave("./PATH/TO/DIRECTORY/scDotplot_level3_arg_LOW_RNA_Normalized.png",scDotplot, height = H, width = W)

# ---
# Feature Plot
# Define Plot Dimension
plot_width <- 10
plot_height <- 6

# Figure 2C
for (gene in arg_metabolism) {
  plot <- FeaturePlot(level3.seurat_L, features = gene, reduction = "umap", label = F, order = T)
  ggsave(filename = paste0("./PATH/TO/DIRECTORY/scFeaturePlots_level3_", gene, "_RNA_Normalized.png"), plot = plot, width = plot_width, height = plot_height, units = "cm")
}

# Figure S2B
for (gene in arg_metabolism) {
  plot <- FeaturePlot_scCustom(level3.seurat_L, split.by = "condition", features = gene, order = T, colors_use = c("lightgrey", "blue"))
  ggsave(paste0("./PATH/TO/DIRECTORY/scFeaturePlot_level3_split_", gene,"_RNA_Norm.png"),
         width = 20,
         height = 6, units = "cm")
}

# Define Genes List
S1Gene <- c("PI16", "CTHRC1", "COL1A1")

# Figure S2C Loop
for (gene in S1Gene) {
  plots <- FeaturePlot(level3.seurat_L, features = gene, order = T)
  ggsave(paste0("./output/2025-03-26/scFeaturePlot_level3_FigS1_",gene,"_RNA_Normalized.png"),plots, width = 10, height = 6, units = "cm")
}

# Figure S2C cut off lower feature expressed cells
plots <- FeaturePlot(level3.seurat_L, features = "COL1A1", order = T, min.cutoff = "q85")
ggsave("./PATH/TO/DIRECTORY/scFeaturePlot_level3_FigS1_COL1A1_RNA_Normalized.png",plots, width = 10, height = 6, units = "cm")
plots <- FeaturePlot(level3.seurat_L, features = "CTHRC1", order = T, min.cutoff = "q70")
ggsave("./PATH/TO/DIRECTORY/scFeaturePlot_level3_FigS1_CTHRC1_RNA_Normalized.png",plots, width = 10, height = 6, units = "cm")
