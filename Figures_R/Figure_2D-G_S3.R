# Load libraries
library(Seurat)
library(ggplot2)

# Load Data
level4.seurat_3 <- readRDS(file = "./PATH/TO/DIRECTORY/level4.seurat_3.rds") # SCT transformed preprocessed fibroblast (MyoFibroblast, Fibroblast except PI16 high)

# Set default assay as "RNA" rds file default assay is SCT
DefaultAssay(level4.seurat_3) <- "RNA"

# Normalize RNA assay
level4.seurat_3 <- NormalizeData(level4.seurat_3)
level4.seurat_3 <- ScaleData(level4.seurat_3)
level4.seurat_3 <- FindVariableFeatures(level4.seurat_3)

# Organize Idents for plot later (Default)
Idents(level4.seurat_3) <- factor(Idents(level4.seurat_3), levels = c("Fibrotic", "Inflammatory", "Alveolar"))

# ---
# Dim Plot
# Figure 2D
dimplot <- DimPlot(level4.seurat_3, cols = c("glasbey"),pt.size = 4) + NoLegend()
ggsave("./PATH/TO/DIRECTORY/scDimplot_level4_new_pt4_noLegend.png", dimplot, width = 3.5, height = 3.5)
dimplot <- DimPlot(level4.seurat_3, cols = c("glasbey"),pt.size = 4)
ggsave("./PATH/TO/DIRECTORY/scDimplot_level4_new_pt4_Legend.png", dimplot, width = 5, height = 3.5)

# ---
# Dot Plot
# Define a genes list
arg_metabolism <- c("ARG2", "ASS1", "ASL", "NOS3")

# Figure 2E
plot <- DotPlot(level4.seurat_3, features = arg_metabolism) + RotatedAxis()
W <- length(arg_metabolism) *1.2
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave(filename = "./PATH/TO/DIRECTORY/scDotPlots_level4_arg_metabolism_RNA_Normalized_HIGH.png", plot = plot, width = W, height = H)
H <- length(unique(level4.seurat_3@active.ident)) * 0.7
ggsave(filename = "./PATH/TO/DIRECTORY/scDotPlots_level4_arg_metabolism_RNA_Normalized_LOW.png", plot = plot, width = W, height = H)

# ---
# Feature Plot
# Define plot dimension
plot_width <- 3.5
plot_height <- 3.5

# Figure 2F
for (gene in arg_metabolism) {
  plot <- FeaturePlot(level4.seurat_3, features = gene, reduction = "umap", label = F, order = T)
  ggsave(filename = paste0("./PATH/TO/DIRECTORY/FeaturePlots_level4_", gene, "_RNA_Normalized.png"), plot = plot, width = plot_width, height = plot_height)
}

# Define the gene list
arg_rebut_genes_comb <- c("ALDH18A1", "PYCR1", "OAT", "ODC", "SRM", "ALDH4A1", "SAT1", "SAT2", "GLUL")

# Figure 2G
dotplot_combined_idents <- DotPlot(level4.seurat_3, features = arg_rebut_genes_comb) + RotatedAxis()
W <- length(arg_rebut_genes_comb) * 0.5
H <- length(unique(level4.seurat_3@active.ident)) *1.4
ggsave("./PATH/TO/DIRECTORY/scDotplot_level4_ARG_UpDown_Combined_RNA_TALL.png", dotplot_combined_idents, width = W, height = H)
H <- length(unique(level4.seurat_3@active.ident)) *0.7
ggsave("./PATH/TO/DIRECTORY/scDotplot_level4_ARG_UpDown_Combined_RNA_LOW.png", dotplot_combined_idents, width = W, height = H)

# Figure S3 cut off lower feature expressed cells
plots <- FeaturePlot(level4.seurat_3, features = "COL1A1", order = T, min.cutoff = "q85")
ggsave("./PATH/TO/DIRECTORY/scFeaturePlot_level4_FigS1_COL1A1_RNA_Normalized.png",plots, width = 3.5, height = 3.5)
plots <- FeaturePlot(level4.seurat_3, features = "CTHRC1", order = T, min.cutoff = "q70")
ggsave("./PATH/TO/DIRECTORY/scFeaturePlot_level4_FigS1_CTHRC1_RNA_Normalized.png",plots, width = 3.5, height = 3.5)
