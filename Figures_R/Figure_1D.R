# Load Necessary Libraries
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # UCSC known gene annotation for Homo sapiens
library(readr)
library(tximport)
library(data.table)
library(tidyverse)
library(biomaRt)

# Initialize variables for human gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- na.omit(tx2gene) # Remove any rows with missing values
tx2gene  # View the resulting tx2gene DataFrame

# Load Meta Data
samples <- read.table("./PATH/TO/DIRECTORY/IPF_meta.txt", header = TRUE)

# Load Data files
files <- paste0(dir1, "/", samples$Run, "/abundance.tsv")

# Ensure all files exist
all(file.exists(files))

# Import Transcript Quantification Data with tximport
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
txi.tx <- tximport(files, type = "kallisto", txOut = TRUE)
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts) # Check equality of counts
attributes(txi) # View attributes

# Write the abundance data to a DataFrame
read_abundance <- data.frame(txi$abundance)
oldcolnames <- colnames(read_abundance)
newcolnames <- samples$Name
setnames(read_abundance, old = oldcolnames, new = newcolnames)
read_abundance <- rownames_to_column(read_abundance, var = "entrezgene_id")

# Map IDs using biomaRt
# For human
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
id_map <- getBM(attributes=c('ensembl_gene_id',
                          'external_gene_name','entrezgene_id'),
             filters = 'entrezgene_id',
             values = read_abundance$entrezgene_id,
             mart = ensembl)

# Merge mapped IDs with abundance data
read_abundance <- read_abundance %>%  mutate(entrezgene_id = as.numeric(entrezgene_id))
read_abundance <- left_join(id_map, read_abundance, by = "entrezgene_id")

# Rename external_gene_name to Symbol
read_abundance <- read_abundance %>% dplyr::rename("Symbol" = "external_gene_name")

# Remove duplicate symbols and set columns as numeric
rep_Symbol <- dim(read_abundance)[1] - length(unique(read_abundance$Symbol))
read_abundance <- read_abundance[!duplicated(read_abundance$ensembl_gene_id), ]
read_abundance <- read_abundance[!duplicated(read_abundance$Symbol), ]
read_abundance[newcolnames] <- sapply(read_abundance[newcolnames],as.numeric)#set data columns as numeric
sapply(read_abundance, class)

# Write "TableOfCounts" data frame
TableOfCounts <- read_abundance[c(-1,-3)]
colnames <- c(samples$Name)
setnames(TableOfCounts, old = oldcolnames, new = colnames)

write_csv(TableOfCounts, "./PATH/TO/DIRECTORY/TableOfCounts_IPF_Male.csv")

# Modify table of counts for further analysis (Remove Symbol Column)
TableOfCounts_noSymbol <- data.frame(TableOfCounts, row.names = 1)

# Check for any missing values
sum(is.na(TableOfCounts_noSymbol))

# Load libarary for pheatmap
library(pheatmap)

# Define data frame for heatmap
counts.in_Clones <- TableOfCounts_noSymbol
cols_clones <- colnames(counts.in_Clones)

# Define the order of column
order <- c(
  grep("^M24_UT", cols_clones),
  grep("^M53_UT", cols_clones),
  grep("^M71_UT", cols_clones),
  grep("^I24_UT", cols_clones),
  grep("^I35_UT", cols_clones),
  grep("^I66_UT", cols_clones),
  grep("^M24_T", cols_clones),
  grep("^M53_T", cols_clones),
  grep("^M71_T", cols_clones),
  grep("^I24_T", cols_clones),
  grep("^I35_T", cols_clones),
  grep("^I66_T", cols_clones)
)

# Reorder the columns
counts.in_Clones <- counts.in_Clones[, order]

# Define the metadata
counts.metadata_Clones <- data.frame(
  dataset= c(colnames(counts.in_Clones)),
  Sample= c("M24_UT","M24_UT","M24_UT",
            "M53_UT","M53_UT","M53_UT",
            "M71_UT","M71_UT","M71_UT",
            "I24_UT","I24_UT","I24_UT",
            "I35_UT","I35_UT","I35_UT",
            "I66_UT","I66_UT","I66_UT",
            "M24_TGFb", "M24_TGFb", "M24_TGFb",
            "M53_TGFb","M53_TGFb","M53_TGFb",
            "M71_TGFb","M71_TGFb","M71_TGFb",
            "I24_TGFb","I24_TGFb","I24_TGFb",
            "I35_TGFb","I35_TGFb","I35_TGFb",
            "I66_TGFb","I66_TGFb", "I66_TGFb"
            ),
  stringsAsFactors = FALSE
)

# Define Groups
group_Clones <- counts.metadata_Clones$Sample

# Create DGEList Object
y_Clones <- DGEList(counts=counts.in_Clones,
                                  genes=row.names.data.frame(counts.in_Clones),
                                  group=group_Clones)

# Filter Low-Expressed Genes
keep_Clones <- rowSums(edgeR::cpm(y_Clones)>1)>=1
table(keep_Clones)
y_Clones<- y_Clones[keep_Clones, , keep.lib.sizes=FALSE]

# Normalize the Data
y_Clones<- calcNormFactors(y_Clones, method = "TMM")
counts_Clones <- as.matrix(y_Clones$counts)
logCPM_Clones <- edgeR::cpm(counts_Clones, prior.count=1, log=TRUE)

# Calculate Z-Score
ZScore_Clones <- t(scale(t(logCPM_Clones)))
ZScore_Clones <- as.data.frame(ZScore_Clones)

# Define "KEGG Arginine and Proline Metabolism" genes lists
kegg_arg_G <- c("ACY1","AGMAT","ALDH18A1","ALDH1B1","ALDH2","ALDH3A2","ALDH4A1",
                "ALDH7A1","ALDH9A1","AMD1","AOC1","ARG1","ARG2","ASL","ASS1",
                "AZIN2","CKB","CKM","CKMT1A","CKMT1B","CKMT2","CPS1","DAO",
                "GAMT","GATM","GLS","GLS2","GLUD1","GLUD2","GLUL","GOT1","GOT2",
                "LAP3","MAOA","MAOB","NAGS","NOS1","NOS2","NOS3","OAT","ODC1",
                "OTC","P4HA1","P4HA2","P4HA3","PRODH","PRODH2","PYCR1","PYCR2",
                "PYCR3","SAT1","SAT2","SMS","SRM")

# Generate Z-Score Matrix for Genes Of Interests
sig.zscore_M_v_I_arg_kg <- ZScore_Clones[kegg_arg_G, ]
sig.zscore.mat_M_v_I_arg_kg <- as.matrix(sig.zscore_M_v_I_arg_kg)

# ---------
# This is done to remove genes like "CKMT2.1" which is duplicate of "CKMT2"
# Convert the matrix to a data frame
df <- as.data.frame(sig.zscore.mat_M_v_I_arg_kg)

# Remove duplicate rows 
df_unique <- df[!duplicated(df), ]

# Convert back to a matrix
sig.zscore.mat_M_v_I_arg_kg <- as.matrix(df_unique)

# ---------
# Ensures the heatmap only plots complete data
sig.zscore.mat_M_v_I_arg_kg <- sig.zscore.mat_M_v_I_arg_kg[complete.cases(sig.zscore.mat_M_v_I_arg_kg),]

# Prepare Heatmap Annotation
heat.annotation_M_v_I_arg_kg <- data.frame(counts.metadata_Clones[,2])
colnames(heat.annotation_M_v_I_arg_kg) <- "Sample"
row.names(heat.annotation_M_v_I_arg_kg) <- counts.metadata_Clones[,1]

# Set Column Names for Heatmap
d.colnames_M_v_I_arg_kg <- c(counts.metadata_Clones[,1])
colnames(sig.zscore.mat_M_v_I_arg_kg) <- d.colnames_M_v_I_arg_kg

# Define Annotation Colors
ann.colors_M_v_I_arg_kg <- list(Sample= c(`M24_UT` = "red",`M53_UT` = "orange", `M71_UT` = "yellow",
                                        `I24_UT` = "gray",`I35_UT` = "black",`I66_UT` = "lightyellow",
                                        `M24_TGFb` = "green",`M53_TGFb` = "blue",`M71_TGFb` = "violet",
                                        `I24_TGFb` = "darkgreen",`I35_TGFb` = "lightblue",`I66_TGFb` = "brown"))

# Italicize Gene Names
newnames_M_v_I_arg_kg <- lapply(
  rownames(sig.zscore.mat_M_v_I_arg_kg),
  function(x) bquote(italic(.(x))))

# Create Heatmap
# Figure 1D
M_v_I_arg_kg_heatmap <- pheatmap(sig.zscore.mat_M_v_I_arg_kg ,
                         annotation_col = heat.annotation_M_v_I_arg_kg, 
                         cluster_cols = FALSE,
                         main="",
                         annotation_colors = ann.colors_M_v_I_arg_kg,
                         show_colnames = FALSE,
                         labels_row = as.expression(newnames_M_v_I_arg_kg),
                         fontsize_row = 15
)

# Save the Plot
ggsave("./PATH/To/DIRECTORY/Male_vs_IPF_heatmap_arg_kg.png", plot = M_v_I_arg_kg_heatmap, width = 10, height = 15, limitsize = FALSE)
