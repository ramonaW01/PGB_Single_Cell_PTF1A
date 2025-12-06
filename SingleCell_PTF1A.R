install.packages("R.utils") 

library(data.table)
library(Seurat)
library(ggplot2)
library(utils)

# 1. Dataset
file_path <- "Documents/PGB/singlecell/Projekt/GSE84133_RAW/GSM2230757_human1_umifm_counts.csv.gz"
full_data <- fread(file_path)

unique_barcodes <- make.unique(full_data$barcode)
cell_metadata <- data.frame(
  row.names = unique_barcodes,
  original_barcode = full_data$barcode,  
  assigned_cluster = full_data$assigned_cluster
)
gene_columns <- colnames(full_data)[4:length(colnames(full_data))]
counts_matrix_wide <- as.matrix(full_data[, ..gene_columns])
counts_matrix_t <- t(counts_matrix_wide)

colnames(counts_matrix_t) <- unique_barcodes 
rownames(counts_matrix_t) <- gene_columns

# Seurat object 
pbmc <- CreateSeuratObject(counts = counts_matrix_t)
pbmc <- AddMetaData(pbmc, metadata = cell_metadata)

print(pbmc)

# 2. QC & filter

# QC-Metriken 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC-Metriken before filter
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")

# Filter
pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 60000 & percent.mt < 10)

print("filtered objec:")
print(pbmc_filtered)

# QC-Metriken after Filterung 
VlnPlot(pbmc_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# 3. STANDARD-ANALYSE-WORKFLOW:
# Normalizing
pbmc_filtered <- NormalizeData(pbmc_filtered)

# find variable features 
pbmc_filtered <- FindVariableFeatures(pbmc_filtered, selection.method = "vst", nfeatures = 2000)

# scaling all genes
all.genes <- rownames(pbmc_filtered)
pbmc_filtered <- ScaleData(pbmc_filtered, features = all.genes)

# PCA
pbmc_filtered <- RunPCA(pbmc_filtered, features = VariableFeatures(object = pbmc_filtered))

# Dim reduction (UMAP & tSNE)
# Test with ElbowPlot(pbmc_filtered) to see what dim value is fitting
ElbowPlot(pbmc_filtered)
pbmc_filtered <- RunUMAP(pbmc_filtered, dims = 1:20)
pbmc_filtered <- RunTSNE(pbmc_filtered, dims = 1:20)

# Clustering
pbmc_filtered <- FindNeighbors(pbmc_filtered, dims = 1:20)
pbmc_filtered <- FindClusters(pbmc_filtered, resolution = 0.5)


# 4. Visualisation
# UMAP Cluster
DimPlot(pbmc_filtered, reduction = "umap", group.by = "assigned_cluster", label = TRUE)

# tSNE Cluster
DimPlot(pbmc_filtered, reduction = "tsne", group.by = "assigned_cluster", label = TRUE)

# TF-Signal
FeaturePlot(pbmc_filtered, features = "PTF1A")
FeaturePlot(pbmc_filtered, reduction = "tsne", features = "PTF1A")

# Markergene-Expression
# Acinar-Markergene - found in Chipseq analysis
FeaturePlot(pbmc_filtered, features = c("CPA1", "AMY2A"))


# Expression in VlnPlot
VlnPlot(pbmc_filtered, features = "PTF1A", group.by = "assigned_cluster")
VlnPlot(pbmc_filtered, features = c("CPA1", "AMY2A"), group.by = "assigned_cluster")

# calculating average Expression
AverageExpression(pbmc_filtered, features = "PTF1A")$RNA


#visualisation of markergene AMY2A correlating with PTF1A

plot_data <- FetchData(pbmc_filtered, vars = c("PTF1A", "AMY2A"))

plot_data$both_expressed <- (plot_data$PTF1A > 0 & plot_data$AMY2A > 0)

ggplot(plot_data, aes(x = PTF1A, y = AMY2A)) +
  geom_point(aes(color = both_expressed), alpha = 0.5) +
  scale_color_manual(
    values = c("grey70", "red"),
    labels = c("no", "yes"),
    name = "Expression of both genes"
  ) +
  ggtitle("Expression PTF1A vs AMY2A") +
  theme_minimal()



#visualisation of markergene CPA1 correlating with PTF1A
plot_data2 <- FetchData(pbmc_filtered, vars = c("PTF1A", "CPA1"))

plot_data2$both_expressed <- (plot_data2$PTF1A > 0 & plot_data2$CPA1 > 0)

ggplot(plot_data2, aes(x = PTF1A, y = CPA1)) +
  geom_point(aes(color = both_expressed), alpha = 0.5) +
  scale_color_manual(
    values = c("grey70", "red"),
    labels = c("no", "yes"),
    name = "Expression of both genes"
  ) +
  ggtitle("Expression PTF1A vs CPA1") +
  theme_minimal()

