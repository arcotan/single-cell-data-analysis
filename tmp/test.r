library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

data <- readRDS(
    "dataset/DuoClustering2018/sce_full/sce_full_SimKumar4hard.rds"
)


counts <- assays(data)$counts
srat <- CreateSeuratObject(counts = counts, project = "kumar4hard", min.cells = 3, min.features = 200)
srat
adj.matrix <- NULL
str(srat)
meta <- srat@meta.data
dim(meta)
head(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)
# srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
VlnPlot(srat,
    features = c("nFeature_RNA", "nCount_RNA"),
    ncol = 2,
    pt.size = 0.1
) &  theme(plot.title = element_text(size=10))
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)


srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(srat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(srat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)


srat <- RunPCA(srat, features = VariableFeatures(object = srat))

# Examine and visualize PCA results a few different ways
print(srat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(srat, dims = 1:2, reduction = "pca")

DimPlot(srat, reduction = "pca")

DimHeatmap(srat, dims = 1, cells = 500, balanced = TRUE)
