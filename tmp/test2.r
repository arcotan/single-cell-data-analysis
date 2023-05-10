library(ExperimentHub)
library(Seurat)
library(ggplot2)

eh <- ExperimentHub()
q <- query(eh, "DuoClustering2018")
kumar4hard = q[["EH1540"]]
#kumar8hard = q[["EH1544"]]
ge = counts(kumar4hard)

srat <- CreateSeuratObject(counts = ge, project = "kumar4hard", min.cells = 3, min.features = 200)

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

srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.5)

# display confusion matrix
clusters_seurat = as.numeric(Idents(srat))
true_clusters = as.numeric(substr(kumar4hard$Group,6,6))
confusion_matrix = matrix(0,4,4) # row predictions, col ground truth
for (i in 1:length(true_clusters)) {
  confusion_matrix[clusters_seurat[i],true_clusters[i]] = confusion_matrix[clusters_seurat[i],true_clusters[i]] + 1
}
confusion_matrix