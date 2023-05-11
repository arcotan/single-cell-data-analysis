
library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)

source("utils.R")


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) # unit variance?

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(object = pbmc, 
          ndims = 50)

# Get clustering data
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.3)

label_df = merge(experiment_data$labels, data.frame(Idents(pbmc)), by.x = "cell", by.y = 0)
names(label_df)[2:3] <- c("true_id", "computed_id") #TODO non usare 2 e 3

label_df$computed_id = as.numeric(label_df$computed_id)
pbmc <- SetIdent(pbmc, cells = label_df$cell, label_df$computed_id)

# rename cluster labels to get best clustering results
clustering_info = align_clusters(label_df, "true_id", "computed_id")
confusion_matrix = clustering_info$confusion_matrix
label_df = clustering_info$label_dataframe
# display confusion matrix
confusion_matrix

# Compare clustering (left) with ground truth
seurat_clustering_plot(pbmc, label_df$cell, label_df$computed_id) + seurat_clustering_plot(pbmc, label_df$cell, label_df$true_id)

# find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.1, logfc.threshold = 0.15) # suggested by Silvia

plot_de(GetAssayData(object = pbmc, assay = "RNA", slot = "data"), pbmc.markers, "gene", "cluster", label_df, "cell", "computed_id", OUT_RES_DIR, RES_FILE_TAG)

write_clustering(OUT_RES_DIR, RES_FILE_TAG, label_df, "cell", "computed_id", "true_id", dist(Embeddings(pbmc[['pca']])[,1:50]))

write_markers(OUT_RES_DIR, RES_FILE_TAG, pbmc.markers, "gene", "cluster", "avg_log2FC", TRUE, TOP_MARKER_NUM)































