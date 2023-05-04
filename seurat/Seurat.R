library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)
library(NMF)

source("utils.R")

IN_DATA_DIR = "./dataset/tabulamuris/droplet/Heart_and_Aorta-10X_P7_4"
OUT_DATA_DIR = "./filtered_dataset/tabulamuris/"
IN_LABEL_DIR = "./dataset/tabulamuris"
OUT_LABEL_DIR = "./filtered_dataset/tabulamuris/"
CHANNEL = "10X_P7_4"


# Loading data
experiment_data = load_data(IN_DATA_DIR, IN_LABEL_DIR, CHANNEL)

# Load the PBMC dataset
pbmc.data <- experiment_data$data

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Write pre-processed data
data_to_write = GetAssayData(object = pbmc, assay = "RNA", slot = "data")
write.csv(data_to_write, paste(OUT_DATA_DIR, "data_", CHANNEL, ".csv", sep=""))
write10xCounts(paste(OUT_DATA_DIR, "data_", CHANNEL, sep=""), data_to_write)
filtered_labels = left_join(data.frame("cell"=colnames(data_to_write)), experiment_data$labels)
write.csv(filtered_labels, paste(OUT_LABEL_DIR, "labels_", CHANNEL, ".csv", sep=""), row.names = FALSE)


# *******  END PREPROCESSING *******


pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Get clustering data
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.4)

label_df = merge(experiment_data$labels, data.frame(Idents(pbmc)), by.x = "cell", by.y = 0)
names(label_df)[2:3] <- c("true_id", "computed_id")

label_df$computed_id = as.numeric(label_df$computed_id)

# rename cluster labels to get best clustering results
clustering_info = align_clusters(label_df)
confusion_matrix = clustering_info$confusion_matrix
label_df = clustering_info$label_dataframe
# display confusion matrix
confusion_matrix

# modifies seurat identities
# returns clustering plot with pca
seurat_clustering_plot = function(seurat_obj, cell_col, label_col) {
  pi = order(label_col)
  cell_col = cell_col[pi]
  label_col = label_col[pi]
  seurat_obj <- SetIdent(seurat_obj, cells = cell_col, label_col)
  return (DimPlot(seurat_obj, reduction = "pca"))
}

# Compare clustering (left) with ground truth
seurat_clustering_plot(pbmc, label_df$cell, label_df$computed_id) + seurat_clustering_plot(pbmc, label_df$cell, label_df$true_id)

# find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

entropy(confusion_matrix)
purity(confusion_matrix)
mean(silhouette(label_df$computed_id, dist(Embeddings(pbmc[['pca']])[,1:50]))[,3])