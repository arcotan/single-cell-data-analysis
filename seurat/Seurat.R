library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)

source("utils.R")

IN_DATA_DIR = "./dataset/tabulamuris/droplet/Heart_and_Aorta-10X_P7_4"
OUT_DATA_DIR = "./filtered_dataset/tabulamuris/"
IN_LABEL_DIR = "./dataset/tabulamuris"
OUT_LABEL_DIR = "./filtered_dataset/tabulamuris/"
CHANNEL = "10X_P7_4"
OUT_RES_DIR = "./results/seurat"
TOP_MARKER_NUM = 20
RES_FILE_TAG = CHANNEL

# Loading data
experiment_data = load_data(IN_DATA_DIR, IN_LABEL_DIR, CHANNEL)

# Load the PBMC dataset
pbmc.data <- experiment_data$data

# plot distribution of amount of cells in which each gene is expressed
ggplot(data.frame("sum" = rowSums(pbmc.data > 0)), aes(x=sum)) + geom_histogram()
# plot distribution of sum of counts for each cell
ggplot(data.frame("sum" = colSums(pbmc.data > 0)), aes(x=sum)) + geom_histogram()

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
Dir10X = paste(OUT_DATA_DIR, "data_", CHANNEL, sep="")
if (!dir.exists(Dir10X)) {
  write10xCounts(Dir10X, data_to_write)
}

filtered_labels = left_join(data.frame("cell"=colnames(data_to_write)), experiment_data$labels)
write.csv(filtered_labels, paste(OUT_LABEL_DIR, "labels_", CHANNEL, ".csv", sep=""), row.names = FALSE)

write.csv(experiment_data$mapping, paste(OUT_LABEL_DIR, "mapping_", CHANNEL, ".csv", sep=""), row.names = FALSE)

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
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)

plot_de(GetAssayData(object = pbmc, assay = "RNA", slot = "data"), pbmc.markers, "gene", "cluster", label_df, "cell", "computed_id", OUT_RES_DIR, RES_FILE_TAG)

write_clustering(OUT_RES_DIR, RES_FILE_TAG, label_df, "cell", "computed_id", "true_id", dist(Embeddings(pbmc[['pca']])[,1:50]))

write_markers(OUT_RES_DIR, RES_FILE_TAG, pbmc.markers, "gene", "cluster", "avg_log2FC", TRUE, TOP_MARKER_NUM)

