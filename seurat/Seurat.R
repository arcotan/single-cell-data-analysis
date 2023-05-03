library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)

IN_DATA_DIR = "./dataset/tabulamuris/droplet/Heart_and_Aorta-10X_P7_4"
OUT_DATA_DIR = "./filtered_dataset/tabulamuris/"
IN_LABEL_DIR = "./dataset/tabulamuris"
OUT_LABEL_DIR = "./filtered_dataset/tabulamuris/"
CHANNEL = "10X_P7_4"

# loads gene expression matrix associated to a channel and gets the cluster label for each cell,
# if metadata for a cell is not found, its cluster id will be NA
load_data <- function(data_dir, label_dir, channel) {
  # Load gene expression matrix
  data = Read10X(data.dir = IN_DATA_DIR, strip.suffix = TRUE)
  
  # Load cell metadata
  metadata = read.csv(file = paste(IN_LABEL_DIR, "annotations_droplet.csv", sep = "/"))
  
  # Filter data to use only data for current dataset
  metadata = metadata[metadata$channel == channel,]
  metadata$cell = substr(metadata$cell, 10, 25)
  
  # Get cluster labels
  cells = colnames(data)
  true_labels = left_join(data.frame("cell"=cells), metadata)[c("cell", "cluster.ids")]
  
  # Print number of cells not mapped to a cluster id
  print(paste("Cells mapped to a cluster: ", sum(!is.na(true_labels$cluster.ids)), "/", nrow(true_labels), sep = ""))
  
  return (list("data" = data, "labels" = true_labels))
} 

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
write10xCounts(paste(OUT_DATA_DIR, "data_", CHANNEL, sep=""), data_to_write)
filtered_labels = left_join(data.frame("cell"=colnames(data_to_write)), experiment_data$labels)
write.csv(filtered_labels, paste(OUT_LABEL_DIR, "labels_", CHANNEL, ".csv", sep=""), row.names = FALSE)

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

label_df = merge(experiment_data$labels, data.frame(Idents(pbmc)), by.x = "cell", by.y = 0) %>% 
  rename(
    true_id = cluster.ids,
    computed_id = Idents.pbmc.
  )

label_df$true_id = label_df$true_id + 1
label_df$computed_id = as.numeric(label_df$computed_id)

label_df_nna = label_df[!is.na(label_df$true_id), ]
sum(label_df_nna$true_id == label_df_nna$computed_id)

# In input: a dataframe with two non NA cluster id columns, named 'true_id' and 'computed_id'
# labels of these columns will be renamed in order to get the best match between clusters.
# The confusion matrix of the clustering is also added to the output:
# rows of the matrix are related to predictions, columns are related to the ground truth
align_clusters = function(label_dataframe) {
  confusion_matrix = table(label_dataframe$computed_id, label_dataframe$true_id) # row predictions, col ground truth
  permutation_computed = c(1:nrow(confusion_matrix))
  permutation_true = c(1:ncol(confusion_matrix))
  it_num = min(nrow(confusion_matrix), ncol(confusion_matrix))
  
  for (i in 1:it_num) {
    #get arg max of sub matrix M[i:end, i:end] as a single value 'max_pos'
    max_pos = unname(which.max(confusion_matrix[i:nrow(confusion_matrix), i:ncol(confusion_matrix)]))
    #get the row and column related to 'max_pos'
    row_max = (max_pos %% (nrow(confusion_matrix) - i + 1))
    if (row_max == 0) {#last row treated as a special case since we are not counting from zero
      row_max = nrow(confusion_matrix)
      col_max = (max_pos %/% (nrow(confusion_matrix) - i + 1)) + i - 1
    }
    else {
      row_max = row_max + i - 1
      col_max = (max_pos %/% (nrow(confusion_matrix) - i + 1)) + i
    }
    
    #apply pivoting to matrix and permutations
    confusion_matrix[, c(i, col_max)] = confusion_matrix[, c(col_max, i)]
    confusion_matrix[c(i, row_max), ] = confusion_matrix[c(row_max, i), ]
    permutation_computed[c(i, row_max)] = permutation_computed[c(row_max, i)]
    permutation_true[c(i, col_max)] = permutation_true[c(col_max, i)]
  }
  
  # apply permutations
  label_dataframe$computed_id = permutation_computed[label_dataframe$computed_id]
  label_dataframe$true_id = permutation_true[label_dataframe$true_id]
  
  # return confusion matrix
  return (list("confusion_matrix" = confusion_matrix, "label_dataframe" = label_dataframe))
}

# rename cluster labels to get best clustering results
clustering_info = align_clusters(label_df_nna)
confusion_matrix = clustering_info$confusion_matrix
label_df_nna = clustering_info$label_dataframe
# display confusion matrix
confusion_matrix

# modifies seurat identities
# returns clustering plot with pca
seurat_clustering_plot = function(seurat_obj, cell_col, label_col) {
  seurat_obj <- SetIdent(seurat_obj, cells = cell_col, label_col)
  return (DimPlot(seurat_obj, reduction = "pca"))
}

# Compare clustering (left) with ground truth
seurat_clustering_plot(pbmc, label_df_nna$cell, label_df_nna$computed_id) + seurat_clustering_plot(pbmc, label_df_nna$cell, label_df_nna$true_id)

# find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
