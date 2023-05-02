library(dplyr)
library(Seurat)
library(patchwork)
library(fossil)

DATA_DIR = "./dataset/tabulamuris/droplet/Heart_and_Aorta-10X_P7_4"
LABEL_DIR = "./dataset/tabulamuris"
CHANNEL = "10X_P7_4"

# loads gene expression matrix associated to a channel and gets the cluster label for each cell,
# if metadata for a cell is not found, its cluster id will be NA
load_data <- function(data_dir, label_dir, channel) {
  # Load gene expression matrix
  data = Read10X(data.dir = DATA_DIR, strip.suffix = TRUE)
  
  # Load cell metadata
  metadata = read.csv(file = paste(LABEL_DIR, "annotations_droplet.csv", sep = "/"))
  
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
experiment_data = load_data(DATA_DIR, LABEL_DIR, CHANNEL)

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
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Get clustering data
pbmc <- FindNeighbors(pbmc, dims = 1:10, reduction="pca")
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

adj.rand.index(label_df_nna$computed_id, label_df_nna$true_id)

# rename cluster labels to get best clustering results
clustering_info = align_clusters(label_df_nna)
confusion_matrix = clustering_info$confusion_matrix
label_df_nna = clustering_info$label_dataframe
# display confusion matrix
confusion_matrix

# [TODO] la funzione ha un bug per qualche motivo
# modifies seurat identities
# returns clustering plot with pca
seurat_clustering_plot = function(seurat_obj, label_dataframe, label_col) {
  print(paste(seurat_obj, label_dataframe, label_col))
  SetIdent(seurat_obj, cells = label_dataframe$cell, label_dataframe$label_col)
  return (DimPlot(seurat_obj, reduction = "pca"))
}


pbmc <- RunTSNE(pbmc, features = VariableFeatures(object = pbmc))


cluster_plot = DimPlot(pbmc, reduction = "pca") # predicted
Idents(pbmc) <- label_df_nna$"true_id"

true_plot = DimPlot(pbmc, reduction = "pca")

cluster_plot + true_plot



# Compare clustering (left) with ground truth
cluster_plot + true_plot

pbmc@ident = label_df_nna$"true_id"

# pippo <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.1, test.use="wilcox", group.by=label_df_nna$"true_id")
pippo <- FindAllMarkers(object = pbmc, thresh.use = 0.25,  min.pct = 0.1,test.use="wilcox")
# print pippo ordered by the p_val column
pippo = pippo[order(pippo$p_val),]
head(pippo)
