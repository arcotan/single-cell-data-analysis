library(dplyr)
library(Seurat)
library(patchwork)

DATA_DIR = "./dataset/Heart_and_Aorta-10X_P7_4"
LABEL_DIR = "./dataset"
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
  metadata = left_join(data.frame("cell"=cells), metadata)
  metadata = metadata[!is.na(metadata$cluster.ids), ]

  data = data[,colnames(data) %in% metadata$cell]

  true_labels = metadata[,c("cell", "cluster.ids")]
  
  # Print number of cells not mapped to a cluster id
  print(paste("Cells mapped to a cluster: ", sum(!is.na(true_labels$cluster.ids)), "/", nrow(true_labels), sep = ""))
  
  return (list("data" = data, "labels" = true_labels))
} 

# Loading data
experiment_data = load_data(DATA_DIR, LABEL_DIR, CHANNEL)

labels = as.list(experiment_data$labels[['cluster.ids']])
names(labels) <- experiment_data$labels[['cell']]
# Load the PBMC dataset
pbmc.data <- experiment_data$data

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


fingerprints <- (pbmc@assays[[1]]@counts@Dimnames[[2]])
for (f in fingerprints) {
   print(c(labels[[f]], f))
}

View(labels)

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
pbmc <- RunTSNE(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["tsne"]], dims = 1:5, nfeatures = 5)
pbmc[["umap"]]

VizDimLoadings(pbmc, dims=1:2, reduction = "umap")

DimHeatmap(pbmc, dims = 1:2, cells = 500, balanced = TRUE)

# Get clustering data
pbmc <- FindNeighbors(pbmc, dims = 1:2, reduction = "umap")
pbmc <- FindClusters(pbmc, resolution = 0.1, reduction.use='umap', n.start = 6)

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

# compute accuracy
library(fossil)
adj.rand.index(label_df_nna$true_id, label_df_nna$computed_id)

plot1 = DimPlot(pbmc, reduction = "umap")

idents <- Idents(pbmc)
for(i in 1:length(idents)) {
  idents[i][[1]] <- as.numeric(labels[names(idents[i])])

}
pbmc[['true_labels']] <- idents
plot2 = DimPlot(pbmc, reduction = "umap", shape.by = "true_labels", pt.size = 3, raster.dpi = c(3440, 1920))
plot2
library(ggplot2)
ggsave('clusters.png', plot2, dpi=400)


# Differential Expression
library(ggplot2)
library(gridExtra)

DE = function(sobj, ident, experiment_data) {
  markers = FindMarkers(sobj, ident.1 = ident, group.by = "true_labels")
  l = experiment_data$labels
  d = experiment_data$data

  ld = merge(t(data.frame(d)), l, by.x = "row.names", by.y = "cell")
  topn = ld[,c(rownames(head(markers,5)), 'cluster.ids')]

  # Define a function for creating each ggplot object
  create_plot <- function(gene_id, topn) {
    title = paste('C:', ident, 'Gene:', names(data.frame(topn[,c(gene_id,6)]))[[1]], sep=' ')
    ggplot(data.frame(topn[,c(gene_id,6)]), aes(y=data.frame(topn[,c(gene_id,6)])[[1]], x=topn$cluster.ids, group=topn$cluster.ids)) +
      geom_boxplot() +
      xlab("Group") + 
      ylab("Counts") +
      ggtitle(title) +
      theme(plot.title = element_text(size = 10))
  }

  # Create the list of ggplots using lapply and the create_plot function
  plots <- lapply(1:(length(topn)-1), function(gene_id) {
    create_plot(gene_id, topn)
  })


  cnt = data.frame()
  for(id in sort(unique(topn$cluster.ids))) {
    a = apply(topn[topn$cluster.ids == id,1:length(topn)-1], 2, sum)
    cnt = rbind(cnt, a)  
  }
  colnames(cnt) = colnames(topn)[1:length(topn)-1]

  return (list("plots" = plots, "counts" = cnt))
}

DE_ANALISYS = lapply(0:5, function(ident) {
  DE(pbmc, ident, experiment_data)
})

plts = lapply(DE_ANALISYS, function(x) {
  x$plots
})

cnts = lapply(DE_ANALISYS, function(x) {
  x$counts
})

box_plots = grid.arrange(grobs = unlist(plts, recursive=FALSE), ncol=5)

ggsave('de.png', box_plots, dpi=400)

cnts