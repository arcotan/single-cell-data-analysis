library(Seurat)

# get ids of clusters bigger than 40 cells
mapping = read.csv('../dataset/PBMC1-Filtered/raw/celltypist_mapping.csv')
counts = read.csv('../dataset/PBMC1-Filtered/raw/celltypist_annotation_counts.csv')
mapping_counts = merge(mapping, counts, by.x = "go", by.y = "cluster.ids")
mapping_counts = subset(mapping_counts, count > 40)
clusters_ids_to_keep = mapping_counts$id

# get barcodes of cells in clusters bigger than 40 cells
labels = read.csv('../dataset/PBMC1-Filtered/raw/celltypist_labels.csv')
labels = subset(labels, cluster.ids %in% clusters_ids_to_keep)
barcodes_to_keep = labels$cell

# read the dataset
dataset = Read10X(data.dir = '../dataset/PBMC1-Filtered/10x', strip.suffix = TRUE)
# keep only cells in clusters bigger than 40 cells
dataset = dataset[, colnames(dataset) %in% barcodes_to_keep]