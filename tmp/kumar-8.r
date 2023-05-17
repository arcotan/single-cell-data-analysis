library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)
source("./libraries/utils.R")

DATASETS_FOLDER = './dataset/'
NAME = 'kumar-8-hard'
RDS = 'sce_full_SimKumar8hard.rds'

inDataDir  = paste(DATASETS_FOLDER, NAME, sep='')
outDataDir = paste(DATASETS_FOLDER, NAME, '-filtered/', sep='')


# Loading data
data = readRDS(paste(inDataDir, RDS, sep="/"))
counts = assays(data)$counts
# Load the PBMC dataset
pbmc.data <- counts

# plot distribution of amount of cells in which each gene is expressed
ggplot(data.frame("sum" = rowSums(pbmc.data > 0)), aes(x=sum)) + geom_histogram()
# plot distribution of sum of counts for each cell
ggplot(data.frame("sum" = colSums(pbmc.data > 0)), aes(x=sum)) + geom_histogram()

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 100)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Set very loose limits for COTAN
pbmc <- subset(pbmc, subset = nFeature_RNA > 0)

# Write pre-processed data
data_to_write = GetAssayData(object = pbmc, assay = "RNA", slot = "data")
write.csv(data_to_write, paste(outDataDir, "data.csv", sep=""))
Dir10X = paste(outDataDir, "10X", sep="")
if (!dir.exists(Dir10X)) {
  write10xCounts(Dir10X, data_to_write)
}

phenoid = data@colData@listData$phenoid
cells = data@colData@rownames
labels = sort(unique(phenoid))
metadata = data.frame("cell"=cells, "cluster.ids"=match(phenoid, labels), "cell_ontology_class"=phenoid)

filtered_labels = merge(data.frame("cell"=colnames(data_to_write)), metadata)
write.csv(filtered_labels[,c('cell', 'cluster.ids')], paste(outDataDir, "labels.csv", sep=""), row.names = FALSE)


mapping = data.frame("go"=labels, "id"=1:length(labels))
write.csv(mapping, paste(outDataDir, "mapping.csv", sep=""), row.names = FALSE)