library(dplyr)
library(Rtsne)
library(Seurat)
library(patchwork)
library(factoextra)
library(DropletUtils)
source("./libraries/utils.R")

DATASETS_FOLDER = './dataset/'
NAME = 'peripheal-blood'

inDataDir  = paste(DATASETS_FOLDER, NAME, sep='')
outDataDir = paste(DATASETS_FOLDER, NAME, '-filtered/', sep='') 
dir.create(outDataDir, recursive = TRUE, showWarnings = FALSE)

# Loading data
allData = Read10X(data.dir = inDataDir, strip.suffix = TRUE)

data = allData[[1]]

antibodyData = data.frame(allData[[2]])

# remove `_TotalSeqB` from column names
rownames(antibodyData) = gsub("_TotalSeqB", "", rownames(antibodyData))

labels = data.frame("cell" = colnames(antibodyData))
labels$go = ''
# labels = vector('list', length = ncol(antibodyData))
# names(labels) = colnames(antibodyData)

for(i in 1:ncol(antibodyData)) {
  cell = colnames(antibodyData)[i]
  cellCounts = antibodyData[[cell]]
  names(cellCounts) <- rownames(antibodyData)
  cellCounts = cellCounts / max(cellCounts)
  cellCounts = sort(cellCounts, decreasing=TRUE)
  for(j in 2:length(cellCounts)) {
    if(cellCounts[j] < 0.95) {
      break
    }
  }
  # labels$go[i] = paste(names(cellCounts)[1:i-1], collapse = '~')
  labels$go[i] = paste(sort(names(cellCounts)[1:j-1]), collapse = '~')
}



# Load the PBMC dataset
pbmc.data <- data

# plot distribution of amount of cells in which each gene is expressed
ggplot(data.frame("sum" = rowSums(pbmc.data > 0)), aes(x=sum)) + geom_histogram()
# plot distribution of sum of counts for each cell
ggplot(data.frame("sum" = colSums(pbmc.data > 0)), aes(x=sum)) + geom_histogram()

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 100)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Set very loose limits for COTAN
pbmc <- subset(pbmc, subset = nFeature_RNA > 0 & nFeature_RNA < 5000 & percent.mt < 5)

# Write pre-processed data
data_to_write = GetAssayData(object = pbmc, assay = "RNA", slot = "data")
write.csv(data_to_write, paste(outDataDir, "data.csv", sep=""))
Dir10X = paste(outDataDir, "10X", sep="")
if (!dir.exists(Dir10X)) {
  write10xCounts(Dir10X, data_to_write)
}

filtered_labels = merge(data.frame("cell"=colnames(data_to_write)), labels)
filtered_labels$cluster.ids = as.numeric(as.factor(filtered_labels$go))

antibodyDataFiltered = antibodyData[,filtered_labels$cell]

antibodyDataFiltered[1:10]
boxplot(t(antibodyDataFiltered)[,1:16])
boxplot(t(antibodyDataFiltered)[,17:32])
fviz_nbclust(t(antibodyDataFiltered),kmeans,"wss",k.max=30)
km=kmeans(t(antibodyDataFiltered),centers=10)
tsne=Rtsne(t(antibodyDataFiltered),perplexity=30)
plot(tsne$Y,pch='*',col=km$cluster)

write.csv(filtered_labels[,c('cell', 'cluster.ids')], paste(outDataDir, "labels.csv", sep=""), row.names = FALSE)

mapping = data.frame("go"=unique(filtered_labels$go), "id"=unique(filtered_labels$cluster.ids))
# sort mapping by id
mapping = mapping[order(mapping$id),]

# Create false GO mapping to ensure compatibility with other datasets
write.csv(mapping, paste(outDataDir, "mapping.csv", sep=""), row.names = FALSE)
