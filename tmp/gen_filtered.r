library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)

source("./libraries/utils.R")

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