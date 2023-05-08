library(monocle3)

source("utils.R")

DATA_DIR = "./filtered_dataset/tabulamuris/data_10X_P7_4"
LABEL_DIR = "filtered_dataset/tabulamuris/"
CHANNEL = "10X_P7_4"

cds <- load_mm_data(
  mat_path = paste(DATA_DIR, "matrix.mtx", sep = "/"), 
  feature_anno_path = paste(DATA_DIR, "genes.tsv", sep = "/"),
  cell_anno_path = paste(DATA_DIR, "barcodes.tsv", sep = "/")
)
true_labels = read.csv(paste(LABEL_DIR, "labels_", CHANNEL, ".csv", sep=""))

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds, reduction_method = "PCA")

cds <- cluster_cells(cds, resolution=0.08, reduction_method = "PCA")

plot_cells(cds, reduction_method = "PCA")

label_df = merge(true_labels, data.frame(clusters(cds, reduction_method = "PCA")), by.x = "cell", by.y = 0)
names(label_df)[2:3] <- c("true_id", "computed_id")

label_df$computed_id = as.numeric(label_df$computed_id)

res = align_clusters(label_df)

res$confusion_matrix

label_df <- res$label_dataframe

marker_test_res <- top_markers(cds,
                               group_cells_by="cluster",
                               reduction_method = "PCA",
                               cores=8)

marker_test_res$cell_group <- res$permutation_computed[as.numeric(marker_test_res$cell_group)]

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(20, pseudo_R2)

#TODO plot markers

write_clustering(LABEL_DIR, paste(CHANNEL, "_Monocle", sep=""), label_df, "cell", "computed_id")
write_markers(LABEL_DIR, paste(CHANNEL, "_Monocle", sep=""), top_specific_markers, "gene_id", "cell_group")

