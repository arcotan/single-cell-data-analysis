# TODO import cose
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

# In input: a dataframe with two cluster id columns, named 'true_id' and 'computed_id'
# labels of these columns will be renamed in order to get the best match between clusters.
# The confusion matrix of the clustering is also added to the output:
# rows of the matrix are related to predictions, columns are related to the ground truth
align_clusters = function(label_dataframe) {
  label_dataframe_nna = label_dataframe[!is.na(label_dataframe$true_id) & !is.na(label_dataframe$computed_id), ]
  confusion_matrix = table(label_dataframe_nna$computed_id, label_dataframe_nna$true_id) # row predictions, col ground truth
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
  return (list("confusion_matrix" = confusion_matrix, 
               "label_dataframe" = label_dataframe, 
               "permutation_true" = permutation_true, 
               "permutation_computed" = permutation_computed))
}



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
  true_labels$cluster.ids <- true_labels$cluster.ids + 1
  
  # Print number of cells not mapped to a cluster id
  print(paste("Cells mapped to a cluster: ", sum(!is.na(true_labels$cluster.ids)), "/", nrow(true_labels), sep = ""))
  
  return (list("data" = data, "labels" = true_labels))
} 

write_clustering = function(outdir, tag, label_df, cell_col, cluster_col) {
  to_write = label_df[c(cell_col, cluster_col)]
  colnames(to_write)[colnames(to_write) == cell_col] = "cell"
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"
  write.csv(to_write, paste(outdir, "/", "clustering_", tag, ".csv", sep=""), row.names = FALSE)
}

write_markers = function(outdir, tag, marker_df, gene_col, cluster_col) {
  to_write = marker_df[c(gene_col, cluster_col)]
  colnames(to_write)[colnames(to_write) == gene_col] = "gene"
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"
  write.csv(to_write, paste(outdir, "/", "markers_", tag, ".csv", sep=""), row.names = FALSE)
}