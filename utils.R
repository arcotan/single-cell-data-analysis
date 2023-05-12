# TODO import cose
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(NMF)

# In input: a dataframe with two cluster id columns
# labels of the computed ids column will be renamed in order to get the best match between clusters.
# The confusion matrix of the clustering is also added to the output:
# rows of the matrix are related to predictions, columns are related to the ground truth
align_clusters = function(label_dataframe, true_id_col, computed_id_col) {
  label_dataframe_nna = label_dataframe[!is.na(label_dataframe[[true_id_col]]) & !is.na(label_dataframe[[computed_id_col]]), ]
  confusion_matrix = table(label_dataframe_nna[[computed_id_col]], label_dataframe_nna[[true_id_col]]) # row predictions, col ground truth
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
  
  # apply permutation to computed ids
  pi = permutation_true[order(permutation_computed)]
  label_dataframe[[computed_id_col]] = pi[label_dataframe[[computed_id_col]]]

  # return confusion matrix, TODO più efficiente senza ricalcolare crosstable
  return (list("confusion_matrix" = table(label_dataframe[[computed_id_col]], label_dataframe[[true_id_col]]), 
               "label_dataframe" = label_dataframe, 
               "permutation_computed" = pi))
}


load_dataset_labels <- function(label_dir, channel) {
  # Load cell metadata
  metadata = read.csv(file = paste(label_dir, "annotations_droplet.csv", sep = "/"))
  
  # Filter data to use only data for current dataset
  metadata = metadata[metadata$channel == channel,]
  metadata$cell = substr(metadata$cell, 10, 25)
  metadata$cluster.ids = metadata$cell_ontology_class
  metadata = metadata[c("cell", "cluster.ids")]
  metadata[metadata$cluster.ids == "",]$cluster.ids <- "unknown cluster"
  labels <- sort(unique(metadata$cluster.ids))
  label_map = list()
  for (i in 1:length(labels)) {
    label_map[[labels[i]]] <- i
  }
  metadata$cluster.ids <- as.numeric(label_map[metadata$cluster.ids])
  
  return (list("metadata" = metadata, "mapping" = data.frame("go" = labels, "id" = c(1:length(labels)))))
}

# loads gene expression matrix associated to a channel and gets the cluster label for each cell,
# if metadata for a cell is not found, its cluster id will be NA
load_data <- function(data_dir, label_dir, channel) {
  # Load gene expression matrix
  data = Read10X(data.dir = data_dir, strip.suffix = TRUE)
  
  # Load cell metadata
  metadata_list = load_dataset_labels(label_dir, channel)
  metadata = metadata_list$metadata
  
  # Get cluster labels
  cells = colnames(data)
  true_labels = left_join(data.frame("cell"=cells), metadata)
  
  
  # Print number of cells not mapped to a cluster id
  print(paste("Cells mapped to a cluster: ", sum(!is.na(true_labels$cluster.ids)), "/", nrow(true_labels), sep = ""))
  
  return (list("data" = data, "labels" = true_labels, "mapping" = metadata_list$mapping))
} 

# returns clustering plot with pca
seurat_clustering_plot = function(seurat_obj, cell_col, label_col) {
  pi = order(label_col)
  cell_col = cell_col[pi]
  label_col = label_col[pi]
  seurat_obj <- SetIdent(seurat_obj, cells = cell_col, label_col)
  return (DimPlot(seurat_obj, reduction = "pca"))
}

clustering_simple_scores = function(label_df, computed_label, true_label) {
  confusion_matrix = table(label_df[[computed_label]], label_df[[true_label]])
  res_overlap = sum(diag(confusion_matrix)) / sum(confusion_matrix)
  res_entropy = entropy(confusion_matrix)
  res_purity = purity(confusion_matrix)
  return (data.frame("entropy" = res_entropy, "purity" = res_purity, "accuracy" = res_overlap))
}

#TODO verifica nan
#clustering_scores = function(label_df, computed_label, true_label, distance_matrix) {
#  simple_scores = clustering_simple_scores(label_df, computed_label, true_label)
#  simple_scores$silhouette = mean(silhouette(label_df[[computed_label]], distance_matrix)[,3])
#  return (simple_scores)
#}

clustering_complex_scores = function(label_df, computed_label, distance_matrix) {
  return (data.frame("silhouette" = mean(silhouette(label_df[[computed_label]], distance_matrix)[,3])))
}

write_clustering = function(outdir, label_df, cell_col, cluster_col, distance_matrix) {
  # compute and write scores
  cscores = clustering_complex_scores(label_df, cluster_col, distance_matrix)
  write.csv(cscores, paste(outdir, "/clustering_scores.csv", sep=""), row.names = FALSE)
  
  # write labels
  to_write = label_df[c(cell_col, cluster_col)]
  colnames(to_write)[colnames(to_write) == cell_col] = "cell"
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"
  write.csv(to_write, paste(outdir, "/clustering_labels.csv", sep=""), row.names = FALSE)
}

write_markers = function(outdir, marker_df, gene_col, cluster_col, order_by_col, is_higher_better, top_k) {
  is_higher_better = is_higher_better*1
  marker_df = marker_df %>%
    group_by(get(cluster_col)) %>%
    slice_max(n = top_k, order_by = (is_higher_better*get(order_by_col))) %>% 
    mutate(rank = row_number(), ties.method = "first")
  marker_df = data.frame(marker_df)
  to_write = marker_df[c(gene_col, cluster_col, "rank")]
  colnames(to_write)[colnames(to_write) == gene_col] = "gene"
  colnames(to_write)[colnames(to_write) == cluster_col] = "cluster"
  write.csv(to_write, paste(outdir, "/markers.csv", sep=""), row.names = FALSE)
}

plot_de = function(expression_matrix, marker_df, gene_col, marker_df_cluster_col, cluster_df, cell_col, cluster_df_cluster_col, out_dir) {
  DE = function(expression_matrix, marker_df, gene_col, marker_df_cluster_col, cur_ident, IDENTS, cell_col, cluster_col) {
    print(paste('Class:', cur_ident, sep=' '))
    markers = marker_df[marker_df[[marker_df_cluster_col]] == cur_ident, ]
    print(head(markers, 5))

    l = IDENTS
    d = expression_matrix
    
    ld = merge(t(data.frame(d)), l, by.x = "row.names", by.y = cell_col)
    colnames(ld)[length(ld)] = 'cluster.ids'
    
    topn = ld[, c(head(markers,5)[[gene_col]], 'cluster.ids')]
    
    # Define a function for creating each ggplot object
    create_plot <- function(gene_id, topn) {
      title = paste('C:', cur_ident, 'Gene:', names(data.frame(topn[,c(gene_id,6)]))[[1]], sep=' ')
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

  DE_ANALISYS = lapply(sort(unique(label_df[[cluster_df_cluster_col]])), function(ident) {
    DE(expression_matrix, marker_df, gene_col, marker_df_cluster_col, ident, label_df, cell_col, cluster_df_cluster_col)
  })
  
  plts = lapply(DE_ANALISYS, function(x) {
    x$plots
  })
  
  cnts = lapply(DE_ANALISYS, function(x) {
    x$counts
  })
  
  box_plots = grid.arrange(grobs = unlist(plts, recursive=FALSE), ncol=5)
  
  ggsave(paste(out_dir, "/de.png", sep=""), box_plots, dpi=400)
}

