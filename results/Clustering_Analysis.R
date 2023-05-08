
library(dplyr)
library(ggplot2)

source("utils.R")

TOOL_TAGS = c("seurat")
DATASET_TAGS = c("10X_P7_4")
LABEL_TAG_TO_LABEL_DIR = list("10X_P7_4" = "./dataset/tabulamuris/")
LABEL_TAG_TO_FILTERED_GE_DIR = list("10X_P7_4" = "./filtered_dataset/tabulamuris/data_10X_P7_4")


read_single_data = function(tool_tag, dataset_tag) {
  label_file = paste("./results/", tool_tag, "/clustering_labels_", dataset_tag, ".csv", sep="")
  score_file = paste("./results/", tool_tag, "/clustering_scores_", dataset_tag, ".csv", sep="")
  marker_file = paste("./results/", tool_tag, "/markers_", dataset_tag, ".csv", sep="")
  if (file.exists(label_file) && file.exists(score_file) && file.exists(marker_file)) {
    label_data = read.csv(label_file)
    score_data = read.csv(score_file)
    marker_data = read.csv(marker_file)
    return (list("labels" = label_data, "scores" = score_data, "markers" = marker_data))
  }
  else {
    return (NULL)
  }
}

read_dataset_data = function(tool_tag_list, dataset_tag) {
  to_init = TRUE
  label_data = NULL
  score_data = NULL
  marker_data = NULL
  for (tool in tool_tag_list) {
    # read data of current tool
    cur_tool_data = read_single_data(tool, dataset_tag)
    if (!is.null(cur_tool_data)) {
      if (to_init) {
        label_data = cur_tool_data$labels
        colnames(label_data)[colnames(label_data) == "cluster"] <- paste(tool ,"_label", sep="")
        score_data = cur_tool_data$scores
        score_data$tool = tool
        marker_data = cur_tool_data$markers
        marker_data$tool = tool
        to_init = FALSE
      }
      else {
        # merge clustering labels
        labels_to_add = cur_tool_data$labels
        colnames(labels_to_add)[colnames(labels_to_add) == "cluster"] <- paste(tool ,"_label", sep="")
        label_data = full_join(label_data, labels_to_add)
        # merge clustering scores
        scores_to_add = cur_tool_data$scores
        scores_to_add$tool = tool
        score_data = full_join(score_data, scores_to_add)# TODO only works due to empty intersections
        # merge markers
        markers_to_add = cur_tool_data$markers
        markers_to_add$tool = tool
        marker_data = full_join(marker_data, markers_to_add)# TODO only works due to empty intersections
      }
    }
  }
  
  # read cluster ids of cells of the dataset and merge 
  # drops cells not used in any clustering
  metadata = load_dataset_labels(LABEL_TAG_TO_LABEL_DIR[[dataset_tag]], dataset_tag)
  true_ids = metadata$metadata
  colnames(true_ids)[colnames(true_ids) == "cluster.ids"] <- "true_labels"
  label_data = left_join(label_data, true_ids)
  
  if (!to_init) {
    # align each clustering with true labels
    for (tool in TOOL_TAGS) {
      label <- paste(tool, "_label", sep="")
      alignment = align_clusters(label_data, "true_labels", label)
      label_data[[label]] <- alignment$label_dataframe[[label]]
      marker_data[marker_data$tool == tool,]$cluster <- alignment$permutation_computed[marker_data[marker_data$tool == tool,]$cluster]
    }
    # compute missing clustering scores (for scanpy and scvi only silhouette should not have NA at this point)
    for (i in 1:nrow(score_data)) {
      cur_info = score_data[i, ]
      scores_to_add = clustering_simple_scores(label_data, paste(cur_info$tool, "_label", sep=""), "true_labels")
      score_data[i, "accuracy"] <- scores_to_add$accuracy
      score_data[i, "entropy"] <- scores_to_add$entropy
      score_data[i, "purity"] <- scores_to_add$purity
    }
  }
  
  return (list("labels" = label_data, "scores" = score_data, "markers" = marker_data, "mapping" = metadata$mapping))
}

collect_data = function(dataset_tag_list, tool_tag_list, write_aggregate = TRUE) {
  datasets_aggregate_data_list = list()
  for (tag in dataset_tag_list) {
    dataset_data = read_dataset_data(tool_tag_list, tag)
    if (write_aggregate) {
      write.csv(dataset_data$labels, paste("./results/aggregate/", tag, "_labels.csv"), row.names = FALSE)
      write.csv(dataset_data$scores, paste("./results/aggregate/", tag, "_scores.csv"), row.names = FALSE)
      write.csv(dataset_data$markers, paste("./results/aggregate/", tag, "_markers.csv"), row.names = FALSE)
    }
    datasets_aggregate_data_list[[tag]] = dataset_data
  }
  return (datasets_aggregate_data_list)
}

# TODO $ non funziona se il tag inizia con un numero
global_data = collect_data(DATASET_TAGS, TOOL_TAGS)
# table(global_data[[DATASET_TAGS[1]]]$labels$scanpy_label, global_data[[DATASET_TAGS[1]]]$labels$scvi_label)
# View(global_data[[DATASET_TAGS[1]]]$scores)

# print NA count
print("NA count")
for (dataset in DATASET_TAGS) {
  print(paste("Dataset: ", dataset))
  print(colSums(is.na(global_data[[dataset]]$labels[,2:ncol(global_data[[dataset]]$labels)])))
}

# plot clustering and save results to eps
for (dataset in DATASET_TAGS) {
  print("--------------------------------------")
  print(paste("Clustering results for dataset ", dataset, sep=""))
  pbmc.data <- Read10X(LABEL_TAG_TO_FILTERED_GE_DIR[[dataset]], strip.suffix = TRUE)
  pbmc <- CreateSeuratObject(counts = pbmc.data)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
  pbmc <- ScaleData(pbmc, features = rownames(pbmc))
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  for (label in colnames(global_data[[dataset]]$labels)[-1]) {
    print(label)
    print(dataset)
    cur_plot <- seurat_clustering_plot(pbmc, global_data[[dataset]]$labels$cell, global_data[[dataset]]$labels[[label]])
    ggsave(filename = paste("./results/aggregate/", dataset, "_", label, ".png", sep=""), cur_plot)
  }
  print("--------------------------------------")
}


# TODO print marker intersection between tools
# library(ggvenn)
# marker_df = read.csv("./results/seurat/markers_10X_P7_4.csv")
# top_5 = marker_df[marker_df$rank <= 5,]
# ggvenn(split(marker_df$gene, marker_df$cluster), show_elements = TRUE)
# any(table(marker_df$gene, marker_df$cluster) > 1)