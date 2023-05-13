
library(dplyr)
library(ggplot2)

source("utils.R")

#TOOL_TAGS = c('scvi', 'scanpy', 'seurat', 'scvitools', 'COTAN')
TOOL_TAGS = c('seurat')
# TODO change tags
#DATASET_TAGS= c('tabula-muris-heart', 'tabula-muris-marrow_P7_3', 'peripheal-blood', 'kumar-4-hard', 'kumar-8-hard')
DATASET_TAGS= c('tabula-muris-heart')

RESULT_DIR = "./results/"
AGGREGATE_RESULT_DIR = paste(RESULT_DIR, "aggregate/", sep="")
DATASET_DIR = "./dataset/"

DATASET_TAG_TO_TRUE_LABEL_DIR = list()
DATASET_TAG_TO_MAPPING_DIR = list()
DATASET_TAG_TO_FILTERED_GE_DIR = list()
for (tag in DATASET_TAGS) {
  DATASET_TAG_TO_MAPPING_DIR[[tag]] = paste(DATASET_DIR, tag, "-filtered/", sep="")
  DATASET_TAG_TO_FILTERED_GE_DIR[[tag]] = paste(DATASET_DIR, tag, "-filtered/10X/", sep="")
}

DATASET_TAG_TO_TRUE_LABEL_DIR = DATASET_TAG_TO_MAPPING_DIR


read_single_data = function(tool_tag, dataset_tag) {
  label_file = paste(RESULT_DIR, dataset_tag, "/", tool_tag, "/clustering_labels", ".csv", sep="")
  score_file = paste(RESULT_DIR, dataset_tag, "/", tool_tag, "/clustering_scores", ".csv", sep="")
  marker_file = paste(RESULT_DIR, dataset_tag, "/", tool_tag, "/markers", ".csv", sep="")
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
        score_data = full_join(score_data, scores_to_add)
        # merge markers
        markers_to_add = cur_tool_data$markers
        markers_to_add$tool = tool
        marker_data = full_join(marker_data, markers_to_add)
      }
    }
  }
  
  # read cluster ids of cells of the dataset and merge 
  # drops cells not used in any clustering
  true_ids = read.csv(paste(DATASET_TAG_TO_TRUE_LABEL_DIR[[dataset_tag]], "labels.csv", sep=""))

  colnames(true_ids)[colnames(true_ids) == "cluster.ids"] <- "true_labels"
  label_data = merge(label_data, true_ids)

  if (!to_init) {
    # align each clustering with true labels
    for (label in colnames(label_data)[-1]) {
      tool = substr(label, 1, nchar(label)-6)
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
  
  return (list("labels" = label_data, "scores" = score_data, "markers" = marker_data))
}

collect_data = function(dataset_tag_list, tool_tag_list, write_aggregate = TRUE) {
  datasets_aggregate_data_list = list()
  for (tag in dataset_tag_list) {
    dataset_data = read_dataset_data(tool_tag_list, tag)
    if (write_aggregate) {
      if (!dir.exists(paste(AGGREGATE_RESULT_DIR, tag, sep=""))) {
        dir.create(paste(AGGREGATE_RESULT_DIR, tag, sep=""))
      }
      write.csv(dataset_data$labels, paste(AGGREGATE_RESULT_DIR, tag, "/labels.csv", sep=""), row.names = FALSE)
      write.csv(dataset_data$scores, paste(AGGREGATE_RESULT_DIR, tag, "/scores.csv", sep=""), row.names = FALSE)
      write.csv(dataset_data$markers, paste(AGGREGATE_RESULT_DIR, tag, "/markers.csv", sep=""), row.names = FALSE)
    }
    datasets_aggregate_data_list[[tag]] = dataset_data
  }
  return (datasets_aggregate_data_list)
}

# read data
global_data = collect_data(DATASET_TAGS, TOOL_TAGS)

# print NA count
print("NA count")
for (dataset in DATASET_TAGS) {
  print(paste("Dataset: ", dataset))
  print(colSums(is.na(global_data[[dataset]]$labels[,2:ncol(global_data[[dataset]]$labels)])))
}

# plot clustering and de and save results in eps format
for (dataset in DATASET_TAGS) {
  # load GO mapping
  go_mapping = read.csv(paste(DATASET_TAG_TO_MAPPING_DIR[[dataset]], "mapping.csv", sep=""))
  go_mapping = go_mapping[order(go_mapping$id),]
  pi = go_mapping$go

  print("--------------------------------------")
  print(paste("Clustering results for dataset ", dataset, sep=""))
  pbmc.data <- Read10X(DATASET_TAG_TO_FILTERED_GE_DIR[[dataset]], strip.suffix = TRUE)
  pbmc <- CreateSeuratObject(counts = pbmc.data)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")
  pbmc <- ScaleData(pbmc, features = rownames(pbmc))
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

  columns = colnames(global_data[[dataset]]$labels)
  for (label in columns[-c(1, length(columns))]) {
    tool = substr(label, 1, nchar(label)-6)
    print(paste("Tool: ", tool, sep=""))
    
    # plot clustering
    cur_plot <- seurat_clustering_plot(pbmc, global_data[[dataset]]$labels$cell, pi[global_data[[dataset]]$labels[[label]]])
    ggsave(filename = paste(AGGREGATE_RESULT_DIR, dataset, "/", label, ".png", sep=""), cur_plot)

    # plot de 
    plot_de(pbmc.data, global_data[[dataset]]$markers[global_data[[dataset]]$markers$tool == tool,], "gene", "cluster", global_data[[dataset]]$labels, "cell", label, paste(RESULT_DIR, dataset, "/", tool, "/", sep=""))
  }
  print("--------------------------------------")
}
