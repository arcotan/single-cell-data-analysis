
library(dplyr)

source("utils.R")

TOOL_TAGS = c("seurat")
DATASET_TAGS = c("10X_P7_4")
LABEL_TAG_TO_DIR = list("10X_P7_4" = "./dataset/tabulamuris/")


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
        full_join(score_data, scores_to_add)# TODO only works due to empty intersections
        # merge markers
        markers_to_add = cur_tool_data$markers
        markers_to_add$tool = tool
        full_join(marker_data, markers_to_add)# TODO only works due to empty intersections
      }
    }
  }
  
  # read cluster ids of cells of the dataset and merge 
  # drops cells not used in any clustering
  true_ids = load_dataset_labels(LABEL_TAG_TO_DIR[[dataset_tag]], dataset_tag)
  colnames(true_ids)[colnames(true_ids) == "cluster.ids"] <- "true_labels"
  label_data = left_join(label_data, true_ids)
  
  return (list("labels" = label_data, "scores" = score_data, "markers" = marker_data))
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


global_data = collect_data(DATASET_TAGS, TOOL_TAGS)


# TODO align clusters
# TODO recompute scores
# TODO plot clustering and save results to png

# TODO print NA count

# TODO print marker intersection between tools
# library(ggvenn)
# marker_df = read.csv("./results/seurat/markers_10X_P7_4.csv")
# top_5 = marker_df[marker_df$rank <= 5,]
# ggvenn(split(marker_df$gene, marker_df$cluster), show_elements = TRUE)
# any(table(marker_df$gene, marker_df$cluster) > 1)