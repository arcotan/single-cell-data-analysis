
library(dplyr)
library(enrichR)
library(ggplot2)

source("./libraries/utils.R")
source("./libraries/venn.R")
source("./libraries/enrichment_lists.R")

TOOL_TAGS = c('monocle', 'scanpy', 'seurat', 'scvi', 'COTAN')
DATASET_TAGS= c('peripheal-blood', 'tabula-muris-heart', 'tabula-muris-marrow_P7_3', 'zheng-4', 'zheng-8')

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

DATASET_TAG_TO_ENRICHR_DB = list("tabula-muris-heart" = "Tabula_Muris",
                                 "tabula-muris-marrow_P7_3" = "Tabula_Muris",
                                 "peripheal-blood" = "Tabula_Sapiens") # TODO: add the correct database

DATASET_TAG_TO_GENES_TO_ENRICH_DIR = list()
DATASET_TAG_TO_ENRICHER_DIR = list()
for (tag in DATASET_TAGS) {
  DATASET_TAG_TO_GENES_TO_ENRICH_DIR[[tag]] = paste(AGGREGATE_RESULT_DIR, "/", tag, "/genes_to_enrich/", sep="")
  DATASET_TAG_TO_ENRICHER_DIR[[tag]] = paste(AGGREGATE_RESULT_DIR, "/", tag, "/enrichr_data/", sep="")
}

read_single_data = function(tool_tag, dataset_tag) {
  label_file = paste(RESULT_DIR, dataset_tag, "/", tool_tag, "/clustering_labels", ".csv", sep="")
  score_file = paste(RESULT_DIR, dataset_tag, "/", tool_tag, "/clustering_scores", ".csv", sep="")
  marker_file = paste(RESULT_DIR, dataset_tag, "/", tool_tag, "/markers", ".csv", sep="")
  res = list()
  if (file.exists(label_file)) {
    res$labels = read.csv(label_file)
  }
  if (file.exists(score_file)) {
    temp_scores = read.csv(score_file)
    if (is.null(temp_scores$accuracy)) {
      temp_scores$accuracy = NA
    }
    if (is.null(temp_scores$entropy)) {
      temp_scores$entropy = NA
    }
    if (is.null(temp_scores$purity)) {
      temp_scores$purity = NA
    }
    if (is.null(temp_scores$silhouette)) {
      temp_scores$silhouette = NA
    }
    res$scores = temp_scores
  }
  else {
    res$scores = data.frame("accuracy" = NA, "entropy" = NA, "purity" = NA, "silhouette" = NA)
  }
  if (file.exists(marker_file)) {
    res$markers = read.csv(marker_file)
  }
  return (res)
}

collect_dataset_data = function(tool_tag_list, dataset_tag, compute_missing_scores = TRUE, filtered_datasets_dir_map = NULL) {
  print(paste("Collecting data for dataset", dataset_tag, "...", sep=" "))
  to_init = TRUE
  label_data = NULL
  score_data = NULL
  marker_data = NULL
  for (tool in tool_tag_list) {
    # read data of current tool
    cur_tool_data = read_single_data(tool, dataset_tag)
    if (!is.null(cur_tool_data$markers) && !is.null(cur_tool_data$labels)) {
      if (to_init) {
        label_data = cur_tool_data$labels
        colnames(label_data)[colnames(label_data) == "cluster"] <- paste(tool ,"_label", sep="")
        if (!is.null(cur_tool_data$scores)) {
          score_data = cur_tool_data$scores
          score_data$tool = tool
        }
        else {
          score_data = data.frame("tool" = tool, "accuracy" = NA, "entropy" = NA, "purity" = NA)
        }
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
        if (!is.null(cur_tool_data$scores)) {
          scores_to_add = cur_tool_data$scores
          scores_to_add$tool = tool
          score_data = full_join(score_data, scores_to_add)
        }
        else {
          score_data = rbind(score_data, data.frame("tool" = tool, "accuracy" = NA, "entropy" = NA, "purity" = NA))
        }
        # merge markers
        markers_to_add = cur_tool_data$markers
        markers_to_add$tool = tool
        marker_data = rbind(marker_data, markers_to_add)
      }
    }
  }

  if (!to_init) {
    # read cluster ids of cells of the dataset and merge 
    # drops cells not used in any clustering
    true_ids = read.csv(paste(DATASET_TAG_TO_TRUE_LABEL_DIR[[dataset_tag]], "labels.csv", sep=""))

    colnames(true_ids)[colnames(true_ids) == "cluster.ids"] <- "true_labels"
    label_data = merge(label_data, true_ids)
    
    # align each clustering with true labels
    for (label in colnames(label_data)[-1]) {
      tool = substr(label, 1, nchar(label)-6)
      alignment = align_clusters(label_data, "true_labels", label)
      label_data[[label]] <- alignment$label_dataframe[[label]]
      marker_data[marker_data$tool == tool,]$cluster <- alignment$permutation_computed[marker_data[marker_data$tool == tool,]$cluster]
    }
    if (compute_missing_scores) {
      # compute missing clustering scores (for scanpy and scvi only silhouette should not have NA at this point)
      for (i in 1:nrow(score_data)) {
        cur_info = score_data[i, ]
        if (is.na(cur_info$accuracy) || is.na(cur_info$entropy) || is.na(cur_info$purity)) {
          scores_to_add = clustering_simple_scores(label_data, paste(cur_info$tool, "_label", sep=""), "true_labels")
          score_data[i, "accuracy"] <- scores_to_add$accuracy
          score_data[i, "entropy"] <- scores_to_add$entropy
          score_data[i, "purity"] <- scores_to_add$purity
        }
        if (is.na(cur_info$silhouette) && !is.null(filtered_datasets_dir_map)) {
          ge = t(Read10X(DATASET_TAG_TO_FILTERED_GE_DIR[[dataset_tag]], strip.suffix = TRUE))
          ge = ge[rownames(ge) %in% label_data$cell,]
          distance_matrix <- dist(ge)
          score_data[i, "silhouette"] <- clustering_complex_scores(label_data, paste(cur_info$tool, "_label", sep=""), distance_matrix)$silhouette
        }
      }
    }
  }
  
  return (list("labels" = label_data, "scores" = score_data, "markers" = marker_data))
}

collect_data = function(dataset_tag_list, tool_tag_list, write_aggregate = TRUE, compute_missing_scores = TRUE, filtered_datasets_dir_map = NULL) {
  datasets_aggregate_data_list = list()
  for (tag in dataset_tag_list) {
    dataset_data = collect_dataset_data(tool_tag_list, tag, compute_missing_scores, filtered_datasets_dir_map)
    if (!is.null(dataset_data$labels) || !is.null(dataset_data$scores) || !is.null(dataset_data$markers)) {
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
  }
  return (datasets_aggregate_data_list)
}

# read data
global_data = collect_data(DATASET_TAGS, TOOL_TAGS, filtered_datasets_dir_map = DATASET_TAG_TO_FILTERED_GE_DIR)

dataset_found = c()
for (dataset in DATASET_TAGS) {
  if (!is.null(global_data[[dataset]])) {
    dataset_found = c(dataset_found, dataset)
  }
}

# print NA count
print("NA count")
for (dataset in dataset_found) {
  print(paste("Dataset: ", dataset))
  print(colSums(is.na(global_data[[dataset]]$labels[,2:ncol(global_data[[dataset]]$labels)])))
}

# plot clustering and de and save results in eps format
for (dataset in dataset_found) {
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
    ggsave(filename = paste(AGGREGATE_RESULT_DIR, dataset, "/", label, ".eps", sep=""), cur_plot)

    # plot de 
    plot_de(pbmc.data, global_data[[dataset]]$markers[global_data[[dataset]]$markers$tool == tool,], "gene", "cluster", global_data[[dataset]]$labels, "cell", label, paste(RESULT_DIR, dataset, "/", tool, "/", sep=""))
  }
  print("--------------------------------------")
}

# Venn diagram
for (dataset in dataset_found) {
  plot_venn(global_data[[dataset]]$markers, paste(AGGREGATE_RESULT_DIR, dataset, "/", sep=""))
}

# write enrichment results
for (dataset in dataset_found) {
  if (!dir.exists(DATASET_TAG_TO_GENES_TO_ENRICH_DIR[[dataset]])) {
    dir.create(DATASET_TAG_TO_GENES_TO_ENRICH_DIR[[dataset]])
  }
  write_markers_to_enrich(global_data[[dataset]]$markers, DATASET_TAG_TO_GENES_TO_ENRICH_DIR[[dataset]])
  cur_enrichr_db <- DATASET_TAG_TO_ENRICHR_DB[[dataset]]

  if (!is.null(cur_enrichr_db)) {
    if (!dir.exists(DATASET_TAG_TO_ENRICHER_DIR[[dataset]])) {
      dir.create(DATASET_TAG_TO_ENRICHER_DIR[[dataset]])
    }
    write_enrichment_result(global_data[[dataset]]$markers, DATASET_TAG_TO_ENRICHER_DIR[[dataset]], cur_enrichr_db)
  }
}

