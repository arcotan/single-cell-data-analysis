# takes the path to the csv with aggregated markers and the cluster id
get_common_markers = function(markers_path, cluster_id) {
  markers = read.csv(markers_path)
  cluster_markers = markers[markers$cluster == cluster_id,]
  common_markers = Reduce(intersect,split(cluster_markers$gene, cluster_markers$tool))
  return (common_markers)
}

# takes the path to the csv with aggregated markers and the cluster id
get_specific_markers = function(markers_path, cluster_id) {
  markers = read.csv(markers_path)
  cluster_markers = markers[markers$cluster == cluster_id,]
  tools = unique(cluster_markers$tool)
  tools_markers = split(cluster_markers$gene, cluster_markers$tool)
  specific_markers <- list()
  for (t1 in tools) {
    t1_markers = tools_markers[[t1]]
    for (t2 in tools) {
      if (t1 != t2) {
        t1_markers = setdiff(t1_markers, intersect(t1_markers, tools_markers[[t2]]))
      }
    }
    if (length(t1_markers) != 0) {
      specific_markers[[t1]] <- t1_markers
    }
  }
  return (specific_markers)
}

# takes the path to the csv with aggregated markers
write_markers_to_enrich = function(markers_path) {
  markers = read.csv(markers_path)
  cluster_ids = unique(markers$cluster)
  for (cid in cluster_ids) {
    common_markers = get_common_markers(markers_path, cid)
    write.table(
      common_markers,
      paste("cluster",cid,"_intersection_to_enrich.txt",sep=""),
      row.names=FALSE,
      quote=FALSE,
      col.names=FALSE)
    specific_markers = get_specific_markers(markers_path, cid)
    for (tool in names(specific_markers)) {
      write.table(
        specific_markers[[tool]],
        paste("cluster",cid,"_",tool,"_to_enrich.txt",sep=""),
        row.names=FALSE,
        quote=FALSE,
        col.names=FALSE)
    }
  }
}