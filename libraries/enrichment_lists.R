library(enrichR)

# takes the path to the csv with aggregated markers and the cluster id
get_common_markers = function(markers, cluster_id) {
  cluster_markers = markers[markers$cluster == cluster_id,]
  common_markers = Reduce(intersect,split(cluster_markers$gene, cluster_markers$tool))
  return (common_markers)
}

# takes the path to the csv with aggregated markers and the cluster id
get_specific_markers = function(markers, cluster_id) {
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

# takes the path to the csv with aggregated markers and the cluster id
get_not_fully_intersected_markers = function(markers, cluster_id) {
  intersection = get_specific_markers(markers, cluster_id)
  cluster_markers = markers[markers$cluster == cluster_id,]
  tools = unique(cluster_markers$tool)
  tools_markers = split(cluster_markers$gene, cluster_markers$tool)
  specific_markers <- list()
  for (t1 in tools) {
    t1_markers = tools_markers[[t1]]
    t1_markers = setdiff(t1_markers, intersect(t1_markers, intersection))
    if (length(t1_markers) != 0) {
      specific_markers[[t1]] <- t1_markers
    }
  }
  return (specific_markers)
}

# takes the path to the csv with aggregated markers
write_markers_to_enrich = function(markers, out_dir) {
  cluster_ids = unique(markers$cluster)
  for (cid in cluster_ids) {
    # Print Cluster Id 
    common_markers = get_common_markers(markers, cid)
    write.table(
      common_markers,
      paste(out_dir,"cluster",cid,"_intersection_to_enrich.txt",sep=""),
      row.names=FALSE,
      quote=FALSE,
      col.names=FALSE)
    specific_markers = get_specific_markers(markers, cid)
    for (tool in names(specific_markers)) {
      write.table(
        specific_markers[[tool]],
        paste(out_dir,"cluster",cid,"_",tool,"_to_enrich.txt",sep=""),
        row.names=FALSE,
        quote=FALSE,
        col.names=FALSE)
    }
    specific_markers = get_not_fully_intersected_markers(markers, cid)
    for (tool in names(specific_markers)) {
      write.table(
        specific_markers[[tool]],
        paste(out_dir,"cluster",cid,"_",tool,"_lone_to_enrich.txt",sep=""),
        row.names=FALSE,
        quote=FALSE,
        col.names=FALSE)
    }
  }
}

get_enriched_cell_type = function(gene_list, enrichr_database, enrichr_site = "Enrichr") {
  websiteLive <- getOption("enrichR.live")
  if (websiteLive) {
      setEnrichrSite(enrichr_site)  
  }
  dbs <- c(enrichr_database)
  if (websiteLive) {
      enriched <- enrichr(gene_list, dbs)
  }

  return (enriched[[enrichr_database]])
}

write_enrichment_result = function(markers, out_dir, enrichr_database) {
  cluster_ids = unique(markers$cluster)
  for (cid in cluster_ids) {
    common_markers = get_common_markers(markers, cid)
    print(paste("Cluster Id: ", cid, "Len: ", length(common_markers)))
    if(!length(common_markers) == 0) {
      enriched <- get_enriched_cell_type(common_markers, enrichr_database)
      if(nrow(enriched) > 0) {
        write.csv(
          enriched,
          paste(out_dir,"cluster",cid,"_enriched_intersection.csv",sep=""),
          row.names=FALSE)
        plotEnrich(enriched, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        ggsave(paste(out_dir,"cluster",cid,"_enriched_intersection_plot.png",sep=""))
      }
    }
    specific_markers = get_specific_markers(markers, cid)
    # View(specific_markers[['scanpy', drop=FALSE]])
    for (tool in names(specific_markers)) {
      print(paste("Tool : ", tool, "Cluster Id: ", cid)) 
      enriched <- get_enriched_cell_type(specific_markers[[tool, drop=FALSE]], enrichr_database)
      print(paste("Enriched: ", length(enriched$Term)))
      if(!length(enriched$Term) == 0) {
        write.csv(
          enriched,
          paste(out_dir,"cluster",cid,"_",tool,"_lone_enriched.csv",sep=""),
          row.names=FALSE)
        print(plotEnrich(enriched, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value"))
        plotEnrich(enriched, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        ggsave(paste(out_dir,"cluster",cid,"_",tool,"_lone_enriched_plot.png",sep=""))
      }
    }
    
    # venn on tool markers minus global intersection
    specific_markers = get_not_fully_intersected_markers(markers, cid)
    for (tool in names(specific_markers)) {
      print(paste("Tool : ", tool, "Cluster Id: ", cid)) 
      enriched <- get_enriched_cell_type(specific_markers[[tool, drop=FALSE]], enrichr_database)
      print(paste("Enriched: ", length(enriched$Term)))
      if(!length(enriched$Term) == 0) {
        write.csv(
          enriched,
          paste(out_dir,"cluster",cid,"_",tool,"_enriched.csv",sep=""),
          row.names=FALSE)
        print(plotEnrich(enriched, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value"))
        plotEnrich(enriched, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
        ggsave(paste(out_dir,"cluster",cid,"_",tool,"_enriched_plot.png",sep=""))
      }
    }
  }
}
