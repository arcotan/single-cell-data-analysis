library(venn)

# takes the path to the csv with aggregated markers and a directory path where to store images
plot_venn = function(markers, out_dir) {
  cluster_ids = unique(markers$cluster)
  for (i in cluster_ids){
    cluster_markers = markers[markers$cluster == i,]
    setEPS()
    postscript(paste(out_dir,"cluster", i, "_markers_size_venn.eps",sep=""))
    venn(
      split(cluster_markers$gene, cluster_markers$tool), 
      zcolor = "style",
      box = FALSE
      )
    dev.off()
    
    cluster_markers_top = cluster_markers[cluster_markers$rank <= 5,]
    setEPS()
    postscript(paste(out_dir,"cluster", i, "_markers_names_venn.eps",sep=""))
    venn(
      split(cluster_markers_top$gene, cluster_markers_top$tool),
      zcolor = "style",
      box = FALSE,
      ilabels = TRUE # TODO: does not work...
    )
    dev.off()
  }
  setEPS()
  postscript(paste(out_dir,"markers_size_venn.eps",sep=""))
  venn(split(markers$gene, markers$tool), box = FALSE)
  dev.off()
}
