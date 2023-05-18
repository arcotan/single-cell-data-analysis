library(venn)

# takes the path to the csv with aggregated markers and a directory path where to store images
plot_venn = function(markers, out_dir) {
  cluster_ids = unique(markers$cluster)
  for (i in cluster_ids){
    cluster_markers = markers[markers$cluster == i,]
    #setEPS()
    #pdf(file=paste(out_dir,"cluster", i, "_markers_size_venn.pdf",sep=""))
    #postscript(paste(out_dir,"cluster", i, "_markers_size_venn.eps",sep=""))
    venn(
      split(cluster_markers$gene, cluster_markers$tool), 
      zcolor = "style",
      box = FALSE
      )
    dev.off()
  }
  # TODO
  #setEPS()
  #postscript(paste(out_dir,"markers_size_venn.eps",sep=""))
  pdf(file=paste(out_dir,"markers_size_venn.pdf",sep=""))
  venn(split(markers$gene, markers$tool), box = FALSE, zcolor = "style")
  dev.off()
}
