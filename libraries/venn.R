library(venn)

# takes the path to the csv with aggregated markers and a directory path where to store images
plot_venn = function(markers, out_dir) {
  cluster_ids = unique(markers$cluster)
  for (i in cluster_ids){
    cluster_markers = markers[markers$cluster == i,]
    pdf(file=paste(out_dir,"cluster", i, "_markers_size_venn.pdf",sep=""))
    venn(
      split(cluster_markers$gene, cluster_markers$tool), 
      zcolor = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"),
      box = FALSE
      )
    dev.off()
  }
  pdf(file=paste(out_dir,"markers_size_venn.pdf",sep=""))
  venn(split(markers$gene, markers$tool), box = FALSE, zcolor = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"))
  dev.off()
}
