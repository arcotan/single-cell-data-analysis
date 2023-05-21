library(venn)

# takes the path to the csv with aggregated markers and a directory path where to store images
plot_venn = function(markers, out_dir, dataset_name = '') {
  COLORS = c("#66c2a5","#e78ac3", "#8da0cb", "#a6d854", "#fc8d62")
  cluster_ids = unique(markers$cluster)
  for (i in cluster_ids){
    cluster_markers = markers[markers$cluster == i,]
    pdf(file=paste(out_dir,"cluster", i, "_markers_size_venn.pdf",sep=""))
    venn(
      split(cluster_markers$gene, cluster_markers$tool), 
      zcolor = COLORS,
      box = FALSE,
      )
    text(150, 1000, dataset_name, cex=1.5)
    text(46, 950, paste("Cluster", i), cex=1)
    dev.off()
  }
  pdf(file=paste(out_dir,"markers_size_venn.pdf",sep=""))
  venn(split(markers$gene, markers$tool), box = FALSE, zcolor = COLORS)
  text(150, 950, paste(dataset_name), cex=1.5)
  labs('tgfsrd')
  dev.off()
}
