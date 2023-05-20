library(venn)

# takes the path to the csv with aggregated markers and a directory path where to store images
plot_venn = function(markers, out_dir) {
  COLORS = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
  cluster_ids = unique(markers$cluster)
  for (i in cluster_ids){
    cluster_markers = markers[markers$cluster == i,]
    pdf(file=paste(out_dir,"cluster", i, "_markers_size_venn.pdf",sep=""))
    venn(
      split(cluster_markers$gene, cluster_markers$tool), 
      zcolor = COLORS,
      box = FALSE
      )
    dev.off()
  }
  pdf(file=paste(out_dir,"markers_size_venn.pdf",sep=""))
  venn(split(markers$gene, markers$tool), box = FALSE, zcolor = COLORS)
  dev.off()
}
