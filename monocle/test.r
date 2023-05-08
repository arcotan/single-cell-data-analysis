rm(list = ls()) # clear the environment 
#load all the necessary libraries 
options(warn=-1) # turn off warning message globally 

# install.packages("flexclust")
# install.packages("mcclust")

suppressMessages(library(reticulate))

# py_install("r-reticulate", "louvain")

suppressMessages(library(devtools))
suppressMessages(load_all("/Users/andyshehmohajeri/Documents/monocle-dev-7-19/monocle-dev"))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))

library(reticulate)
import("louvain")

MCA <- readRDS("/Users/andyshehmohajeri/Downloads/MCA_merged_mat.rds")
cell.meta.data <- read.csv("/Users/andyshehmohajeri/Downloads/MCA_All-batch-removed-assignments.csv", row.names = 1)

overlapping_cells <- intersect(row.names(cell.meta.data), colnames(MCA))
gene_ann <- data.frame(gene_short_name = row.names(MCA), row.names = row.names(MCA))

pd <- new("AnnotatedDataFrame",data=cell.meta.data[overlapping_cells, ])
fd <- new("AnnotatedDataFrame",data=gene_ann)

MCA_cds <- newCellDataSet(MCA[, overlapping_cells], phenoData = pd,featureData =fd,
                           expressionFamily = negbinomial.size(),
                           lowerDetectionLimit=1)


# MCA_cds <- MCA_cds[, Matrix::colSums(exprs(MCA_cds)) > 1000 & Matrix::colSums(exprs(MCA_cds) > 0) > 500]
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
MCA_cds <- estimateSizeFactors(MCA_cds)
MCA_cds <- estimateDispersions(MCA_cds)


library(dplyr)

disp_table = dispersionTable(MCA_cds)
disp_table = disp_table %>% mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>% arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 2500)$gene_id)

MCA_cds = setOrderingFilter(MCA_cds, top_subset_genes)
MCA_cds <- preprocessCDS(MCA_cds,  method = 'PCA', 
                         norm_method = 'log', 
                         num_dim = 50, 
                         verbose = T)
MCA_cds <- reduceDimension(MCA_cds, max_components = 2,   
                       reduction_method = 'UMAP',
                       metric="correlation", 
                       min_dist = 0.75, 
                       n_neighbors = 50, 
                       verbose = T)

MCA_cds <- clusterCells(MCA_cds, method = 'louvain', res = 1e-6, louvain_iter = 1, verbose = T)

col_vector_origin <- c("#db83da",
                "#53c35d",
                "#a546bb",
                "#83b837",
                "#a469e6",
                "#babb3d",
                "#4f66dc",
                "#e68821",
                "#718fe8",
                "#d6ac3e",
                "#7957b4",
                "#468e36",
                "#d347ae",
                "#5dbf8c",
                "#e53e76",
                "#42c9b8",
                "#dd454a",
                "#3bbac6",
                "#d5542c",
                "#59aadc",
                "#cf8b36",
                "#4a61b0",
                "#8b8927",
                "#a24e99",
                "#9cb36a",
                "#ca3e87",
                "#36815b",
                "#b23c4e",
                "#5c702c",
                "#b79add",
                "#a55620",
                "#5076af",
                "#e38f67",
                "#85609c",
                "#caa569",
                "#9b466c",
                "#88692c",
                "#dd81a9",
                "#a35545",
                "#e08083",
                "#17becf",
                "#9edae5")
col_vector <- col_vector_origin[1:length(unique(as.character(pData(MCA_cds)$Tissue)))]
names(col_vector) <- unique(as.character(pData(MCA_cds)$Tissue))
options(repr.plot.width = 11)
options(repr.plot.height = 8)
plot_cell_clusters(MCA_cds, color_by = 'Tissue', cell_size = 0.1, show_group_id = T)  + 
theme(legend.text=element_text(size=6)) + #set the size of the text 
theme(legend.position="right") #put the color legend on the right


options(repr.plot.width = 11)
options(repr.plot.height = 8)
cluster_col_vector <- col_vector_origin[1:length(unique(as.character(pData(MCA_cds)$Cluster)))]
names(cluster_col_vector) <- unique(as.character(pData(MCA_cds)$Cluster))
plot_cell_clusters(MCA_cds, color_by = 'Cluster', cell_size = 0.1, show_group_id = T) + 
scale_color_manual(values = cluster_col_vector) + 
theme(legend.text=element_text(size=6)) + #set the size of the text 
theme(legend.position="right") #put the color legend on the right

Cluster_tissue_stat <- table(pData(MCA_cds)[, c('Cluster', 'Tissue')])

Cluster_tissue_stat <- apply(Cluster_tissue_stat, 1, function(x) x / sum(x))

Cluster_tissue_stat_ordered <- t(Cluster_tissue_stat)

options(repr.plot.width=10, repr.plot.height=7)
    
order_mat <- t(apply(Cluster_tissue_stat_ordered, 1, order))
max_ind_vec <- c()
    
for(i in 1:nrow(order_mat)) {
  tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
  max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
}
    
max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]

max_ind_vec <- c(max_ind_vec, setdiff(1:ncol(Cluster_tissue_stat), max_ind_vec))
Cluster_tissue_stat_ordered <- Cluster_tissue_stat_ordered[, row.names(Cluster_tissue_stat)[max_ind_vec]]
                             
pheatmap::pheatmap(Cluster_tissue_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10)) # annotation_row = annotation_row, annotation_colors = annotation_colors


Cluster_tissue_stat <- table(pData(MCA_cds)[, c('Cluster', 'Tissue')])

Cluster_tissue_stat <- apply(Cluster_tissue_stat, 2, function(x) x / sum(x))

Cluster_tissue_stat_ordered <- t(Cluster_tissue_stat)

options(repr.plot.width=10, repr.plot.height=7)
order_mat <- t(apply(Cluster_tissue_stat_ordered, 2, order))
max_ind_vec <- c()
    
for(i in 1:nrow(order_mat)) {
  tmp <- max(which(!(order_mat[i, ] %in% max_ind_vec)))
  max_ind_vec <- c(max_ind_vec, order_mat[i, tmp])
}

max_ind_vec <- max_ind_vec[!is.na(max_ind_vec)]

max_ind_vec <- c(max_ind_vec, setdiff(1:ncol(Cluster_tissue_stat), max_ind_vec))
Cluster_tissue_stat_ordered <- Cluster_tissue_stat_ordered[colnames(Cluster_tissue_stat)[max_ind_vec], ]
                             
pheatmap::pheatmap(Cluster_tissue_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10)) # annotation_row = annotation_row, annotation_colors = annotation_colors

options(repr.plot.width = 11)
options(repr.plot.height = 8)
MCA_cds_res_resolution <- clusterCells(MCA_cds, use_pca = F, method = 'louvain', res = 5e-5, louvain_iter = 1, verbose = T)
cluster_col_vector <- col_vector_origin[1:length(unique(as.character(pData(MCA_cds_res_resolution)$Cluster)))]
names(cluster_col_vector) <- unique(as.character(pData(MCA_cds_res_resolution)$Cluster))

plot_cell_clusters(MCA_cds_res_resolution, color_by = 'Cluster', cell_size = 0.1, show_group_id = T) + 
theme(legend.text=element_text(size=6)) + #set the size of the text 
theme(legend.position="right") #put the color legend on the right

options(repr.plot.width = 11)
options(repr.plot.height = 8)
MCA_cds_res_resolution <- clusterCells(MCA_cds, use_pca = F, method = 'louvain', res = 5e-7, louvain_iter = 1, verbose = T)
cluster_col_vector <- col_vector_origin[1:length(unique(as.character(pData(MCA_cds_res_resolution)$Cluster)))]
names(cluster_col_vector) <- unique(as.character(pData(MCA_cds_res_resolution)$Cluster))

plot_cell_clusters(MCA_cds_res_resolution, color_by = 'Cluster', cell_size = 0.1, show_group_id = T) + 
theme(legend.text=element_text(size=6)) + #set the size of the text 
theme(legend.position="right") #put the color legend on the right

start <- Sys.time()
spatial_res <- principalGraphTest(MCA_cds, relative_expr = TRUE, k = 25, cores = detectCores() - 2, verbose = FALSE)
end <- Sys.time()
end - start


cluster_marker_res <- find_cluster_markers(MCA_cds, spatial_res, group_by = 'Cluster', morans_I_threshold = 0.25)

genes <- (cluster_marker_res %>% dplyr::filter(mean > 0.5, percentage > 0.1) %>% dplyr::group_by(Group) %>% dplyr::slice(which.max(specificity)))
options(repr.plot.width=22, repr.plot.height=12)
plot_markers_by_group(MCA_cds, genes$gene_short_name, group_by = 'Cluster', ordering_type = 'maximal_on_diag')

genes <- (cluster_marker_res %>% dplyr::filter(mean > 0.5, percentage > 0.1) %>% dplyr::group_by(Group) %>% dplyr::top_n(5, wt = specificity))
plot_markers_cluster(MCA_cds, as.character(genes$gene_short_name), minimal_cluster_fraction = 0.05)

options(repr.plot.width=22, repr.plot.height=12)
genes <- cluster_marker_res %>% dplyr::filter(mean > 0.5, percentage > 0.1, specificity > 0.7) %>% dplyr::group_by(Group) %>% dplyr::arrange(Group, dplyr::desc(specificity))
plot_markers_cluster(MCA_cds, as.character(genes$gene_short_name), minimal_cluster_fraction = 0.05, show_rownames = F)