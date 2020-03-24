suppressPackageStartupMessages({library(monocle3)
  library(dplyr)
  library(viridis)
  library(reshape2)
  library(ggplot2)
  library(data.table)})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)

# this notebook has been updated to do the analysis with the newest version of monocle 3 (v0.2.0)

# Load data ---------------------------------------------------------------

cds = readRDS("R_objects/yeast_del_strain_monocle3_cds.RDS")

# Preprocess and run UMAP dimensionality reduction ------------------------

cds = detect_genes(cds)
cds = preprocess_cds(cds, num_dim = 100, norm_method = "log")
cds = reduce_dimension(cds,
                       preprocess_method = "PCA",
                       reduction_method = "UMAP", 
                       max_components = 2, 
                       umap.min_dist = 0.05, 
                       umap.n_neighbors = 10L,
                       umap.metric = "cosine")

cds = cluster_cells(cds, resolution = 1e-4, k = 3)

saveRDS(cds, "R_objects/yeast_del_strain_monocle3_cds.RDS")

# plot strain umap
plot_cells(cds, color_cells_by = "cluster_50", 
           cell_size = 0.75,
           label_groups_by_cluster = F,
           label_cell_groups = T, 
           group_label_size = 4)

# plot functional groups
plot_cells(cds, color_cells_by = "cluster", 
           label_cell_groups = F, 
           group_label_size = 4) +
  theme(legend.position = "none") +
  facet_wrap(~kemmeren_functional_group)

# plot expression of specific genes
plot_cells(cds, genes = "YCH1", 
           cell_size = 0.75, norm_method = "size_only") +
  scale_color_viridis_c()


# Visualize sub-clusters ---------------------------------------------------------

# extract large main cluster (clusters 1 - 10) and re-analyze
cds_sub = cds[,colData(cds)$cluster_50 == 3,]
cds_sub = detect_genes(cds_sub)
cds_sub = preprocess_cds(cds_sub, num_dim = 20, norm_method = "log")
cds_sub = reduce_dimension(cds_sub,
                       preprocess_method = "PCA",
                       reduction_method = "UMAP", 
                       max_components = 2, 
                       umap.min_dist = 0.05, 
                       umap.n_neighbors = 10L,
                       umap.metric = "cosine")

cds_sub = cluster_cells(cds_sub, resolution = 1e-4, k = 3)

# plot UMAP with sub-cluster designations from paper
plot_cells(cds_sub, color_cells_by = "sub_cluster", 
           cell_size = 1,
           label_groups_by_cluster = F,
           label_cell_groups = T, 
           group_label_size = 4)

# plot specific gene expression
plot_cells(cds_sub, genes = "ATP8", 
           cell_size = 1.25, norm_method = "size_only") +
  scale_color_viridis_c()




