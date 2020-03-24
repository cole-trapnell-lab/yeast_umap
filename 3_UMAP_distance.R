library(dplyr)
library(tidyr)
library(stringr)
library(monocle)
library(ggplot2)
library(reshape2)
library(fields)
library(PRROC)
library(pROC)


setwd("~/Dropbox/manuscripts/yeast_txnClusters/yeast_umap/")


# Load in protein complex IP-MS data --------------------------------------

# generates a vector of protein complexes associated with genes (from Benschop, et al 2010)
ppi.ipms <- read.csv("Kemmeren_Consensus.csv", stringsAsFactors = FALSE)
ppi.ipms <- ppi.ipms %>% dplyr::select(Complex = complex_name, sys_gene = complex_members_sys)
ppi_list <- separate_rows(ppi.ipms, sys_gene, sep = ";")
ppi_list$sys_gene <- str_trim(ppi_list$sys_gene)


# Load in pairwise distances for each metric  -----------------------------------

dist.pca = read.table('pca_dist.txt')
dist.hd = read.table('hd_dist.txt')
dist.pcc = read.table('pcc_dist.txt')

# Calculate pairwise distances in UMAP space for downstream analyses
# Read in meta-data of umap coordinates
meta_umap = read.table('umap_metadata_clust50.txt')

# Filter for genes in the ppi dataset
pdat_filt <- meta_umap[meta_umap$strain.sys %in% intersect(ppi_list$sys_gene, meta_umap$strain.sys),]

# Calculate pair-wise distances in UMAP space; remove distal members
umap_mat <- pdat_filt %>% dplyr::select(umap_1, umap_2)
rownames(umap_mat) <- pdat_filt$strain.sys
umap_mat <- as.matrix(umap_mat)
udist <- as.matrix(dist(umap_mat, method='euclidean'))
udist_melt = melt(udist)
colnames(udist) <- rownames(umap_mat)
rownames(udist) <- rownames(umap_mat)
udist_melt <- melt(udist) %>% dplyr::select(sys_gene = Var1, gene2 = Var2, dist = value)
udist_melt <- udist_melt[udist_melt$dist <= 0.2,]

# add column showing whether gene pair occurs in the same complex
dist.umap <- left_join(udist_melt, ppi_list, by = "sys_gene") %>% 
  dplyr::select(gene1 = sys_gene, gene2, gene1.complex = Complex, dist)
dist.umap <- left_join(dist.umap, ppi_list, by = c("gene2" = "sys_gene")) %>% 
  dplyr::select(gene1, gene2, gene1.complex, gene2.complex = Complex, dist)
dist.umap$gene_pair <- paste(dist.umap$gene1, dist.umap$gene2, sep = "_")
dist.umap.ppi <- dist.umap %>% filter(gene1 != gene2)
dist.umap.ppi$pair.in.complex <- ifelse(dist.umap.ppi$gene1.complex == dist.umap.ppi$gene2.complex, TRUE, FALSE)


# Generate tables that identify gene pairs that encode members of  --------

# PCC
dist.pcc.filt <- dist.pcc[dist.pcc$sys_gene %in% intersect(ppi_list$sys_gene, dist.pcc$sys_gene),]

dist.pcc.ppi <- left_join(dist.pcc.filt, ppi_list, by = "sys_gene") %>% 
  dplyr::select(gene1 = sys_gene, gene2, gene1.complex = Complex, dist)
dist.pcc.ppi <- left_join(dist.pcc.ppi, ppi_list, by = c("gene2" = "sys_gene")) %>% 
  dplyr::select(gene1, gene2, gene1.complex, gene2.complex = Complex, dist)
dist.pcc.ppi$gene_pair <- paste(dist.pcc.ppi$gene1, dist.pcc.ppi$gene2, sep = "_")
dist.pcc.ppi <- dist.pcc.ppi %>% filter(gene1 != gene2)
dist.pcc.ppi$pair.in.complex <- ifelse(dist.pcc.ppi$gene1.complex == dist.pcc.ppi$gene2.complex, TRUE, FALSE)


# PCA
dist.pca.filt <- dist.pca[dist.pca$sys_gene %in% intersect(ppi_list$sys_gene, dist.pca$sys_gene),]

dist.pca.ppi <- left_join(dist.pca.filt, ppi_list, by = "sys_gene") %>% 
  dplyr::select(gene1 = sys_gene, gene2, gene1.complex = Complex, dist)
dist.pca.ppi <- left_join(dist.pca.ppi, ppi_list, by = c("gene2" = "sys_gene")) %>% 
  dplyr::select(gene1, gene2, gene1.complex, gene2.complex = Complex, dist)
dist.pca.ppi$gene_pair <- paste(dist.pca.ppi$gene1, dist.pca.ppi$gene2, sep = "_")
dist.pca.ppi <- dist.pca.ppi %>% filter(gene1 != gene2)
dist.pca.ppi$pair.in.complex <- ifelse(dist.pca.ppi$gene1.complex == dist.pca.ppi$gene2.complex, TRUE, FALSE)

# High-dimensional distance
dist.hd.filt <- dist.hd[dist.hd$sys_gene %in% intersect(ppi_list$sys_gene, dist.hd$sys_gene),]

dist.hd.ppi <- left_join(dist.hd.filt, ppi_list, by = "sys_gene") %>% 
  dplyr::select(gene1 = sys_gene, gene2, gene1.complex = Complex, dist)
dist.hd.ppi <- left_join(dist.hd.ppi, ppi_list, by = c("gene2" = "sys_gene")) %>% 
  dplyr::select(gene1, gene2, gene1.complex, gene2.complex = Complex, dist)
dist.hd.ppi$gene_pair <- paste(dist.hd.ppi$gene1, dist.hd.ppi$gene2, sep = "_")
dist.hd.ppi <- dist.hd.ppi %>% filter(gene1 != gene2)
dist.hd.ppi$pair.in.complex <- ifelse(dist.hd.ppi$gene1.complex == dist.hd.ppi$gene2.complex, TRUE, FALSE)

# ROC analysis ------------------------------------------------------------

# Calculate AUC
roc(dist.pcc.ppi$pair.in.complex, dist.pcc.ppi$dist, auc=TRUE)
roc(dist.pca.ppi$pair.in.complex, dist.pca.ppi$dist, auc=TRUE)
roc(dist.hd.ppi$pair.in.complex, dist.hd.ppi$dist, auc=TRUE)
roc(dist.umap.ppi$pair.in.complex, dist.umap.ppi$dist, auc=TRUE)

plot(roc(dist.pcc.ppi$pair.in.complex, dist.pcc.ppi$dist), col="darkgrey", lwd=3)
lines(roc(dist.pca.ppi$pair.in.complex, dist.pca.ppi$dist), col="black", lwd=3)
lines(roc(dist.hd.ppi$pair.in.complex, dist.hd.ppi$dist), col="grey", lwd=3)
lines(roc(dist.umap.ppi$pair.in.complex, dist.umap.ppi$dist), col="forestgreen", lwd=3)


