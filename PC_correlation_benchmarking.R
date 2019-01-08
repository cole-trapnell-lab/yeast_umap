library(dplyr)
library(tidyr)
library(stringr)
library(monocle)
library(ggplot2)
library(reshape2)
library(fields)

setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")

cds <- readRDS("deltxn_UMAP_2D_171clust_V2.rds")

# generate a list of protein complexes associated with genes (from Benschop, et al 2010)
PCs <- read.csv("Kemmeren_Consensus.csv", stringsAsFactors = FALSE)
colnames(PCs)
PCs <- PCs %>% select(Complex = X..Complex.name, sys_gene = Complex.members..systematic.name.)

pc_list <- separate_rows(PCs, sys_gene, sep = ";")
pc_list$sys_gene <- str_trim(pc_list$sys_gene)

table(duplicated(pc_list$sys_gene))
length(intersect(pc_list$sys_gene, pdat$strain.sys))


pdat <- pData(cds)

# filter pData for genes in complexes
pdat_filt <- pdat[pdat$strain.sys %in% intersect(pc_list$sys_gene, pdat$strain.sys),]

# make a matrix for distance calculations
umap_mat <- pdat_filt %>% select(strain.sys, umap_1, umap_2)
###umap_mat[umap_mat$strain.sys == "YAL056W",1][1] <- "YAL056W-2" #remove duplicated name
rownames(umap_mat) <- umap_mat$strain.sys
umap_mat <- Matrix::as.matrix(umap_mat[,2:3])
umap_mat[1:5,]

# generate pairwise distance matrix
udist <- rdist(umap_mat)

colnames(udist) <- rownames(umap_mat)
rownames(udist) <- rownames(umap_mat)



udist_melt <- melt(udist) %>% select(sys_gene = Var1, gene2 = Var2, dist = value)
head(udist_melt)
dim(udist_melt)

cmplx_dist <- left_join(udist_melt, pc_list, by = "sys_gene") %>% 
              select(gene1 = sys_gene, gene2, gene1.complex = Complex, dist)
cmplx_dist <- left_join(cmplx_dist, pc_list, by = c("gene2" = "sys_gene")) %>% 
              select(gene1, gene2, gene1.complex, gene2.complex = Complex, dist)

# filter for things in the same complex (UMAP)
pc_dists <- cmplx_dist %>% filter(gene1.complex == gene2.complex & gene1 != gene2) %>%
        select(gene1, gene2, dist)
pc_dists <- pc_dists[!duplicated(pc_dists[,c("gene1", "gene2")]),]
#pc_dists <- pc_dists %>% group_by(gene1) %>% summarize(dist_mean = mean(dist))
pc_dists$group <- "PC_umap" # add group name for combining lists to plot

head(pc_dists,40)
dim(pc_dists)
head(pc_dists)

ggplot(pc_dists, aes(x = dist)) + geom_density() + 
  monocle:::monocle_theme_opts() +
  xlim(0,2)

# filter for things NOT in the same complex (UMAP)
nonpc_dists <- cmplx_dist %>% filter(gene1.complex != gene2.complex & gene1 != gene2) %>%
  select(gene1, gene2, dist)
nonpc_dists <- nonpc_dists[!duplicated(nonpc_dists[,c("gene1", "gene2")]),]
#nonpc_dists <- nonpc_dists %>% group_by(gene1) %>% summarize(dist_mean = mean(dist))

head(nonpc_dists,40)
dim(nonpc_dists)
head(nonpc_dists)

#sample the same number of genes from the nonPC list because the numbers are really different
nonpc_2k <- nonpc_dists[sample(nrow(nonpc_dists), 2200),]

nonpc_2k$group <- "nonPC_umap"

comb_df <- rbind(pc_dists, nonpc_2k)
comb_df$neg_dist <- -(comb_df$dist)

ggplot(comb_df, aes(x = neg_dist, color = group)) + 
  geom_density(size = 1.2) + 
  scale_color_manual(values = c("#7c7c7c", "#ea3327")) +
  monocle:::monocle_theme_opts() +
  xlim(-1.5,0)

# stats
wilcox.test(pc_dists$dist, nonpc_2k$dist)


####### Compare to published correlation metric ##########

expr_mat <- Biobase::exprs(cds)
expr_mat[1:5, 1:5]
identical(pData(cds)$strain, colnames(Biobase::exprs(cds)))
colnames(expr_mat) <- pData(cds)$strain.sys
pc_expr <- expr_mat[,pdat_filt$strain.sys] #filter for strains where genes are from protein complexes
dim(pc_expr) # should be 619 cols

cor_mat <- cor(pc_expr)
cor_mat[1:5, 1:5]

# assess gene by gene correlation for those within complexes or not
cor_melt <- melt(cor_mat) %>% select(sys_gene = Var1, gene2 = Var2, correlation = value)
head(cor_melt)
dim(cor_melt)

cmplx_cor <- left_join(cor_melt, pc_list, by = "sys_gene") %>% 
  select(gene1 = sys_gene, gene2, gene1.complex = Complex, correlation)
cmplx_cor <- left_join(cmplx_cor, pc_list, by = c("gene2" = "sys_gene")) %>% 
  select(gene1, gene2, gene1.complex, gene2.complex = Complex, correlation)

# filter for things in the same complex
pc_cor <- cmplx_cor %>% filter(gene1.complex == gene2.complex & gene1 != gene2) %>%
  select(gene1, gene2, correlation)
pc_cor <- pc_cor[!duplicated(pc_cor[,c("gene1", "gene2")]),]
#pc_cor <- pc_cor %>% group_by(gene1) %>% summarize(dist_mean = mean(dist))
pc_cor$group <- "PC_cor"

head(pc_cor,10)
dim(pc_cor)
head(pc_cor)

ggplot(pc_cor, aes(x = correlation)) + geom_density() + 
  monocle:::monocle_theme_opts()

# filter for things NOT in the same complex (correlation)
nonpc_cor <- cmplx_cor %>% filter(gene1.complex != gene2.complex & gene1 != gene2) %>%
  select(gene1, gene2, correlation)
nonpc_cor <- nonpc_cor[!duplicated(nonpc_cor[,c("gene1", "gene2")]),]
#nonpc_cor <- nonpc_cor %>% group_by(gene1) %>% summarize(dist_mean = mean(dist))

head(nonpc_cor,10)
dim(nonpc_cor)
head(nonpc_cor)

#sample the same number of genes from the nonPC list because the numbers are really different
nonpc_cor_2k <- nonpc_cor[sample(nrow(nonpc_cor), 2200),]

nonpc_cor_2k$group <- "nonPC_cor"

comb_cor <- rbind(pc_cor, nonpc_cor_2k)

ggplot(comb_cor, aes(x = correlation, color = group)) + 
  geom_density(size = 1.2) + 
  scale_color_manual(values = c("#7c7c7c", "#ea3327")) +
  monocle:::monocle_theme_opts()

# stats
wilcox.test(pc_cor$correlation, nonpc_cor_2k$correlation)


# #combine all data sets (may not be ok to do this because of differences between cor and dist metrics)
# comb_df2 <- rbind(pc_dists, nonpc_2k, pc_cor, nonpc_cor_2k)
# comb_df2$neg_dist <- -(comb_df2$dist)
# 
# my_colors <- c()
# 
# ggplot(comb_df2, aes(x = neg_dist, color = group)) + 
#   geom_density(size = 1.2) + 
#   #scale_color_manual(values = c("#7c7c7c", "#ea3327")) +
#   monocle:::monocle_theme_opts() +
#   xlim(-1.5,0)





