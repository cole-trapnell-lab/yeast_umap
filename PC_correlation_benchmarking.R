library(dplyr)
library(tidyr)
library(monocle)
library(ggplot2)
library(reshape2)
library(viridis)
library(fields)
library(stringr)

setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")

cds <- readRDS("deltxn_UMAP_2D_171clust.rds")

# generate a list of protein complexes associated with genes (from Benschop, et al 2010)
PCs <- read.csv("Kemmeren_Consensus.csv", stringsAsFactors = FALSE)
colnames(PCs)
PCs <- PCs %>% select(Complex = X..Complex.name, sys_gene = Complex.members..systematic.name.)

pc_list <- separate_rows(PCs, sys_gene, sep = ";")
pc_list$sys_gene <- str_trim(pc_list$sys_gene)

table(duplicated(pc_list$sys_gene))
length(intersect(pc_list$sys_gene, pdat$strain.sys))


pdat <- pData(cds)

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

cmplx_dist <- left_join(udist_melt, pc_list, by = "sys_gene") %>% select(gene1 = sys_gene, gene2, gene1.complex = Complex, dist)
cmplx_dist <- left_join(cmplx_dist, pc_list, by = c("gene2" = "sys_gene")) %>% select(gene1, gene2, gene1.complex, gene2.complex = Complex, dist)

# filter for things in the same complex
pc_dists <- cmplx_dist %>% filter(gene1.complex == gene2.complex & gene1 != gene2) %>%
        select(gene1, gene2, dist)
pc_dists <- pc_dists[!duplicated(pc_dists[,c("gene1", "gene2")]),]
#pc_dists <- pc_dists %>% group_by(gene1) %>% summarize(dist_mean = mean(dist))

head(pc_dists,40)
dim(pc_dists)
head(pc_dists)

pc_dists$group <- "PCs"

ggplot(pc_dists, aes(x = dist)) + geom_density() + 
  monocle:::monocle_theme_opts() +
  xlim(0,2)

# filter for things NOT in the same complex
nonpc_dists <- cmplx_dist %>% filter(gene1.complex != gene2.complex & gene1 != gene2) %>%
  select(gene1, gene2, dist)
nonpc_dists <- nonpc_dists[!duplicated(nonpc_dists[,c("gene1", "gene2")]),]
#nonpc_dists <- nonpc_dists %>% group_by(gene1) %>% summarize(dist_mean = mean(dist))

head(nonpc_dists,40)
dim(nonpc_dists)
head(nonpc_dists)

tmp3 <- nonpc_dists[sample(nrow(nonpc_dists), 2200),]

nonpc_dists$group <- "nonPCs"

comb_df <- rbind(pc_dists, nonpc_dists)

ggplot(comb_df, aes(x = dist, color = group)) + 
  geom_density(size = 1.2) + 
  scale_color_manual(values = c("#7c7c7c", "#ea3327")) +
  monocle:::monocle_theme_opts() +
  xlim(0,1.5)

wilcox.test(tmp$dist, tmp3$dist)

# maybe something like this will work once we have things in the right format? 
cmplx_dist %>% group_by(Complex) %>% mutate(Dist = colMeans())
