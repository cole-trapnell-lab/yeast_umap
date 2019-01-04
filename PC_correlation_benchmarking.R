library(dplyr)
library(tidyr)
library(monocle)
library(ggplot2)
library(reshape2)
library(viridis)
library(fields)

setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")

cds <- readRDS("deltxn_UMAP_2D_171clust.rds")

# generate a list of protein complexes associated with genes (from Benschop, et al 2010)
PCs <- read.csv("Kemmeren_Consensus.csv", stringsAsFactors = FALSE)
colnames(PCs)
PCs <- PCs %>% select(Complex = X..Complex.name, sys_gene = Complex.members..systematic.name.)

pc_list <- separate_rows(PCs, sys_gene, sep = ";")


table(duplicated(pc_list$sys_gene))

pdat <- pData(cds)

# make a matrix for distance calculations
umap_mat <- pdat %>% select(strain.sys, umap_1, umap_2)
#umap_mat[umap_mat$strain.sys == "YAL056W",1][1] <- "YAL056W-2" #remove duplicated name
rownames(umap_mat) <- umap_mat$strain.sys
umap_mat <- Matrix::as.matrix(umap_mat[,2:3])
umap_mat[1:5,]

# generate pairwise distance matrix
udist <- rdist(umap_mat)


# maybe something like this will work once we have things in the right format? 
umap_mat %>% group_by(Complex) %>%
  mutate(Dist = colMeans(as.matrix(dist(cbind(x, y)))))
