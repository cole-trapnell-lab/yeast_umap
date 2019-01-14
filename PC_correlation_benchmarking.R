library(dplyr)
library(tidyr)
library(stringr)
library(monocle)
library(ggplot2)
library(reshape2)
library(fields)

#setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")
setwd("~/Desktop/yeast_txnClusters/useful_files/")


cds <- readRDS("deltxn_UMAP_2D_171clust_V2.rds")

# generate a list of protein complexes associated with genes (from Benschop, et al 2010)
PCs <- read.csv("Kemmeren_Consensus.csv", stringsAsFactors = FALSE)
colnames(PCs)
PCs <- PCs %>% select(Complex = X..Complex.name, sys_gene = Complex.members..systematic.name.)

pc_list <- separate_rows(PCs, sys_gene, sep = ";")
pc_list$sys_gene <- str_trim(pc_list$sys_gene)
write.table(pc_list, "protein_complex_list.tsv", sep = "\t")

table(duplicated(pc_list$sys_gene))

pdat <- pData(cds)
length(intersect(pc_list$sys_gene, pdat$strain.sys))

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

cmplx_dist$gene_pair <- paste(cmplx_dist$gene1, cmplx_dist$gene2, sep = "_")

# filter for things in the same complex (UMAP)
pc_dists <- cmplx_dist %>% filter(gene1.complex == gene2.complex & gene1 != gene2) %>%
        select(gene1, gene2, dist, gene_pair)
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
nonpc_dists <- cmplx_dist %>% filter(!(gene_pair %in% pc_dists$gene_pair) & gene1 != gene2) %>%
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

# add a column of gene pairs
cmplx_cor$gene_pair <- paste(cmplx_cor$gene1, cmplx_cor$gene2, sep = "_")

# filter for things in the same complex
pc_cor <- cmplx_cor %>% filter(gene1.complex == gene2.complex & gene1 != gene2) %>%
  select(gene1, gene2, correlation, gene1.complex, gene_pair)
pc_cor <- pc_cor[!duplicated(pc_cor[,c("gene1", "gene2")]),]
#pc_cor <- pc_cor %>% group_by(gene1) %>% summarize(cor_mean = mean(correlation))
pc_cor$group <- "PC_cor"

head(pc_cor,10)
dim(pc_cor)
head(pc_cor)

identical(pc_dists$gene1, pc_cor$gene1)

df <- cbind(pc_dists, pc_cor %>% select(correlation))

ggplot(df, aes(dist, correlation)) +
    geom_point() +
    geom_density_2d() +
    xlim(0,0.02) +
    geom_smooth(method = "lm") +
    monocle:::monocle_theme_opts()

ggplot(df, aes(dist, correlation)) +
  geom_point() +
  geom_density_2d() +
  xlim(0,1) +
  geom_smooth(method = "loess") +
  monocle:::monocle_theme_opts()

summary(lm(correlation ~ dist, data = df))
summary(lm(correlation ~ dist, data = df2))

ggplot(pc_cor, aes(x = cor_mean)) + geom_density() + 
  monocle:::monocle_theme_opts()

# filter for things NOT in the same complex (correlation)
nonpc_cor <- cmplx_cor %>% filter(!(gene_pair %in% pc_cor$gene_pair) & gene1 != gene2) %>%
  select(gene1, gene2, correlation)
nonpc_cor <- nonpc_cor[!duplicated(nonpc_cor[,c("gene1", "gene2")]),]

nonpc_cor %>% filter(grepl("YNL330C", gene1)) %>% filter(grepl("YOL004W", gene2)) %>% head(10)

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


identical(pc_dists$gene1, pc_cor$gene1)

nonpc_all <- cbind(nonpc_dists, nonpc_cor %>% select(correlation))

ggplot(nonpc_all, aes(dist, correlation)) +
  geom_point(size=0.1) +
  geom_density_2d() +
  xlim(0,.05) +
  geom_smooth(method = "lm") +
  monocle:::monocle_theme_opts()

#dude = nonpc_all[dfnonpc_all2$dist > 0,]
#head(dude[order(-(dude$correlation + -dude$dist)),],20)

nonpc_all %>% filter(dist > 0) %>% arrange(-correlation, -dist)


# stats
wilcox.test(pc_cor$correlation, nonpc_cor_2k$correlation)




###### ROC ANALYSIS #######

cmplx_dist_filt <- cmplx_dist %>% filter(gene1 != gene2)
dim(cmplx_dist_filt)

cmplx_dist_filt$pair.in.complex <- ifelse(cmplx_dist_filt$gene1.complex == cmplx_dist_filt$gene2.complex, TRUE, FALSE)
head(cmplx_dist_filt)


cmplx_dist_filt$dist_bin = cut(cmplx_dist_filt$dist, seq(0,0.2,0.01))

cmplx_dist_filt %>% group_by(dist_bin, pair.in.complex) %>% summarize(count = n()) %>% head()

# add a column with cumulative bins
dude_true <- dude %>% filter(pair.in.complex == TRUE)
dude_true <- within(dude_true, acc_dist <- cumsum(count))

dude_false <- dude %>% filter(pair.in.complex == TRUE)
dude_false <- within(dude_true, acc_dist <- cumsum(count))

dude_final <- rbind(dude_true, dude_false)



# dist hist
hist(cmplx_dist_filt$dist, breaks=1000, ylim=c(0,25000), xlim=c(0,1))

#
tpr = dude[dude$pair.in.complex == FALSE,]$count / sum(dude[dude$pair.in.complex == FALSE,]$count)
fnr = dude[dude$pair.in.complex == TRUE,]$count / sum(dude[dude$pair.in.complex == TRUE,]$count)

plot(tpr, 1-fnr, type='l', xlim=c(0,1), ylim=c(0,1))

### one way to do it, maybe ###

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

tmp <- simple_roc(for_roc$pair.in.complex, for_roc$dist)

tmp %>% 
  ggplot(aes(x = TPR, y = 1-FPR)) +
      geom_point() +
      xlim(0,1) +
      ylim(0,1) +
      theme_classic()

### another way, with a different result ###

library(plotROC)

for_roc %>% 
  ggplot(aes(d = pair.in.complex, m = dist)) +
  geom_roc() +
  theme_bw()


### annnnnd a third way, with a third result ###
library(pROC)

plot(roc(for_roc$pair.in.complex, for_roc$dist), col="blue", lwd=3)
