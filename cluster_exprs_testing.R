library(dplyr)
library(monocle)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")

#cds <- readRDS("deltxn_UMAP_2D_50clust.rds")
cds <- readRDS("deltxn_UMAP_2D_171clust.rds")

# # plot clusters with overlaid labels 
# tmp.df <- pData(cds) %>%
#   group_by(Cluster) %>%
#   summarize(median_umap_1 = median(umap_1), median_umap_2 = median(umap_2))
# 
# ggplot(pData(cds), aes(x = umap_1, y = umap_2, color = Cluster)) +
#     geom_point() +
#     geom_text(aes(x = median_umap_1, y = median_umap_2, label = Cluster), data = tmp.df,
#             size = 6, colour = "black") +
#     ylim(13.8,15) +
#     xlim(13.5,15) +
#     theme(legend.position = "none") +
#     monocle:::monocle_theme_opts()
# 
# # modify data distribution label
# cds = newCellDataSet(
#   log2(exprs(cds)/100),
#   phenoData = new("AnnotatedDataFrame", pData(cds)),
#   featureData = new("AnnotatedDataFrame", fData(cds)),
#   expressionFamily = gaussianff()
# )
# 
# sub.dat <- read.table("mainclust1-10_subclust.txt", header = TRUE)
# sub.dat$sub.clust <- paste(sub.dat$main_cluster, sub.dat$Cluster, sep = "-")
# 
# # add column to pData with subcluster information
# tmp.pd <- merge(pData(cds), sub.dat %>% select(strain, sub.clust), by = "strain", all = TRUE)
# tmp.pd$sub.clust[is.na(tmp.pd$sub.clust)] <- as.character(tmp.pd$Cluster[is.na(tmp.pd$sub.clust)])
# rownames(tmp.pd) = tmp.pd$strain
# tmp.pd = tmp.pd[rownames(pData(cds)),]
# identical(tmp.pd$strain, pData(cds)$strain) #check
# pData(cds)$sub.cluster <- tmp.pd$sub.clust
# length(unique(pData(cds)$sub.cluster)) #there are 171

#saveRDS(cds, "deltxn_UMAP_2D_171clust.rds")

# run DE test between 2 groups
# you'll need to pull in the two.set.DE test function at the bottom of this script
dea_df = two.set.differential.gene.test(cds, 
                    pData(cds)$Cluster == 11, 
                    pData(cds)$Cluster == 22, formal = T)

# check expr
dea_df %>% filter(higher.expr == "Set 1") %>% 
  select(gene, symbol, set.1.mean.expr, set.2.mean.expr, higher.expr, log2.ratio) %>% head(20)

# assess lambda GC (genomic inflation factor)
dea_df$test.stat = qchisq(1.0 - dea_df$p.val, df = 1)

lambda.gc = median(dea_df$test.stat) / qchisq(0.5, 1)
lambda.gc

dea_df$inflation.corrected.p.val = pchisq(dea_df$test.stat / lambda.gc, df = 1, lower.tail = F)
dea_df$inflation.corrected.q.val = p.adjust(dea_df$inflation.corrected.p.val, method = "fdr")

# number of DEGs with and without correction
sum(dea_df$q.val < 0.05) / nrow(dea_df)
sum(dea_df$inflation.corrected.q.val < 0.05) / nrow(dea_df)

# look at top genes
dea_df %>% filter(dea_df$inflation.corrected.q.val < 0.05) %>% head(20)



## DEG test for all clusters against the rest of the strains to determine cluster-specific DEGs
cluster_indices = sort(as.character(unique(pData(cds)$sub.cluster)))

cluster_indices

# Instead of re-running below you can use the above cell...
DEG.results = lapply(cluster_indices,
    function(this.cluster) {
  message("Finding markers for cluster ", this.cluster)
    cbind(
        two.set.differential.gene.test(
            cds,
            pData(cds)$sub.cluster == this.cluster,
            pData(cds)$sub.cluster != this.cluster,
            formal = T,
            cores = 10),
        data.frame(cluster = this.cluster))
})

saveRDS(DEG.results, file = "~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/yeast_txn_all-subclust_DEGs.rds")

#DEG.results <- readRDS("yeast_txn_all-clust_DEGs.rds")
DEG.results <- readRDS("yeast_txn_all-subclust_DEGs.rds")

# this function combines DEG results from each cluster against all others 
cluster.markers.50 =
  do.call(rbind, lapply(DEG.results, head, n = 50L)) %>%
  select(Cluster = cluster, gene, symbol,
    p.val, q.val, everything()) %>%
  arrange(q.val)


cluster.markers.all =
  do.call(rbind, DEG.results) %>%
  select(Cluster = cluster, gene, symbol,
         p.val, q.val, everything()) %>%
  arrange(q.val)

write.csv(cluster.markers.50, 
          file = "~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/cluster_marker_top50.csv",
          row.names = FALSE)

cluster.markers.50 %>% filter(Cluster == "1-1") %>% head()

DEG %>% select(gene = id,
               p.val = pval, q.val = qval,
               set.1.mean.expr, 
               set.2.mean.expr, diff,
               inflation.corrected.p.val,
               inflation.corrected.q.val)

write.csv(cluster.markers, 
     file = "~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/cluster_marker_top50.csv",
    haeder = FALSE)

# label label complexes
cluster.markers$Complex = rep("unknown", nrow(cluster.markers))
cluster.markers$Complex = with(cluster.markers, ifelse(Cluster == 44, "Ribosomal", Complex))


cluster.degs <- cluster.markers %>% 
                filter(inflation.corrected.q.val < 0.001) %>% group_by(Cluster) %>% 
                arrange(inflation.corrected.q.val)


write.csv(cluster.degs, file = "~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/cluster_marker_degs.csv", row.names = FALSE)





#### wrapper function for differential gene expression testing (modified from J. Packer)
two.set.differential.gene.test = function(cds, set.1.filter, set.2.filter, formal = F, cores = 1, thresh = 1) {
  message(paste("# of cells in set 1:", sum(set.1.filter)))
  message(paste("# of cells in set 2:", sum(set.2.filter)))
  
  s1.cds = cds[, set.1.filter]
  s2.cds = cds[, set.2.filter]
  
  s1.mean = apply(Biobase::exprs(s1.cds), 1, mean)
  s2.mean = apply(Biobase::exprs(s2.cds), 1, mean)
  
  diff = s1.mean - s2.mean
  
  res = data.frame(
    gene = fData(cds)$id,
    symbol = fData(cds)$gene_short_name,
    set.1.mean.expr = s1.mean,
    set.2.mean.expr = s2.mean,
    diff = diff
  ) %>% arrange(-diff)
  
  if (formal) {
    pData(cds)$tmp = ifelse(set.1.filter, 1, ifelse(set.2.filter, 2, NA))
    
    cds.subset = cds[, set.1.filter | set.2.filter]
    
    message("Computing differential expression p-values")
    
    DEG = differentialGeneTest(cds.subset,
                               fullModelFormulaStr = "~ tmp", cores = cores, relative_expr = FALSE)
    
    DEG$test.stat = qchisq(1.0 - DEG$pval, df = 1)
    
    infl = median(qchisq(1.0 - DEG$pval, df = 1))
    
    lambda.gc = median(DEG$test.stat) / qchisq(0.5, 1)
    
    message(paste("inflation value is", infl))
    
    DEG$inflation.corrected.p.val = pchisq(DEG$test.stat / lambda.gc, df = 1, lower.tail = F)
    DEG$inflation.corrected.q.val = p.adjust(DEG$inflation.corrected.p.val, method = "fdr")
    
    res = inner_join(res,
                     DEG %>% select(gene = id,
                                    p.val = pval, q.val = qval, 
                                    inflation.corrected.p.val,
                                    inflation.corrected.q.val),
                     by = "gene") %>% arrange(inflation.corrected.q.val)
    
    pData(cds)$tmp = NULL
  }
  
  return(res)
}
