library(dplyr)
library(monocle)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")

cds <- readRDS("deltxn_UMAP_2D_50clust.rds")


tmp.df <- pData(cds) %>%
  group_by(Cluster) %>%
  summarize(median_umap_1 = median(umap_1), median_umap_2 = median(umap_2))

ggplot(pData(cds), aes(x = umap_1, y = umap_2, color = Cluster)) +
    geom_point() +
    geom_text(aes(x = median_umap_1, y = median_umap_2, label = Cluster), data = tmp.df,
            size = 6, colour = "black") +
    ylim(13.8,15) +
    xlim(13.5,15) +
    theme(legend.position = "none") +
    monocle:::monocle_theme_opts()


# run DE test between 2 groups
dea_df <- two.set.differential.gene.test(new.cds, 
                      pData(cds)$Cluster == 11, 
                      pData(cds)$Cluster == 22, formal = TRUE)

# check expr
dea_df %>% filter(higher.expr == "Set 1") %>% 
  select(gene, symbol, set.1.mean.expr, set.2.mean.expr, higher.expr, log2.ratio) %>% head(20)

de_res <- differentialGeneTest(cds, relative_expr = FALSE, fullModelFormulaStr = ~Cluster, expressionFamily = "gaussianff")


new.cds = newCellDataSet(cds@assayData$exprs,
  phenoData = new("AnnotatedDataFrame", pData(cds)),
  featureData = new("AnnotatedDataFrame", fData(cds)),
                    expressionFamily = gaussianff())



### perform DEA between two sets of cells (from JP)
two.set.differential.gene.test = function(cds, set.1.filter, set.2.filter, formal = F, cores = 1, thresh = 1) {
  message(paste("# of cells in set 1:", sum(set.1.filter)))
  message(paste("# of cells in set 2:", sum(set.2.filter)))
  
  s1.cds = cds[, set.1.filter]
  s2.cds = cds[, set.2.filter]
  
  s1.norm.expr = Biobase::exprs(s1.cds)
  s2.norm.expr = Biobase::exprs(s2.cds)
  
  s1.mean = apply(s1.norm.expr, 1, mean)
  s2.mean = apply(s2.norm.expr, 1, mean)
  
  #s1.n.umi = apply(Biobase::exprs(s1.cds), 1, sum)
  #s2.n.umi = apply(Biobase::exprs(s2.cds), 1, sum)
  
  higher.expr = ifelse(s1.mean > s2.mean, "Set 1", "Set 2")
  
  s1.ratio = s1.mean / s2.mean
  s2.ratio = s2.mean / s1.mean
  log2.ratio = ifelse(
    s1.mean == 0 & s2.mean == 0, 0, ifelse(
      higher.expr == "Set 1", log2(s1.ratio), log2(s2.ratio)))
  
  s1.n.expr = apply(Biobase::exprs(s1.cds), 1, function(x) mean(x >= thresh)) #sum?
  s2.n.expr = apply(Biobase::exprs(s2.cds), 1, function(x) mean(x >= thresh)) #sum?
  
  s1.precision = s1.n.expr / (s1.n.expr + s2.n.expr)
  s1.recall = s1.n.expr / ncol(s1.cds)
  s2.precision = s2.n.expr / (s1.n.expr + s2.n.expr)
  s2.recall = s2.n.expr / ncol(s2.cds)
  
  precision = ifelse(higher.expr == "Set 1", s1.precision, s2.precision)
  recall = ifelse(higher.expr == "Set 1", s1.recall, s2.recall)
  
  f.score = 2 * precision * recall / (precision + recall)
  
  res = data.frame(
    gene = fData(cds)$id,
    symbol = fData(cds)$gene_short_name,
    set.1.mean.expr = s1.mean,
    set.2.mean.expr = s2.mean,
    higher.expr = higher.expr,
    log2.ratio = log2.ratio,
    precision = precision,
    recall = recall,
    f.score = f.score
  ) %>% arrange(-log2.ratio)
  
  
  if (formal) {
    pData(cds)$tmp = ifelse(set.1.filter, 1, ifelse(set.2.filter, 2, NA))
    
    cds.subset = cds[, set.1.filter | set.2.filter]
  
    message("Computing differential expression p-values")
    
    DEG = differentialGeneTest(cds.subset,
                        fullModelFormulaStr = "~ tmp", cores = cores, relative_expr = FALSE)
    
    print(head(DEG))
    
    res = inner_join(res,
                     DEG %>% select(gene = id,
                                    p.val = pval, q.val = qval),
                     by = "gene") %>% arrange(q.val)
    
    pData(cds)$tmp = NULL
  }
  
  return(res)
}







### perform DEA between two sets of cells (from JP)
two.set.differential.gene.test.short = function(cds, set.1.filter, set.2.filter, cores = 1, thresh = 1) {
  message(paste("# of cells in set 1:", sum(set.1.filter)))
  message(paste("# of cells in set 2:", sum(set.2.filter)))
  
  s1.cds = cds[, set.1.filter]
  s2.cds = cds[, set.2.filter]
  
  s1.norm.expr = Biobase::exprs(s1.cds)
  s2.norm.expr = Biobase::exprs(s2.cds)
  
  s1.mean = apply(s1.norm.expr, 1, mean)
  s2.mean = apply(s2.norm.expr, 1, mean)
  
  #s1.n.umi = apply(Biobase::exprs(s1.cds), 1, sum)
  #s2.n.umi = apply(Biobase::exprs(s2.cds), 1, sum)
  
  higher.expr = ifelse(s1.mean > s2.mean, "Set 1", "Set 2")
  
  s1.ratio = s1.mean / s2.mean
  s2.ratio = s2.mean / s1.mean
  log2.ratio = ifelse(
    s1.mean == 0 & s2.mean == 0, 0, ifelse(
      higher.expr == "Set 1", log2(s1.ratio), log2(s2.ratio)))
  
  s1.n.expr = apply(Biobase::exprs(s1.cds), 1, function(x) mean(x >= thresh)) #sum?
  s2.n.expr = apply(Biobase::exprs(s2.cds), 1, function(x) mean(x >= thresh)) #sum?
  
  s1.precision = s1.n.expr / (s1.n.expr + s2.n.expr)
  s1.recall = s1.n.expr / ncol(s1.cds)
  s2.precision = s2.n.expr / (s1.n.expr + s2.n.expr)
  s2.recall = s2.n.expr / ncol(s2.cds)
  
  precision = ifelse(higher.expr == "Set 1", s1.precision, s2.precision)
  recall = ifelse(higher.expr == "Set 1", s1.recall, s2.recall)
  
  f.score = 2 * precision * recall / (precision + recall)
  
  res = data.frame(
    gene = fData(cds)$id,
    symbol = fData(cds)$gene_short_name,
    set.1.mean.expr = s1.mean,
    set.2.mean.expr = s2.mean,
    higher.expr = higher.expr,
    log2.ratio = log2.ratio,
    precision = precision,
    recall = recall,
    f.score = f.score
  ) %>% arrange(-log2.ratio)
}

### EXAMPLE #####

##### note: run with formal = T to get q values (or FALSE if you want other info but faster)


# filter for expressed genes
expressed_genes <-  row.names(subset(fData(cds), num_cells_expressed >= 10))
filtered_cds <- cds[expressed_genes,]

# run DE test between 2 groups
dea_df <- two.set.differential.gene.test(cds, 
                                         pData(cds)$Cluster == 1, 
                                         pData(cds)$Cluster == 2,
                                         formal = T)

# name groups by what they are
dea_df$higher.cluster = rep("unknown", nrow(dea_df))
dea_df$higher.cluster <- with(dea_df, ifelse(higher.expr == "Set 1", "Cluster 8", higher.cluster))
dea_df$higher.cluster <- with(dea_df, ifelse(higher.expr == "Set 2", "Cluster 16", higher.cluster))

