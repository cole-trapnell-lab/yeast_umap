library(dplyr)
library(monocle)
library(ggplot2)
library(reshape2)

setwd("~/Dropbox (Cole Trapnell's Lab)/yeast_txnClusters/useful_files/")

cds <- readRDS("deltxn_UMAP_2D_50clust.rds")

# plot clusters with overlaid labels 
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

# modify data distribution label
cds = newCellDataSet(
  log2(exprs(cds)/100),
  phenoData = new("AnnotatedDataFrame", pData(cds)),
  featureData = new("AnnotatedDataFrame", fData(cds)),
  expressionFamily = gaussianff()
)

# run DE test between 2 groups
dea_df = two.set.differential.gene.test(cds, 
                    pData(cds)$Cluster == 11, 
                    pData(cds)$Cluster == 22, formal = T)

# check expr
dea_df %>% filter(higher.expr == "Set 1") %>% 
  select(gene, symbol, set.1.mean.expr, set.2.mean.expr, higher.expr, log2.ratio) %>% head(20)

# look at p-val distribution
hist(dea_df$p.val)

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
                     by = "gene") %>% arrange(q.val)
    
    pData(cds)$tmp = NULL
  }
  
  return(res)
}
