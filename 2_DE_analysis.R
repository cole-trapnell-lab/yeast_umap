suppressPackageStartupMessages({
  library(dplyr)
  library(viridis)
  library(reshape2)
  library(ggplot2)
  library(data.table)})

#library(monocle3)
library(devtools)
load_all("/Users/laurensaunders/Documents/Software/monocle3")

# This notebook is to do differential expression testing using the most updated version of monocle 3

# Load data ---------------------------------------------------------------

cds = readRDS("R_objects/yeast_del_strain_monocle3_cds.RDS")

#### wrapper function for differential gene expression testing (modified from J. Packer)
two.set.differential.expression = function(cds, set.1.filter, set.2.filter, formal = F, cores = 1, thresh = 1) {
  message(paste("# of cells in set 1:", sum(set.1.filter)))
  message(paste("# of cells in set 2:", sum(set.2.filter)))
  
  s1.cds = cds[, set.1.filter]
  s2.cds = cds[, set.2.filter]
  
  s1.mean = apply(assay(s1.cds), 1, mean)
  s2.mean = apply(assay(s2.cds), 1, mean)
  
  diff = s1.mean - s2.mean
  
  res = data.frame(
    gene = rowData(cds)$id,
    symbol = rowData(cds)$gene_short_name,
    set.1.mean.expr = s1.mean,
    set.2.mean.expr = s2.mean,
    diff = diff
  ) %>% arrange(-diff)
  
  if (formal) {
    colData(cds)$tmp = ifelse(set.1.filter, 1, ifelse(set.2.filter, 2, NA))
    
    cds.subset = cds[, set.1.filter | set.2.filter]
    
    message("Computing differential expression p-values")
    
    fits = fit_models(cds.subset, model_formula_str = "~tmp", cores = cores, expression_family = "gaussian")
    DEG = coefficient_table(fits)
    
    DEG$test.stat = qchisq(1.0 - DEG$p_value, df = 1)
    
    infl = median(qchisq(1.0 - DEG$p_value, df = 1))
    
    lambda.gc = median(DEG$test.stat) / qchisq(0.5, 1)
    
    message(paste("inflation value is", infl))
    
    DEG$inflation.corrected.p.val = pchisq(DEG$test.stat / lambda.gc, df = 1, lower.tail = F)
    DEG$inflation.corrected.q.val = p.adjust(DEG$inflation.corrected.p.val, method = "fdr")
    
    res = inner_join(res,
                     DEG %>% select(gene = id,
                                    p_value, q_value, 
                                    inflation.corrected.p.val,
                                    inflation.corrected.q.val),
                     by = "gene") %>% arrange(inflation.corrected.q.val)
    
    colData(cds)$tmp = NULL
  }
  
  return(res)
}

dea_df = two.set.differential.expression(cds, 
                                        colData(cds)$cluster_50 == 11, 
                                        colData(cds)$cluster_50 == 22, formal = T)

## DEG test for all clusters against the rest of the strains to determine cluster-specific DEGs
cluster_indices = sort(as.character(unique(colData(cds)$cluster_50)))

cluster_indices

# Instead of re-running below you can use the above cell...
DEG.results = lapply(cluster_indices,
                     function(this.cluster) {
                       message("Finding markers for cluster ", this.cluster)
                       cbind(
                         two.set.differential.gene.test.m3(
                           cds,
                           colData(cds)$cluster_50 == this.cluster,
                           colData(cds)$cluster_50 != this.cluster,
                           formal = T,
                           cores = 10),
                         data.frame(cluster = this.cluster))
                     })

DEG.results %>% filter(cluster_50 == 11) %>% arrange(q_value) 

degs = fread("~/Desktop/test_fits.txt")[,-1]
degs %>% filter(term != "(Intercept)") %>% arrange(q_value) %>% head

# check differences via plotting, may want to subset cds
plot_cells(cds, genes = "RPC34", 
           cell_size = 0.75, norm_method = "size_only") +
  scale_color_viridis_c()
