setwd("~/Dropbox/R")
source("ignition.txt")
options(stringsAsFactors = FALSE)
library(RColorBrewer)
library(tidyverse)

# install the SGD gene mappings
source("https://bioconductor.org/biocLite.R")
library(org.Sc.sgd.db)
library(GO.db)
library(topGO)
library(GOstats)

setwd('~/Dropbox/q_lab/projects/35_yeastTranscriptome_pathways/180913_monocle_lauren/')
data = read.csv('deltxn_UMAP_coords_50clust.csv')
strain.names = read.csv('../181003_allStrains_ordered_wSystematic.csv', header=FALSE)
#data$strain = strain.names$V2

setwd('~/Dropbox/q_lab/projects/35_yeastTranscriptome_pathways/180913_monocle_lauren/')
data.main = read.csv('deltxn_UMAP_coords_50clust.csv')
data.sub = read.table('../181001_subClustering/181003_50clust_1-10subclust.txt', header=TRUE)
data.sub$Cluster = paste(data.sub$main_cluster, data.sub$Cluster, sep='-')
data.sub = data.sub[,!names(data.sub) %in% c('num_genes_expressed','main_cluster')]

# replace rows of clusters 1-10 with the subcluster rows
data.all = data.main[data.main$Cluster > 10,]
data.all = rbind(data.all, data.sub)
data.all = data.all[match(data$strain, data.all$strain),]
data.all$strain = strain.names$V2
data.all$strain.common = strain.names$V1
head(data.all)
data = data.all


#for reference's sake
#how to get GO terms from an Entrez gene id following the manual
#http://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf

#org.Sc.sgdALIAS is an R object that provides
#mappings between entrez gene identifers and the GO
#identifers that they are directly associated with


x <- org.Sc.sgdALIAS
# Get the probe identifiers that are mapped to alias names
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])


#as the universal list, I will use all the genes with GO terms

universe <- mapped_genes
gene_subset = data[data$Cluster == '1-2',]$strain
length(gene_subset)
length(universe)


params <- new('GOHyperGParams',
              geneIds=gene_subset,
              universeGeneIds=universe,
              ontology='BP',
              pvalueCutoff=0.001,
              conditional=F,
              testDirection='over',
              annotation="org.Sc.sgd.db"
)
hgOver <- hyperGTest(params)
hgOver

result <- summary(hgOver)
head(result,20)


strain.names[strain.names$V2 %in% gene_subset,]



# write a function that just cycles through all the subclusters, then prints the results into a table
unique(data$Cluster)[order(unique(data$Cluster))]
#cluster_num = paste(1, seq(1,35), sep='-')
#cluster_num = paste(2, seq(1,24), sep='-')
cluster_num = paste(3, seq(1,19), sep='-')

goRunner = function(cluster_num) {
  gene_subset = data[data$Cluster == cluster_num,]$strain
  params <- new('GOHyperGParams',
                geneIds=gene_subset,
                universeGeneIds=universe,
                ontology='BP',
                pvalueCutoff=0.001,
                conditional=F,
                testDirection='over',
                annotation="org.Sc.sgd.db"
  )
  hgOver <- hyperGTest(params)
  result <- summary(hgOver)
  result$subclust = rep(cluster_num,length(result$Term))
  return(result)
}


goList = lapply(cluster_num, goRunner)
df.GO = do.call(rbind.data.frame, goList)

write.table(df.GO, file='181006_subClust_main3_GO_list', quote=FALSE, row.names=FALSE)



# genes of beliveau interest
# Rif1, Rif2, Rap1

data[data$strain.common %in% c('RIF1', 'RIF2', 'RAP1', 'SIR2', 'SIR3', 'SIR4', 'STN1', 'CDC13', 'TEN1', 'YKU70', 'YKU80',
                               'EST1','EST2','EST3', 'TLC1'),]

data[data$Cluster == '1-15',]


write.table(data[data$Cluster == '3-2',]$strain, quote=FALSE, row.names=FALSE)
write.table(data[data$Cluster == '3-9',]$strain, quote=FALSE, row.names=FALSE)
write.table(data[data$Cluster == '3-20',]$strain, quote=FALSE, row.names=FALSE)
