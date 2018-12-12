setwd("~/Dropbox/R")
source("ignition.txt")
options(stringsAsFactors = FALSE)
library(DESeq)
library(monocle)
library(RColorBrewer)
library(tidyverse)


setwd('~/Dropbox/q_lab/projects/35_yeastTranscriptome_pathways/180913_monocle_lauren/')
data.main = read.csv('deltxn_UMAP_coords_50clust.csv')
data.sub = read.table('../181001_subClustering/181003_50clust_1-10subclust.txt', header=TRUE)
data.sub$sub_clust = paste(data.sub$main_cluster, data.sub$Cluster, sep='-')

setwd('~/Dropbox/q_lab/projects/35_yeastTranscriptome_pathways/yeast_annotations/')
pcs = read.delim('CYC2008_complex.tab.txt')
#pcs = pcs[pcs$Method == 'Affinity Capture-MS',]

setwd('~/Dropbox/q_lab/projects/35_yeastTranscriptome_pathways/')
#gene.names = read.table('181002_allStrains_reNamed.txt')
gene.names = read.csv('181003_allStrains_ordered_wSystematic.csv', header=FALSE)
gene.names = gene.names[,-3]
colnames(gene.names) = c('common', 'systematic')


# can rename colnames of cds now
data.main$strain = gene.names$common
data.main$strain.sys = gene.names$systematic


data = merge(data.main, pcs[,c('ORF','Complex')], by.x='strain.sys', by.y='ORF', all.x=TRUE)


# first pass of plotting complexes on umap
data.comps = as.data.frame(table(data$Complex))
head(data.comps[order(-data.comps$Freq),],20)

complex_of_interest = 'HOPS complex'

write.table(na.omit(data[data$Complex == complex_of_interest,]), quote=FALSE, row.names=FALSE, sep='\t')

#data = data.main
ggplot(data,aes(x=umap_1,y=umap_2,fill=factor(Cluster))) + geom_point(size=3, shape=21, colour='#888888') +
  scale_fill_manual(values=rep('gray',length(data$strain))) +
  geom_point(data=data[data$Complex == complex_of_interest,], aes(x=umap_1,y=umap_2, color=factor(Cluster)), fill='red',
             size=3, shape=21, colour='#888888') + 
  xlim(13.75,14.8) +
  ylim(13.75,14.95)


write.table(na.omit(data[data$Complex == complex_of_interest,]), quote=FALSE, row.names=FALSE, sep='\t')


cluster_of_interest = 11
data[data$Cluster == cluster_of_interest,]

ggplot(data,aes(x=umap_1,y=umap_2,fill=factor(Cluster))) + geom_point(size=3, shape=21, colour='#888888') +
  scale_fill_manual(values=rep('gray',length(data$strain))) +
  geom_point(data=data[data$Cluster == cluster_of_interest,], aes(x=umap_1,y=umap_2, color=factor(Cluster)), fill='red',
             size=3, shape=21, colour='#888888') + 
  xlim(13.75,14.8) +
  ylim(13.75,14.95)





dude = merge(gene.names, pcs, by.x='systematic', by.y='ORF')




dude[order(dude$Complex),]
