##reset
rm(list=ls())

library(rEnrich)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## build cluster membership data.frame
ids  = igraph::get.vertex.attribute(gg,"name",V(gg))
coms = igraph::get.vertex.attribute(gg,"louvain",V(gg))

membership = as.data.frame(cbind(ids,coms))

## load annotation flat file
anno1 = read.delim("flatfile.go.BP.csv",skip=1,sep="\t",header=F)
anno1 = as.data.frame(anno1)

anno2 = read.delim("flatfile_human_gene2HDO.csv",skip=1,sep="\t",header=F)
anno2 = as.data.frame(anno2)



## load clustering and annotation data
rEnrich::load(x=membership,
              anno1=anno1,
              anno2=anno2,
              anno3=NULL)

## run enrichment analysis on loaded data
rEnrich::run()

## get enrichment values 
res = rEnrich::getResults()

