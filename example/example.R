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
anno = read.delim("flatfile_human_gene2HDO.csv",skip=1,sep="\t",header=F)
anno = as.data.frame(anno)

## load clustering and annotation data
rEnrich::load(x=membership,anno=anno)

## run enrichment analysis on loaded data
rEnrich::run()

## get enrichment values 
res = rEnrich::getResults()

