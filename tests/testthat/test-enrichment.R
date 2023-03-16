library(igraph)
file <- system.file("extdata", "PPI_Presynaptic.gml", package = "rEnrich")
gg = read.graph(file,format="gml")
## build cluster membership data.frame
ids  = get.vertex.attribute(gg,"name",V(gg))
coms = get.vertex.attribute(gg,"louvain",V(gg))
membership = as.data.frame(cbind(ids,coms))
## load annotation flat file
file <- system.file("extdata", "flatfile_human_gene2HDO.csv", package = "rEnrich")
anno = read.delim(file,skip=1,sep="\t",header=F)
anno = as.data.frame(anno)
wm<-membership
wa<-anno
names(membership)<-c('names','membership')
names(anno)<-c('termID','termName','names')
test_that("Proper column names are submitted",{
    expect_error(rEnrich::run_enrichment(wm, wa, "name"),'.*Membership.*')
    expect_error(rEnrich::run_enrichment(membership, wa, "name"),'.*Annotation*')
})

test_that("Shape is correct",{
    r<-rEnrich::run_enrichment(membership, anno, "name")
    expect_equal(c(180,17),dim(r))
    r<-rEnrich::run_enrichment(membership, anno, "name",usePrintAlt = FALSE)
    expect_equal(c(180,15),dim(r))
})

