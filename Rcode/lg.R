load("/HDD8T/eunji/proj/multi/tgdata.RData")

lgdata <- list()
for ( i in c(names(lgdata))){
  lgdata[[i]]$tpmcount<- tgdata[[i]]$tpmcount 
  lgdata[[i]]$coldata <- tgdata[[i]]$coldata}
i="ACC"
for(i in c(names(lgdata))){
  n <- lgdata[[i]]$coldata %>% dplyr::filter(group=="T")
  lgdata[[i]]$coldata <- n
  lgdata[[i]]$tpmcount <- lgdata[[i]]$tpmcount[,c(rownames(lgdata[[i]]$coldata))]}
save(lgdata,file = "/HDD8T/eunji/proj/lg/lgdata.RData")

h="ACC"
library(Hmisc)
load("/HDD8T/eunji/proj/2D3D/ID.RData")
pc <- ID %>% dplyr::filter(gene_biotype == "protein_coding")
pc <- pc[!duplicated(pc$gene_name),]
pc <- data.frame(row.names = pc$gene_name,gene_id=pc$gene_id)
for(h in c(names(lgdata))){
  a <- lgdata[[h]]$tpmcount
  a <- merge(pc,a,by="row.names")
  a <- data.frame(row.names = a$Row.names,a[,-c(1:2)])
  b <- as.data.frame(t(a))
  res <- rcorr(as.matrix(b))
  flattenCorrMatrix <- function(cormat, pmat){
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  res2 <- flattenCorrMatrix(res$r, res$P)
  lgdata[[h]][['pccor']] <- res2}
save(lgdata,file = "/HDD8T/eunji/proj/lg/lgdata.RData")
