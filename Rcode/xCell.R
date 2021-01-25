# https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_19
BiocManager::install("pheatmap")
devtools::install_github('dviraran/xCell')
install.packages("psych")
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xCell)
library(psych)


setwd("/home/eunji/proj/2D3D/")
count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
geneLength <- as.data.frame(fread("./hg38_gene_length.txt",fill = T,header = F))
names(geneLength)[1] <- 'gene_id'
names(geneLength)[2] <- 'gene_length'

#TPM normalization 
rownames(count) <- count[,1]
count <- count[,-1]
geneLength <- geneLength[order(geneLength$gene_id),]
rownames(geneLength) <- NULL

tpm <- function(counts,lengths){
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk)/1e6
  rpk/coef
} 

tpms <- apply(count,2,function(x)tpm(x,geneLength$gene_length))
tpms <- as.data.frame(tpms)
tpms <- data.frame(gene_id=rownames(tpms),tpms)
rownames(tpms) <- NULL
merge <- merge(ID,tpms,by="gene_id")
merge <- merge[!duplicated(merge$gene_name),]
merge <- merge[,-1]
rownames(merge) <- merge[,1]
merge <- merge[,-1]
merge <- merge[,c(1:3,7:9,14,4:6,10:12,13)]
names(merge) <- rownames(Allcol)


#xCell 

xcell <- xCellAnalysis(merge,rnaseq = T)
fcs <- as.data.frame(xcell)
paad <- list(expr=merge,fcs=fcs,metadata=Allcol)
fcs = paad$fcs[rownames(xcell),colnames(xcell)]
res = corr.test(t(xcell),t(fcs),adjust='none')
qplot(x=rownames(res$r),y=diag(res$r),
      fill=diag(res$p)<0.1,geom='col',
      main='association with immunoprofiling',
      ylab='Pearson R', xlab=``) + labs(fill = "p-value<0.1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(pheatmap)
scaled.scores = t(scale(t(xcell)))
scaled.scores[scaled.scores > 3] = 2
scaled.scores[scaled.scores < -3] = -2
pheatmap(scaled.scores,show_colnames = F,
         annotation_col=paad$metadata[,
                                     c('RESPONSE','SAMPLE'),drop=F],
         clustering_method = 'ward.D2')
