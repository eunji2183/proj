rm(list = ls())
setwd("/home/eunji/proj/200612_microdust/")
path <- "./data/featurecounts/con_vs_upm"

source("./code/DEG_functions.R")
library(pheatmap)
library(magrittr) #%>%
library(dplyr)
library(ggplot2) #ggplot
library(ggrepel) #geom_label_repel
library(TCC)
library(DESeq2)
library(limma)
options(scipen = 999) #e remove


gtf <- rtracklayer::import("/home/eunji/proj/0_sh/ref/rna/Homo_sapiens.GRCh38.96.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()
save(ID,file="./data/ID.RData")

group <- rep(c("con","upm"),c(2,2))
count <- data(path = path)

#DESeq2
resdata <- DESeq2.1(count = count)
write.csv(resdata,file="./result/res.csv")

ggplot(data = resdata, aes(x=log2FoldChange,y=-log10(padj),col = DEG)) +
  geom_point(alpha=0.5)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept=2 ,linetype=4) +
  geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "blue", "black"), limits = c("Up", "Down", "not")) +
  geom_label_repel(aes(label=sign), fontface="bold", color="grey50", box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), segment.colour = "grey50")

ggplot(resdata, aes(-log10(padj),log2FoldChange, col = DEG)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(name = "", values = c("red", "blue", "black"), limits = c("Up", "Down", "not")) +
  theme_bw() +
  ggtitle("DEG(DESeq2) plot")


# m.value=log2FC  q.value=FDR(padj)
tccRES <- DESeq2.2(count=count)
write.csv(tccRES, "./result/DESeq2_table.csv")

ggplot(tccRES, aes(a.value, m.value, col = estimatedDEG)) +
  geom_point(alpha = 0.5) +
  scale_color_gradientn(colours = c("gray", 'red')) +
  theme_bw() +
  ggtitle("DEG(DESeq2) plot")

gene <- c("EIF2AK3","HSPA5","ATF6","ERN1",
          "BANF1","EXOSC10","EXOSC4","EIF4A1","NOLC1","NHP2","TTC37","CXXC1","ATF3","DNAJC3",
          "PDIA6","DDIT4","HYOU1","FUS","ALDN18A1","IMP3","GEMIN4",
          "MEF2C","NFATC2","SIK1","CAMK2A","IL12RB1","CYP1A1","NQO1",
          "TNF","F2RL3","CXCR2","CCL3","XBP1","GRIA4","FPR1",
          "GCLM","HSP90B1","TXNRD1","CALR","CRELD2")
heatmap(df=resdata)


geneset <- GSA.read.gmt("./data/GSEA/geneset.gmt")
UPR <- geneset$genesets 
UPR <- UPR[[2]]
heat3 <- heat %>%
  dplyr::filter(rownames(heat) %in% UPR)

BiocManager::install("GSEABase")
library(GSEABase)
library(GSA)
gmt <- getGmt("./data/GSEA/msigdb.v7.1.symbols.gmt")
S <- GSA.read.gmt("./data/GSEA/geneset.gmt")
len_vec=c()
len_vec[1]=3
for(i in 1:length(S$genesets)){len_vec[i] <- c(length(S$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(S$geneset.names,len_vec),use.names = F)
S7_1 <- as.data.frame(cbind(geneset=pathway_vec,symbol=unlist(S$genesets,use.names = F)))

geneset <- S7_1[grep("OXIDATIVE",S7_1$geneset),]
geneset <- geneset[grep("STRESS",geneset$geneset),]

gsea <- resdata %>%
  dplyr::select(gene_name,log2FoldChange)
write.table(gsea,file = "./result/gsearank.rnk",sep = "\t",row.names = F,col.names = F,quote = F)



gsea <- gsea[order(gsea$log2FoldChange,decreasing = T),]
gsea <- gsea %>%
  dplyr::mutate(rank = rank(log2FoldChange, ties.method = "random"))

gsea <- gsea[order(gsea$rank),]
gene_list <- gsea$log2FoldChange
names(gene_list) <- as.character(gsea$gene_name)
library(clusterProfiler)
GSEA <- GSEA(gene_list,exponent = 1,TERM2GENE = S7_1,TERM2NAME = NA,nPerm = 1000,minGSSize = 15,maxGSSize = 500,pvalueCutoff = 0.25,pAdjustMethod = "BH",verbose = T,seed = F)
