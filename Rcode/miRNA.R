setwd("/home/eunji/proj/200709_micro/8.mirna/")
BiocManager::install("multiMiR")
BiocManager::install("miRNAtap")
BiocManager::install("miRNAtap.db")
library(miRNAtap.db)
library(miRNAtap)
library(multiMiR)
library(org.Hs.eg.db)
library(topGO)
library(DESeq2)
library(magrittr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(stringr)

mature <- read.table("./2.bam/all.counts.mature.txt",sep = "\t",header = T)
mature <- mature[,-c(2:6)]
mature <- mature[,-c(4:5)]
mature <- mature[,c(1,2,4,3,5)]
names(mature)[2] <- 'Control1'
names(mature)[3] <- 'Control2'
names(mature)[4] <- 'PM1'
names(mature)[5] <- 'PM2'
rownames(mature) <- mature[,1]
mature <- mature[,-1]
group <- rep(c("con","upm"),c(2,2))
coldata <- data.frame(row.names = colnames(mature),group)
dds <- DESeqDataSetFromMatrix(countData = mature,colData = coldata,design = ~group)
dds <- DESeq(dds)
res <- results(dds) 
resultsNames(dds) 
res <- as.data.frame(res)
res <- res %>% dplyr::filter(!is.na(log2FoldChange)) %>%
  dplyr::filter(padj < 0.25)
res <- res[order(res$pvalue,decreasing = F),]
write.csv(res,file = "./miRNA.csv")
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)

db.ver = multimir_dbInfoVersions()
db.ver[,1:3]
mirna <- rownames(res)
multimir_results <- get_multimir(org     = 'hsa',
                                 mirna   = mirna,
                                 table   = 'validated',
                                 summary = TRUE)
result <- multimir_results@data
result1 <- result[,3:5] %>% unique()
write.csv(result1,file = "./target.csv")
names(result1)[2] <- 'gene_name'

DEG <- read.csv("../2.DEG/featurecounts/count.csv",sep = ",",header = T,row.names = 1)
DEG <- DEG[,c(1,2,3,5,6)]
rownames(DEG) <- DEG[,1]
DEG <- DEG[,-1]
names(DEG)[1] <- 'Control1'
names(DEG)[2] <- 'Control2'
names(DEG)[3] <- 'PM1'
names(DEG)[4] <- 'PM2'
DEG_res <- DESeq2.1(count = DEG,group = group,ID=ID)
DEG_res2 <- DEG_res %>% 
  dplyr::filter(pvalue < 0.1)

GO <- merge(DEG_res,result1,by='gene_name') 
gene <- data.frame(gene_id=result1[,2]) %>% unique()
gtf <- rtracklayer::import("/home/eunji/ref//Homo_sapiens.GRCh38.96.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()
gene <- merge(ID,gene,by="gene_id")

GO_gene <- GO$gene_name
GO_gene <- bitr(GO_gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- GO_gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all",pvalueCutoff = 0.1)
EGG <- enrichKEGG(gene = de,organism = 'hsa',pvalueCutoff = 1,qvalueCutoff = 5)
KEGG <- EGG@result
dotplot(go,showCategory=10,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free") + ggtitle("dotplot")
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")
dotplot(EGG)
barplot(go)
