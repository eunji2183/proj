#PAAD DEG - 2D 3D 
# drug - gemcitabine 
setwd("/home/eunji/proj/PAAD/")

install.packages("data.table") 
install.packages("gplots")
install.packages("corrplot")
BiocManager::install("TCC")
BiocManager::install("apeglm")
install.packages('readxl')

library(readxl)
library(DESeq2)
library(magrittr)
library(dplyr)
library(survminer)
library(survival)
library(RTCGAToolbox) 
library(ggplot2)
library(readr)
library(ffbase)
library(UCSCXenaTools)
library(asaur)
library(reshape)
library(ggpubr)
library(stringr)
library(survminer)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(GenomicDataCommons)
library(CePa)
library(readr)
library(tidyr)
library(tibble)
library(pheatmap)
library(VennDiagram)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)
library(numDeriv)
library(reshape2)
library(UCSCXenaTools)
library(rtracklayer) #gtf file 
library(utils)
library(data.table) #fread
library(ggrepel) #geom_label_repel
library(TCC)
library(limma)
library(apeglm)
library(gplots)
library(corrplot)


options(scipen = 999) #e remove
options(scipen = -100)
STAR_count <- as.data.frame(fread("./data/RNA/count/STAR/count.txt",fill = T,header = T))
featurecount <- as.data.frame(fread("./data/RNA/count/featurecount/count.txt",fill = T,header = T))
hisat <- as.data.frame(fread("./data/RNA/count/HISAT2_F/count.txt",fill = T,header = T))
sample <- read.table("./data/RNA/config2")
tsample <- as.data.frame(t(sample))
write.table(tsample,file = "./data/RNA/count/featurecount/sample.txt",sep = "\t",row.names = F,col.names = F,quote = F)

#GTF
gtf <- rtracklayer::import("./data/RNA/ref/Homo_sapiens.GRCh38.102.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()
save(ID,file="./data/ID.RData")
group <- as.factor(rep(c("NR","R"),c(3,3)))


names(featurecount)[1] <- 'gene_id'
All <- merge(ID,featurecount,by="gene_id")
TCGA <- data.frame(gene_name=rownames(merge3),merge3)
rownames(TCGA) <- NULL
save(TCGA,file = "./TCGA.RData")

All <- merge(All,TCGA,by="gene_name")
All <- All[,-2]
All <- All[!duplicated(All$gene_name),]
save(All,file = "./All.RData")

meta <-as.data.frame(read_excel("./data/paper/mmc2.xlsx",
                                sheet = "B6",
                                col_names = T,na="NA"))
meta$Gene.Symbol <- toupper(meta$Gene.Symbol)
All$SYMBOL <- toupper(All$SYMBOL)
names(meta)[1] <- 'SYMBOL'



rownames(All) <- All[,1]
All <- All[,-1]





gene <- read.table("./geneset/Gempathway.txt")

annotation_row = data.frame(gene_name = gene$V1,
                            GeneGroup=factor(rep(c("KEGG_HEDGEHOG_SIGNALING_PATHWAY","KEGG_WNT_SIGNALING_PATHWAY",
                                                   "KEGG_NOTCH_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY",
                                                   "BIOCARTA_AKT_PATHWAY","BIOCARTA_NFKB_PATHWAY",
                                                   "PID_HIF1A_PATHWAY","KEGG_APOPTOSIS",
                                                   "KEGG_ABC_TRANSPORTERS","BIOCARTA_EGF_PATHWAY",
                                                   "BIOCARTA_TGFB_PATHWAY","KEGG_VEGF_SIGNALING_PATHWAY",
                                                   "REACTOME_RAF_ACTIVATION","BIOCARTA_ERK_PATHWAY",
                                                   "REACTOME_PI3K_CASCADE","KEGG_MTOR_SIGNALING_PATHWAY",
                                                   "BIOCARTA_MTOR_PATHWAY","BIOCARTA_TEL_PATHWAY",
                                                   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_KRAS_SIGNALING_UP",
                                                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
                                                 c(56,151,47,267,22,21,19,87,44,27,19,76,34,27,68,52,22,17,200,200,200))))
annotation_row <- annotation_row[!duplicated(annotation_row$gene_name),]

gene <- data.frame(gene_name=annotation_row$gene_name)



# 2D vs 3D R DEG 
R23 <- All[,c(1:3,7:9)]
R23col <- PCA_col[c(1:3,7:9),]
R23I <- as.matrix(sapply(R23,as.numeric))
R23I[is.na(R23I)] <- 0
row.names(R23I) <- rownames(R23)
dds <- DESeqDataSetFromMatrix(countData = R23I,colData = R23col,design = ~SAMPLE)
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="SAMPLE_3D_R_vs_2D_R", type="apeglm")
res <- as.data.frame(resLFC)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),
                 by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 
                                & abs(resdata$log2FoldChange) > 1,
                                ifelse(resdata$log2FoldChange > 1 , 
                                       'UP','DOWN'),'NOT'))
names(resdata)[1] <- 'gene_name'
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1.5,
                       resdata$gene_name,NA)
R23DEG <- resdata
save(R23DEG,file = "./R23DEG.RData")
ggplot(data = R23DEG, aes(x=log2FoldChange,y=-log10(padj),col = DEG)) +
  geom_point(alpha=0.5)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept=2 ,linetype=4) +
  geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "blue", "black"), 
                     limits = c("UP", "DOWN", "NOT")) +
  geom_label_repel(aes(label=sign), fontface="bold", color="grey50", 
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"),
                   segment.colour = "grey50")

R23DEGs <- R23DEG %>%
  dplyr::filter(DEG == "UP" | DEG == "DOWN")
gene <- R23DEGs$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")



# 2D vs 3D NR 

NR23 <- All[,c(4:6,10:12)]
NR23col <- PCA_col[c(4:6,10:12),]
NR23I <- as.matrix(sapply(NR23,as.numeric))
NR23I[is.na(NR23I)] <- 0
row.names(NR23I) <- rownames(NR23)
dds <- DESeqDataSetFromMatrix(countData = NR23I,colData = NR23col,design = ~SAMPLE)
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="SAMPLE_3D_NR_vs_2D_NR", type="apeglm")
res <- as.data.frame(resLFC)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),
                 by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1,
                                ifelse(resdata$log2FoldChange > 1 , 
                                       'UP','DOWN'),'NOT'))
names(resdata)[1] <- 'gene_name'
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1.5,
                       resdata$gene_name,NA)
NR23DEG <- resdata
save(NR23DEG,file = "./res2/NR23DEG.RData")

ggplot(data = NR23DEG, aes(x=log2FoldChange,y=-log10(padj),col = DEG)) +
  geom_point(alpha=0.5)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept=2 ,linetype=4) +
  geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "blue", "black"), 
                     limits = c("UP", "DOWN", "NOT")) +
  geom_label_repel(aes(label=sign), fontface="bold", color="grey50", 
                   box.padding=unit(0.35, "lines"), 
                   point.padding=unit(0.5, "lines"), segment.colour = "grey50")

NR23DEGs <- NR23DEG %>%
  dplyr::filter(DEG == "UP" | DEG == "DOWN")
gene <- NR23DEGs$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")

overlapgene <- rbind(R23DEGs[,1:6],NR23DEGs[,1:6])
overlapgene <- overlapgene[duplicated(overlapgene$gene_name),]
overlapgene <- data.frame(gene_name=overlapgene$gene_name)
gene <- overlapgene$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
overlapres <- as.data.frame(go@result)

names(All)[1] <- 'SYMBOL'

meta <-as.data.frame(read_excel("./data/paper/mmc2.xlsx",
                                sheet = "B6",
                                col_names = T,na="NA"))
meta$Gene.Symbol <- toupper(meta$Gene.Symbol)
All$SYMBOL <- toupper(All$SYMBOL)
names(meta)[1] <- 'SYMBOL'

selected_gene <- data.frame(SYMBOL=
                              c("SLC1A5","SLC38A1","GLS","GOT1",
                                "SLC38A5","SLC7A11","GPT2","PPAT","GLS2",
                                "HIF1A","VEGFA","MYC","CCND1","EPAS1","PDK1",
                        "SNAI2","VIM","ZEB1","CDH2","CD44","TWIST1","ALDH1A3"),
                        GeneGroup=factor(rep(c("Glutaminolysis","HYPOXIA","EMT"),
                                             c(9,6,7)))
                            )
filterAll <- filter(All,(SYMBOL %in% c(selected_gene$SYMBOL)) == TRUE)

rownames(filterAll) <- filterAll[,1]
filterAll <- filterAll[,-1]
colnames(filterAll) <- gsub('.','-',colnames(filterAll),fixed = T)



Allcol <- data.frame(row.names = colnames(filterAll),
                      RESPONSE=factor(rep(c("NR","R","NR","R","R","NR"),
                                          c(3,3,3,3,17,24),levels=c("R","NR"))),
                      SAMPLE=factor(rep(c("2D_NR","2D_R","3D_NR","3D_R","TCGA_R","TCGA_NR"),
                                        c(3,3,3,3,17,24))))

#PCA
AllI <- as.matrix(sapply(filterAll,as.numeric))
AllI[is.na(AllI)] <- 0
row.names(AllI) <- rownames(filterAll)
dds <- DESeqDataSetFromMatrix(countData = AllI,colData = Allcol,design = ~RESPONSE)
dds <- DESeq(dds)
vsd <- vst(dds, blind=TRUE,nsub=nrow(dds))
vsd <- vst(dds)
pcaData <- plotPCA(vsd,intgroup="RESPONSE",returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))
ggplot(pcaData,aes(PC1,PC2,color=RESPONSE,shape=RESPONSE)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  geom_label_repel(aes(label=name), 
                   fontface="bold", color="grey50", box.padding=unit(0.35, "lines"), 
                   point.padding=unit(0.2, "lines"), segment.colour = "grey50",size=02)

norm <- as.data.frame(counts(dds,normalized=TRUE))


# 2D R vs NR 

RNR2D <- filterAll[,1:6]
RNR2Dcol <- Allcol[1:6,]
RNR2DI <- as.matrix(sapply(RNR2D,as.numeric))
RNR2DI[is.na(RNR2DI)] <- 0
row.names(RNR2DI) <- rownames(RNR2D)
dds <- DESeqDataSetFromMatrix(countData = RNR2DI,colData = RNR2Dcol,design = ~RESPONSE)
dds <- DESeq(dds)
RNR2Dnorm <- as.data.frame(counts(dds,normalized=TRUE))
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="RESPONSE_R_vs_NR", type="apeglm")
res <- as.data.frame(resLFC)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),
                 by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1,
                                ifelse(resdata$log2FoldChange > 1 , 
                                       'UP','DOWN'),'NOT'))
names(resdata)[1] <- 'gene_name'
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1.5,
                       resdata$gene_name,NA)
RNR2DEG <- resdata
save(RNR2DEG,file = "./res4/RNR2DEG.RData")

RNR2DEGs <- RNR2DEG %>%
  dplyr::filter(padj < 0.05)
gene <- RNR2DEGs$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")

anno <- read.table(file = "./data/mart_export.txt",sep = "\t",header = T,fill = T,quote = "")
anno <- anno %>%
  dplyr::select(Gene.name,Gene.description,Gene.type) %>% distinct()
names(anno)[1] <- 'gene_name'
save(anno,file = "./res3/anno.RData")

RNR2DEGs <- merge(RNR2DEGs,anno,by="gene_name")
RNR2DEGs <- RNR2DEGs[order(RNR2DEGs$pvalue),]
write.table(RNR2DEGs,file = "./res4/RNR2DEGs.csv",sep = ",",quote = F,row.names = F)


# 3D R vs NR 
RNR3D <- filterAll[,7:12]
RNR3Dcol <- Allcol[7:12,]
RNR3DI <- as.matrix(sapply(RNR3D,as.numeric))
RNR3DI[is.na(RNR3DI)] <- 0
row.names(RNR3DI) <- rownames(RNR3D)
dds <- DESeqDataSetFromMatrix(countData = RNR3DI,colData = RNR3Dcol,design = ~RESPONSE)
dds <- DESeq(dds)
RNR3Dnorm <- as.data.frame(counts(dds,normalized=TRUE))
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="RESPONSE_R_vs_NR", type="apeglm")
res <- as.data.frame(resLFC)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),
                 by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1,
                                ifelse(resdata$log2FoldChange > 1 , 
                                       'UP','DOWN'),'NOT'))
names(resdata)[1] <- 'gene_name'
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1.5,
                       resdata$gene_name,NA)
RNR3DEG <- resdata
save(RNR3DEG,file = "./res4/RNR3DEG.RData")

RNR3DEGs <- RNR3DEG %>%
  dplyr::filter(padj < 0.05)
gene <- RNR3DEGs$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
RNR3DEGs <- merge(RNR3DEGs,anno,by="gene_name")
RNR3DEGs <- RNR3DEGs[order(RNR3DEGs$pvalue),]
write.table(RNR3DEGs,file = "./res4/RNR3DEGs.csv",sep = ",",quote = F,row.names = F)

#TCGA R vs NR 
RNRTCGA <- filterAll[,13:53]
RNRTCGAcol <- Allcol[13:53,]
RNRTCGAI <- as.matrix(sapply(RNRTCGA,as.numeric))
RNRTCGAI[is.na(RNRTCGAI)] <- 0
row.names(RNRTCGAI) <- rownames(RNRTCGA)
dds <- DESeqDataSetFromMatrix(countData = RNRTCGAI,colData = RNRTCGAcol,design = ~RESPONSE)
dds <- DESeq(dds)
RNRTCGAnorm <- as.data.frame(counts(dds,normalized=TRUE))
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="RESPONSE_R_vs_NR", type="apeglm")
res <- as.data.frame(resLFC)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),
                 by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1,
                                ifelse(resdata$log2FoldChange > 1 , 
                                       'UP','DOWN'),'NOT'))
names(resdata)[1] <- 'gene_name'
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)> 1.5,
                       resdata$gene_name,NA)
RNRTCGADEG <- resdata
save(RNRTCGADEG,file = "./res4/RNRTCGADEG.RData")

RNRTCGADEGs <- RNRTCGADEG %>%
  dplyr::filter(padj < 0.05)
gene <- RNRTCGADEGs$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
RNRTCGADEGs <- merge(RNRTCGADEGs,anno,by="gene_name")
RNRTCGADEGs <- RNRTCGADEGs[order(RNRTCGADEGs$pvalue),]
write.table(RNRTCGADEGs,file = "./res4/RNRTCGADEGs.csv",sep = ",",quote = F,row.names = F)


gene <- data.frame(gene_name=rownames(Allnorm))
annorow <- merge(gene,annotation_row,by="gene_name")
annorow <- data.frame(row.names = selected_gene$SYMBOL,GeneGroup=selected_gene$GeneGroup)

#heatmap 
Allnorm <- merge(RNR2Dnorm,RNR3Dnorm,by="row.names")
rownames(Allnorm) <- Allnorm[,1]
Allnorm <- Allnorm[,-1]
Allnorm <- merge(Allnorm,RNRTCGAnorm,by="row.names")
rownames(Allnorm) <- Allnorm[,1]
Allnorm <- Allnorm[,-1]
Allnorm <- merge(annorow,Allnorm,by="row.names")
Allnorm <- Allnorm[order(Allnorm$GeneGroup),]
Allnorm <- Allnorm[,-2]
rownames(Allnorm) <- Allnorm[,1]
Allnorm <- Allnorm[,-1]
pheatmap(Allnorm,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize_col = 10, fontsize_row=10,
         show_rownames = T,cluster_rows = F,annotation_col = Allcol,
         annotation_row = annorow,legend = T)

===============================================================================================================================



#TPM normalization - FC & log2FC - p.value(2D 3D) two side t-test
# FC 가 2D <1 , 3D & T > 1 인 gene 들이 어떤 pathway 에서 enrichment 되였는지 (p <0.1)
# 

BiocManager::install("pheatmap")
BiocManager::install("sva")
BiocManager::install("genefilter")
devtools::install_github('dviraran/xCell')
install.packages("psych")
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xCell)
library(psych)
library(sva)
library(genefilter)



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

merge$FC2DNR <-  apply(ss[,1:6],1,function(x){(x-mean(x))/sd(x)})

#FC for 2D 3D patient (R / NR) 
for (i in 1:nrow(merge)) {
  a <- merge[i,8:10]
  b <- merge[i,1:3]
  c <- mean(as.matrix(a))
  d <- mean(as.matrix(b))
  merge$FC2D[i] <- c/d 
}

for (i in 1:nrow(merge)) {
  a <- merge[i,11:13]
  b <- merge[i,4:6]
  c <- mean(as.matrix(a))
  d <- mean(as.matrix(b))
  merge$FC3D[i] <- c/d 
}

for (i in 1:nrow(merge)) {
  a <- merge[i,14]
  b <- merge[i,7]
  c <- mean(as.matrix(a))
  d <- mean(as.matrix(b))
  merge$FCT[i] <- c/d 
}

merge$log2FC2D = apply( merge, 1, function(t) {log2( mean(t[8:10])/ mean(t[1:3]))})
merge$log2FC3D = apply( merge, 1, function(t) {log2( mean(t[11:13])/ mean(t[4:6]))})
merge$log2FCT = apply( merge, 1, function(t) {log2( mean(t[14])/ mean(t[7]))})


#p value

tpmdf <- merge[,c(1:14)]
log2tpmdf <- log2(tpmdf)
log2tpmdf <- 
  
  pvalue = padj = log2FoldChange = matrix(0, nrow(tpmdf), 1)
for(i in 1:nrow(tpmdf)){
  x=as.numeric(tpmdf[i, 8:10])
  y=as.numeric(tpmdf[i, 1:3])
  tpmdf$p.2D[i]=t.test(x, y,alternative = "two.sided")$p.value}

for(i in 1:nrow(tpmdf)){
  x=as.numeric(tpmdf[i, 11:13])
  y=as.numeric(tpmdf[i, 4:6])
  tpmdf$p.3D[i]=t.test(x, y,alternative = "two.sided")$p.value}


#protein_coding gene 

coding <- merge(coding_ID,tpms,by="gene_id")
coding <- coding[,-1]
coding <- coding[!duplicated(coding$gene_name),]
rownames(coding) <- coding[,1]
coding <- coding[,-1]
coding <- coding[,c(1:3,7:9,14,4:6,10:12,13)]
names(coding) <- rownames(Allcol)
merge <- coding






metabolic_gene <- read_excel("./2D3D/geneset/1-s2.0-S1550413120305453-mmc2.xlsx",sheet = 1)

names(metabolic_gene) <- c('gene_name','GENESET')

save(metabolic_gene,file = "./metabolic_gene.RData")

save(merge,file = "./merge.RData")

#2D 

two <- merge %>%
  dplyr::select(`2D_NR1`,`2D_NR2`,`2D_NR3`,
                `2D_R1`,`2D_R2`,`2D_R3`,FC2D,log2FC2D)

two <- data.frame(gene_name=rownames(two),two)
rownames(two) <- NULL

two_coding_GSEA <- two %>%
  dplyr::filter(FC2D != Inf) %>%
  dplyr::filter(!is.nan(FC2D)) %>%
  dplyr::select(gene_name,FC2D) 

two <- merge(metabolic_gene,two,by="gene_name")

two_up <- two %>%
  dplyr::filter(log2FC2D > 0 & log2FC2D != Inf) %>%
  dplyr::select(gene_name)

two_down <- two %>%
  dplyr::filter(log2FC2D < 0 & log2FC2D != -Inf )%>%
  dplyr::select(gene_name)


#3D 
three <- merge %>%
  dplyr::select(`3D_NR1`,`3D_NR2`,`3D_NR3`,
                `3D_R1`,`3D_R2`,`3D_R3`,FC3D,log2FC3D)
three <- data.frame(gene_name=rownames(three),three)
rownames(three) <- NULL 

three <- merge(metabolic_gene,three,by="gene_name")

three_up <- three %>%
  dplyr::filter(log2FC3D > 0 & log2FC3D != Inf) %>%
  dplyr::select(gene_name)

three_down <- three %>%
  dplyr::filter(log2FC3D < 0 & log2FC3D != -Inf ) %>%
  dplyr::select(gene_name)

#patient 

p <- merge %>%
  dplyr::select(Y5_NR,Y26_R,FCT,log2FCT)
p<- data.frame(gene_name=rownames(p),p)
rownames(p) <- NULL 

p <- merge(metabolic_gene,p,by="gene_name")

p_up <- p %>%
  dplyr::filter(log2FCT > 0 & log2FCT != Inf) %>%
  dplyr::select(gene_name)

p_down <- p %>%
  dplyr::filter(log2FCT < 0 & log2FCT != -Inf ) %>%
  dplyr::select(gene_name)


#venn diagram 


two_up <- two_up$gene_name 
three_up <- three_up$gene_name
p_up <- p_up$gene_name

two_down <- two_down$gene_name
three_down <- three_down$gene_name
p_down <- p_down$gene_name
x = list(two_up,three_up,p_up)

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

VennDiagram::draw.triple.venn(
  x,
  category.names = c("2D" , "3D " , "patient"),
  filename = 'all_up.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(list(twoD=two_up,threeD=three_up,patient=p_up), filename=NULL, fill=rainbow(3))  

#FC 



FC <- merge %>%
  dplyr::select(log2FC2D,log2FC3D,log2FCT,FC2D,FC3D,FCT)

FC_up <- FC %>%
  dplyr::filter(log2FC3D > 0 & log2FCT > 0) %>%
  dplyr::filter(log2FC2D < 0) %>%
  dplyr::filter(log2FCT != Inf & log2FC3D != Inf) %>%
  dplyr::filter(log2FC2D !=  -Inf)

FC_up <- data.frame(gene_name=rownames(FC_up),FC_up)
rownames(FC_up) <- NULL 

FC_up_met <- merge(metabolic_gene,FC_up,by="gene_name")

FC_up_met <- FC_up_met %>%
  dplyr::filter(log2FC3D > 1 & log2FCT > 1)
  
FC_down <- FC %>%
  dplyr::filter(log2FC3D < 0 & log2FCT < 0) %>%
  dplyr::filter(log2FC2D > 0) %>%
  dplyr::filter(log2FCT != -Inf & log2FC3D != -Inf) %>%
  dplyr::filter(log2FC2D !=  Inf)

FC_down <- data.frame(gene_name=rownames(FC_down),FC_down)
rownames(FC_down) <- NULL

FC_down_met <- merge(metabolic_gene,FC_down,by="gene_name")

FC_down_met <- FC_down_met %>%
  dplyr::filter(log2FC3D < -1 & log2FCT < -1)

FC_merge <- rbind(FC_up,FC_down)


save(FC_up,file = "./metabolic_FC_up.RData")
save(FC_down,file = "./metabolic_FC_down.RData")
save(FC_merge,file = "./metabolic_FC_merge.RData")



GOgenedown <- sig_down$gene_name
GOgeneup <- sig_up$gene_name

GOgenedown <- bitr(GOgenedown,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
GOgeneup <- bitr(GOgeneup,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

de1 <- GOgenedown$ENTREZID
go_down <- enrichGO(gene = de1, OrgDb = "org.Hs.eg.db",ont = "all")
de2 <- GOgeneup$ENTREZID
go_up <- enrichGO(gene = de2, OrgDb = "org.Hs.eg.db",ont = "all")

p1 <- barplot(go_down,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free") +ggtitle("GO_DOWN")

p2 <- barplot(go_up,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free") +ggtitle("GO_UP")
plot_grid(p1, p2, ncol=2)


#gene 골라서 heatmap 
GOgene <- c("AXL","CLK1","DYRK4","FES","JAK1")
GOheat <- filter(sig_down,(gene_name %in% GOgene == TRUE))
GOheat <- GOheat[,c(1,21,22,23)]
rownames(GOheat) <- GOheat[,1]
GOheat <- GOheat[,-1]
names(GOheat) <- c('2D','3D','Patient')

paletteLength <- 30
myColor <- colorRampPalette(c('navy','white','red'))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(GOheat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(GOheat)/paletteLength, max(GOheat), length.out=floor(paletteLength/2)))

pheatmap(GOheat, cluster_cols = F,cluster_rows = F,
         clustering_distance_row = "correlation",
         fontsize_col=12,fontsize_row=11, cellwidth = 30, cellheight = 15,
         color=myColor, breaks=myBreaks)

edo <- enrichDGN(de)
library(enrichplot)
barplot(edo, showCategory=10)

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

xCellgene <- data.frame(gene_name=xCell.data$genes)


#ssGSEA
count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
rownames(count) <- count[,1]
count <- count[,-1]
countData <- count[apply(count, 1, sum) > 1 , ]
names(countData)[13] <-'Y26_R' 
names(countData)[14] <- 'Y5_NR'
countData <- countData[,c(1:3,7:9,14,4:6,10:12,13)]
countData <- data.frame(gene_id=rownames(countData),countData)
rownames(countData) <- NULL
length <- merge(countData,geneLength,by="gene_id")
geneLength <- length[,c(1,16)]
rownames(countData) <- countData[,1]
countData <- countData[,-1]
names(countData) <- rownames(Allcol)
tpms <- apply(countData,2,function(x)tpm(x,geneLength$gene_length))
tpms <- as.data.frame(tpms)
tpms <- data.frame(gene_id=rownames(tpms),tpms)
rownames(tpms) <- NULL
merge <- merge(ID,tpms,by="gene_id")
merge <- merge[!duplicated(merge$gene_name),]
merge <- merge[,-1]
rownames(merge) <- merge[,1]
merge <- merge[,-1]
names(merge) <- rownames(Allcol)
boxplot(merge)
logTPM <- log2(merge+1)
boxplot(logTPM)
save(logTPM,file = "./logTPM.RData")

geneSet <- read.csv2("./2D3D/geneset/h.all.v7.2.symbols.gmt",header = F,sep = "\t")
geneSet <- geneSet[,-2]
geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
geneSet <- geneSet[,-c(4,5,8,11,12,15,19,20,22,30,34,43,51,55,58,61,65)]
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
save(l,file = "./gene_set.Rdata")

library(GSVA)
library(limma)
dat <- as.matrix(logTPM)

FCAll <- merge[,c(15:17)]
FCAll2 <- FCAll %>%
  dplyr::filter(FC2D != Inf) %>%
  dplyr::filter(!is.nan(FC2D)) %>%
  dplyr::filter(FC3D != Inf) %>%
  dplyr::filter(!is.nan(FC3D)) %>%
  dplyr::filter(FCT != Inf) %>%
  dplyr::filter(!is.nan(FCT))
dat <- as.matrix(FCAll2)


FCM <- rbind(FC_up,FC_down)
rownames(FCM) <- FCM[,1]
FCM <- FCM[,-1]
FCM <- FCM[,-c(1:3)]
dat <- as.matrix(FC_up)

sigM <- sig_up %>%
  dplyr::select(gene_name,FC2D,FC3D,FCT)
rownames(sigM) <- sigM[,1]
sigM <- sigM[,-1]
dat <- as.matrix(sigM)

sigM <- sig_down %>%
  dplyr::select(gene_name,FC2D,FC3D,FCT)
rownames(sigM) <- sigM[,1]
sigM <- sigM[,-1]
dat <- as.matrix(sigM)

FCm <- FC_merge
rownames(FCm) <- FCm[,1]
FCm <- FCm[,-1]
dat <- as.matrix(FCm)


ssgsea<- gsva(dat, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ss <- as.data.frame(ssgsea)
names(ss) <- c('2D','3D','Patient')

ss <- ss %>%
  rowwise() %>%
  mutate(mdNR = median(c(`2D_NR1`,`2D_NR2`,`2D_NR3`,`3D_NR1`,`3D_NR2`,`3D_NR3`))) %>%
  mutate(mNR = mean(c(`2D_NR1`,`2D_NR2`,`2D_NR3`,`3D_NR1`,`3D_NR2`,`3D_NR3`))) %>%
  mutate(mdR = median(c(`2D_R1`,`2D_R2`,`2D_R3`,`3D_R1`,`3D_R2`,`3D_R3`))) %>%
  mutate(mR = mean(c(`2D_R1`,`2D_R2`,`2D_R3`,`3D_R1`,`3D_R2`,`3D_R3`)))

ss <- as.data.frame(ss)
rownames(ss) <- rownames(ssgsea)  
  
ss1 <- apply(ss[,1:6],1,function(x){(x-mean(x))/sd(x)})  
ss2 <- apply(ss[,8:13],1,function(x){(x-mean(x))/sd(x)})
ss3 <- rbind(ss1,ss2)

ss4 <- apply(ss[,1:6],1,function(x){(x-median(x))/(sqrt(sum((x - median(x))^2) / 6))})  
ss5 <- apply(ss[,8:13],1,function(x){(x-median(x))/(sqrt(sum((x - median(x))^2) / 6))})
ss6 <- rbind(ss4,ss5)

ss$Y5_z <- (ss$Y5_NR-ss$mdNR)/ss$sdNR
ss$Y26_z <- (ss$Y26_R-ss$mdR)/(sqrt(sum((ss$Y26_R-ss$mdR)^2)/6))


for (i in rownames(ss)) {
  a <- ss[i,1:6]
  b <- ss[i,8:13]
  c <- apply(a[,1:6],function(x){sqrt(sum((x - median(x))^2) / 6)})
  d <- apply(b[,1:6],function(x){sqrt(sum((x - median(x))^2) / 6)})
  ss$sdNR <- c
  ss$sdR <- d
}

ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}

zmerge <- apply(ssgsea,1,function(x){(x-mean(x))/sd(x)})
zmerge <- t(ss6)
pheat <- data.frame(row.names = rownames(zmerge),zmerge[,1:6],Y5_NR=ss$Y5_z,zmerge[,7:12],Y26_R=ss$Y26_z)
names(pheat) <- rownames(Allcol)
library(pheatmap)

bk <- c(seq(-0.6,-0.001,by=0.02),seq(0,0.5,by=0.02))
bk2 <- c(seq(-0.6,-0.000001,by=0.04),seq(0,0.6,by=0.01))
pheatmap(pheat,
         annotation_col = Allcol,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         fontsize=10)
annocol <- data.frame(row.names=colnames(ss),SAMPLE=colnames(ss))
pheatmap(ss,
         
         show_colnames = T,
         cluster_rows = T,
         cluster_cols = F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk2)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk2)/2)),
         fontsize_col=10,fontsize_row=8, cellwidth = 20, cellheight = 9)
ssGSEA <- as.data.frame(ss)
write.csv(ssGSEA,file="./ssGSEA.csv")
heat <- read.csv(file = "./heat.csv",sep = ",",header = T,row.names = 1)
names(GOheat) <- c('2D','3D','Patient')


#p value

tpmdf <- merge[,c(1:14)]
log2tpmdf <- log2(tpmdf)
log2tpmdf <- 

pvalue = padj = log2FoldChange = matrix(0, nrow(tpmdf), 1)
for(i in 1:nrow(tpmdf)){
  x=as.numeric(tpmdf[i, 8:10])
  y=as.numeric(tpmdf[i, 1:3])
  tpmdf$pvalue[i]=t.test(x, y,alternative = "two.sided")$p.value}

for(i in 1:nrow(tpmdf)){
  x=as.numeric(tpmdf[i, 11:13])
  y=as.numeric(tpmdf[i, 4:6])
  tpmdf$p.3D[i]=t.test(x, y,alternative = "two.sided")$p.value}

siggene <- tpmdf %>%
  dplyr::filter(p.2D < 0.1 & p.3D < 0.1)
sig <- data.frame(gene_name=rownames(siggene),siggene)
rownames(sig) <- NULL
colnames(FC)[4:6]
name <- c('gene_name',rownames(Allcol),'p.2D','p.3D',colnames(FC)[4:6],colnames(FC)[1:3])
names(sig) <- name
study_gene <- data.frame(gene_name=
                              c("SLC1A5","SLC38A1","GLS","GOT1",
                                "SLC38A5","SLC7A11","GPT2","PPAT","GLS2",
                                "HIF1A","VEGFA","MYC","CCND1","EPAS1","PDK1",
                                "SNAI2","VIM","ZEB1","CDH2","CD44","TWIST1","ALDH1A3"),
                            GeneGroup=factor(rep(c("Glutaminolysis","HYPOXIA","EMT"),
                                                 c(9,6,7))))
study <- merge(study_gene,siggene,by="gene_name") 
rownames(study) <- study[,1]
study <- study[,-1]
study <- merge(study,merge,by="row.names")
study <- study[,-c(19:32)]

siggene <- merge(siggene,merge,by="row.names")
siggene <- siggene[,-c(18:31)]
names(siggene) <- name
siggene <-siggene %>%
  dplyr::filter(FC2D != Inf & log2FC2D != -Inf) %>%
  dplyr::filter(!is.nan(FC2D)) %>%
  dplyr::filter(FC3D != Inf & log2FC3D != -Inf) %>%
  dplyr::filter(!is.nan(FC3D)) %>%
  dplyr::filter(FCT != Inf & log2FCT != -Inf) %>%
  dplyr::filter(!is.nan(FCT))
sig_up <- siggene %>%
  dplyr::filter(log2FC3D > 0 & log2FCT > 0) %>%
  dplyr::filter(log2FC2D < 0)
sig_down <- siggene %>%
  dplyr::filter(log2FC3D < 0 & log2FCT < 0) %>%
  dplyr::filter(log2FC2D > 0)

sigmerge <- rbind(sig_up,sig_down)

GOgenedown <- sig_down$gene_name
GOgeneup <- sig_up$gene_name

GOgenedown <- bitr(GOgenedown,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
GOgeneup <- bitr(GOgeneup,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

de1 <- GOgenedown$ENTREZID
go_down <- enrichGO(gene = de1, OrgDb = "org.Hs.eg.db",ont = "all")
de2 <- GOgeneup$ENTREZID
go_up <- enrichGO(gene = de2, OrgDb = "org.Hs.eg.db",ont = "all")

p1 <- dotplot(go_down,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free") +ggtitle("GO_DOWN")

p2 <- barplot(go_up,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free") +ggtitle("GO_UP")
plot_grid(p1, p2, ncol=2)

GOmerge <- sigmerge$gene_name
GOmerge <- bitr(GOmerge,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de3 <- GOmerge$ENTREZID
go_merge <- enrichGO(gene = de3, OrgDb = "org.Hs.eg.db",ont = "all")
barplot(go_merge,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free") +ggtitle("GO")

#GSEA 

sigmerge2GSEA <- sigmerge %>%
  dplyr::select(gene_name,log2FC2D)
write.table(sigmerge2GSEA,file = "./sigmerge2GSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

sigmerge3GSEA <- sigmerge %>%
  dplyr::select(gene_name,log2FC3D)
write.table(sigmerge3GSEA,file = "./sigmerge3GSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

sigmergeTGSEA <- sigmerge %>%
  dplyr::select(gene_name,log2FCT)
write.table(sigmergeTGSEA,file = "./sigmergeTGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)




fcup2GSEA <- FC_up %>%
  dplyr::select(gene_name,log2FC2D)
write.table(fcup2GSEA,file = "./Fcup2GSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)
fcup3GSEA <- FC_up %>%
  dplyr::select(gene_name,log2FC3D)
write.table(fcup3GSEA,file = "./Fcup3GSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

fcupTGSEA <- FC_up %>%
  dplyr::select(gene_name,log2FCT)
write.table(fcupTGSEA,file = "./FcupTGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

fcdown2GSEA <- FC_down %>%
  dplyr::select(gene_name,log2FC2D)
write.table(fcdown2GSEA,file = "./Fcdown2GSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)
fcdown3GSEA <- FC_down %>%
  dplyr::select(gene_name,log2FC3D)
write.table(fcdown3GSEA,file = "./Fcdown3GSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

fcdownTGSEA <- FC_down %>%
  dplyr::select(gene_name,log2FCT)
write.table(fcdownTGSEA,file = "./FcdownTGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

#GSEA result pheatmap 
fcmerge <- read.table("/home/eunji/gsea_home/output/feb24/hallmark_fcdiff/fcmerge.tsv",
                      sep = "\t",header = T,row.names = 1)

names(fcmerge) <- c('2D','2D_p','3D','3D_p','Patient','Patient_p')
p <- fcmerge %>%
  dplyr::select(`2D`,`3D`,Patient)

p1 <- p %>%
  dplyr::filter(`3D` < 0 & Patient < 0) %>%
  dplyr::filter(`2D` > 0)
p2 <- p %>%
  dplyr::filter(`3D` > 0 & Patient > 0) %>%
  dplyr::filter(`2D` < 0)  

p <- rbind(p1,p2)

paletteLength <- 30
myColor <- colorRampPalette(c('navy','white','red'))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(p), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(p)/paletteLength, max(p), length.out=floor(paletteLength/2)))

pheatmap(p, cluster_cols = T,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=12,fontsize_row=11, cellwidth = 30, cellheight = 15,
         color=myColor, breaks=myBreaks)

fcmerge_kegg <- read.table("/home/eunji/gsea_home/output/feb24/kegg_fcdiff/merge.tsv",
                      sep = "\t",header = T,row.names = 1)

names(fcmerge_kegg) <- c('2D','2D_p','3D','3D_p','Patient','Patient_p')
k <- fcmerge_kegg %>%
  dplyr::select(`2D`,`3D`,Patient)

k1 <- k %>%
  dplyr::filter(`3D` < 0 & Patient < 0) %>%
  dplyr::filter(`2D` > 0)
k2 <- k %>%
  dplyr::filter(`3D` > 0 & Patient > 0) %>%
  dplyr::filter(`2D` < 0)  

k <- rbind(k1,k2)

paletteLength <- 30
myColor <- colorRampPalette(c('navy','white','red'))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(k)/paletteLength, max(k), length.out=floor(paletteLength/2)))

pheatmap(k, cluster_cols = T,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=12,fontsize_row=8, cellwidth = 20, cellheight = 8,
         color=myColor, breaks=myBreaks)

fcdown <- read.table("/home/eunji/gsea_home/output/feb24/hallmark_fcdiff/fcdown.tsv",
                      sep = "\t",header = T,row.names = 1)

names(fcdown) <- c('2D','2D_p','3D','3D_p','Patient','Patient_p')
d <- fcdown %>%
  dplyr::select(`2D`,`3D`,Patient)

d1 <- d %>%
  dplyr::filter(`3D` < 0 & Patient < 0) %>%
  dplyr::filter(`2D` > 0)
d2 <- d %>%
  dplyr::filter(`3D` > 0 & Patient > 0) %>%
  dplyr::filter(`2D` < 0)  

d <- rbind(d1,d2)

paletteLength <- 30
myColor <- colorRampPalette(c('navy','white','red'))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(d), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(d)/paletteLength, max(d), length.out=floor(paletteLength/2)))

pheatmap(d, cluster_cols = T,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=12,fontsize_row=11, cellwidth = 30, cellheight = 15,
         color=myColor, breaks=myBreaks)









twoGSEA <- sig_up %>%
  dplyr::select(gene_name,log2FC2D)
write.table(twoGSEA,file = "./twoGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

twoGSEA <- twoGSEA[order(twoGSEA$log2FC2D,decreasing = T),]
twoGSEA <- twoGSEA %>%
  dplyr::mutate(rank = rank(log2FC2D, ties.method = "random"))

twoGSEA <- twoGSEA[order(twoGSEA$rank,decreasing = T),]
gene_list <- twoGSEA$log2FC2D
names(gene_list) <- as.character(twoGSEA$gene_name)

threeGSEA <- sig_up %>%
  dplyr::select(gene_name,log2FC3D)
write.table(threeGSEA,file = "./threeGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

threeGSEA <- threeGSEA[order(threeGSEA$log2FC3D,decreasing = T),]
threeGSEA <-threeGSEA %>%
  dplyr::mutate(rank = rank(log2FC3D, ties.method = "random"))

threeGSEA <-threeGSEA [order(threeGSEA$rank,decreasing = T),]
gene_list <- threeGSEA$log2FC3D
names(gene_list) <- as.character(threeGSEA$gene_name)

pGSEA <- sig_up %>%
  dplyr::select(gene_name,log2FCT)
write.table(threeGSEA,file = "./threeGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

pGSEA <- pGSEA[order(pGSEA$log2FCT,decreasing = T),]
pGSEA <-pGSEA %>%
  dplyr::mutate(rank = rank(log2FCT, ties.method = "random"))

pGSEA <-pGSEA [order(pGSEA$rank,decreasing = T),]
gene_list <- pGSEA$log2FCT
names(gene_list) <- as.character(pGSEA$gene_name)

geneset <- GSA.read.gmt("./2D3D/geneset/msigdb.v7.2.symbols.gmt")
S <-geneset
len_vec=c()
len_vec[1]=3
for(i in 1:length(S$genesets)){len_vec[i] <- c(length(S$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(S$geneset.names,len_vec),use.names = F)
S7_1 <- as.data.frame(cbind(geneset=pathway_vec,symbol=unlist(S$genesets,use.names = F)))

GSEA <- GSEA(gene_list,exponent = 1,TERM2GENE = S7_1,TERM2NAME = NA,nPerm = 1000,minGSSize = 15,maxGSSize = 500,pvalueCutoff = 0.25,pAdjustMethod = "BH",verbose = T,seed = F)

sig2GSEAall <- GSEA@result
sig3GSEAall <- GSEA@result
sigpGSEAall 
UPR <- geneset$genesets 
UPR <- UPR[[2]]
install.packages("GSA")
library(GSA)





