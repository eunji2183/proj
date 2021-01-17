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
