##GEO&TCGA
rm(list = ls())
working_dir="/home/eunji/R/lab0610/"
setwd(working_dir)
options(stringsAsFactors = F)
save(list = ls(),file="lab0610.RData")
BiocManager::install("minfi")
BiocManager::install("AnnotationDbi")
BiocManager::install("phyloseq")
BiocManager::install("S4Vectors")
BiocManager::install("impute")
BiocManager::install("wateRmelon")
library(GEOquery) 
library(devtools)
library(UCSCXenaTools)
library(dplyr)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(data.table)
library(reshape)
library(purrr)
library(vsn)
library(hexbin)
library(apeglm)
library(grid)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(minfi)
library(stringr)
library(DT)
library(AnnotationDbi)
library(S4Vectors)
library(phyloseq)
library(impute)
library(wateRmelon)
##############################################################################
##GEO data 
gset3 <- getGEO('GSE78754',destdir = ".", AnnotGPL = F, getGPL = F) 
gset2 <- getGEO('GSE78751',destdir = ".", AnnotGPL = F, getGPL = F)
gset <- getGEO('GSE78758',destdir = ".", AnnotGPL = F, getGPL = F)
gset=gset[[1]]
gset2=gset2[[1]]
gset3=gset3[[1]]
pdata=pData(gset)
pdata2=pData(gset2)
pdata3=pData(gset3)
##############################################################################
#TCGA data : UCSC Xena, GDC TCGA BRCA - HTseq-count(log2(count+1)), phenotype
#phe<-read.csv('TCGA-BRCA.GDC_phenotype.tsv',sep = '\t')
expr<-read.csv('TCGA-BRCA.htseq_counts.tsv',sep = '\t')
colnames(expr)
rownames(expr)<-expr$Ensembl_ID
expr<-expr[,-1]
colnames(expr)<-gsub('.','-',colnames(expr),fixed = T) 
pheno=getTCGAdata(project = 'BRCA',download = TRUE,forceDownload = TRUE)
clinical=XenaPrepare(pheno) 
phe <- clinical[,c(54,8,10,22)]
#TNBC(n=107)
phe_tnbc <- phe %>% filter(!is.na(phe$bcr_sample_barcode) & 
                                  !is.na(phe$ER_Status_nature2012) &
                                  !is.na(phe$HER2_Final_Status_nature2012) &
                                  !is.na(phe$PR_Status_nature2012)) %>%
  filter((ER_Status_nature2012 == "Negative" & 
            HER2_Final_Status_nature2012 == "Negative" & 
            PR_Status_nature2012 == "Negative"))
phe_tpbc <- phe %>% filter(!is.na(phe$bcr_sample_barcode) & 
                                  !is.na(phe$ER_Status_nature2012) &
                                  !is.na(phe$HER2_Final_Status_nature2012) &
                                  !is.na(phe$PR_Status_nature2012)) %>%
  filter(!(ER_Status_nature2012 == "Negative" & 
            HER2_Final_Status_nature2012 == "Negative" & 
            PR_Status_nature2012 == "Negative"))
#TPBC (n=57)
phe_tpbc2 <- phe %>% filter(!is.na(phe$bcr_sample_barcode) & 
                             !is.na(phe$ER_Status_nature2012) &
                             !is.na(phe$HER2_Final_Status_nature2012) &
                             !is.na(phe$PR_Status_nature2012)) %>%
  filter((ER_Status_nature2012 == "Positive" & 
             HER2_Final_Status_nature2012 == "Positive" & 
             PR_Status_nature2012 == "Positive")) 
sample <- rbind(phe_tnbc,phe_tpbc2)
exprs <- t(expr)
write.csv(exprs,"exprs.csv")
exprs <- read.csv("exprs.csv")
exprs <- rename(exprs,c(X = "bcr_sample_barcode"))
count <- merge(sample,exprs,by="bcr_sample_barcode") #sample(n=161)
count <- count[c(order(count$ER_Status_nature2012)),]
count <- t(count)
write.csv(count,"count.csv")
count <- read.csv("count.csv")
count <- rename(count,c(X = "Gene.stable.ID.version"))
ID <- read.table("ID.txt",sep=",",header=TRUE,stringsAsFactors=FALSE,fill = TRUE,quote="")
HDcount <- merge(count,ID,by="Gene.stable.ID.version")
HDcount <- HDcount[,c(163,2:162)]
rownames(HDcount) <- make.names(HDcount[,1],unique = TRUE)
HDcount <- HDcount[,-1]
write.csv(HDcount,"HDcount.csv")
group <- factor(c(rep("TNBC",104),rep("TPBC",57)),levels=c("TNBC","TPBC"))
type <- c(rep("paired-end",161))
sampleID <- data.frame(sample=colnames(HDcount))
sampleID$sample <-gsub('.','-',sampleID$sample,fixed = T) 
sampleID$sample <- str_sub(sampleID$sample,1,12)

##############################################################################
# DNA methylation data (450K)
meth <- getTCGAdata(project = 'BRCA',download = T,Methylation = T,MethylationType = c("450K")) 
meth = XenaPrepare(meth) 
meth <- meth$HumanMethylation450.gz
write.csv(meth,"meth.csv")
meth <- read.table("meth.csv",sep = ",")
meth <- rename(meth,c(X = "sample"))
meth <- meth[,-1]
tmeth <- t(meth)
save(list = ls(),file="lab0610.RData") 
tmeth2 <- tmeth[-1,]
tmeth3 <- as.data.frame(tmeth2)
colnames(tmeth2) <- tmeth2[1,]
tmeth3 <- tmeth2[-1,]
tmeth4 <- as.data.frame(tmeth3)
tmeth4$sample <- str_sub(tmeth4$sample,1,12)
meth_merge <- merge(sampleID,tmeth4,by="sample")
meth_merge2 <- t(meth_merge)
colnames(meth_merge2) <- meth_merge2[1,]
meth_merge3 <- meth_merge2[-1,]
sampleID2 <- data.frame(sample=colnames(meth_merge3))
write.csv(meth_merge3,"meth_merge.csv")
save(meth_merge3,file = "./meth.Rdata")
HDcount2 <- t(HDcount)
HDcount2 <- as.data.frame(HDcount2)
write.csv(HDcount2,"HDcount2.csv")
HDcount2 <- read.csv("HDcount2.csv")
HDcount2 <- rename(HDcount2,c(X = "sample"))
HDcount2$sample <-gsub('.','-',HDcount2$sample,fixed = T) 
HDcount2$sample <- str_sub(HDcount2$sample,1,12)
sampleID_TNBC <- data.frame(sample=HDcount2$sample[1:104])
sampleID_TPBC <- data.frame(sample=HDcount2$sample[105:161])
sample_TNBC <- merge(sampleID_TNBC,sampleID2,by="sample")
sample_TPBC <- merge(sampleID_TPBC,sampleID2,by="sample")
HDcount_TNBC <- merge(sample_TNBC,HDcount2)
HDcount_TPBC <- merge(sample_TPBC,HDcount2)
HDcount <- rbind(HDcount_TNBC,HDcount_TPBC)
HDcount <- unique(HDcount)
HDcount <- t(HDcount)
HDcount <- as.data.frame(HDcount)
colnames(HDcount) <- HDcount[1,]
HDcount <- HDcount[-1,]
write.csv(HDcount,"HDcount.csv")
HDcount2 <- read.csv("HDcount.csv")
rownames(HDcount2) <- HDcount2[,1]
HDcount3 <- HDcount2[,-1]
#htseq-count : DEG - DESeq2
group <- factor(c(rep("TNBC",62),rep("TPBC",29)),levels=c("TNBC","TPBC"))
type <- c(rep("paired-end",91))
coldata <- data.frame(row.names = colnames(HDcount3),group,type)
HDcount <- 2^HDcount - 1  # log2(x+1) to interger 
HDcount <- ceiling(HDcount) 
countDataMatrix <- as.matrix(HDcount)
dds <- DESeqDataSetFromMatrix(countData = HDcount3,colData = coldata,design = ~group)
keep <- rowSums(counts(dds)) >= 20 # count first filtering 
dds <- dds[keep,]
dds <- DESeq(dds) #estimating size factors , dispersions, model and testing
res <- results(dds,contrast = c("group","TNBC","TPBC")) #result (TNBC vs TPBC)
resultsNames(dds) 
summary(res)
vsd <- vst(dds, blind=TRUE)
png('./DEG_DESeq2/SDplot.png')
meanSdPlot(assay(vsd))
dev.off()
pcaData <- plotPCA(vsd,intgroup=c("group","type"),returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))
png('./DEG_DESeq2/PCA.png')
ggplot(pcaData,aes(PC1,PC2,color=group.1,shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()
save(HDcount3,file = "./HDcount3.Rdata")
table(res$padj<0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
##ggplot2 : volcano plot 
deseq2 <- as.data.frame(res)
write.csv(deseq2,file="./DEG_DESeq2/deseq2.csv")
deseq2 <- read.csv("./DEG_DESeq2/deseq2.csv")
#padj(0.05) & log2FC(2) group 
for (i in 1:nrow(deseq2)) {
  if (abs(deseq2[i,'log2FoldChange']) >= 1) deseq2[i,'select_change'] <- 'y' else deseq2[i,'select_change'] <- 'n'
  if (deseq2[i,'padj'] %in% NA | abs(deseq2[i,'padj']) >= 0.05) deseq2[i,'select_pvalue'] <- 'n' else deseq2[i,'select_pvalue'] <- 'y'
  deseq2[i,'select'] <- paste(deseq2[i,'select_change'], deseq2[i,'select_pvalue'], sep = '')
}
deseq2 <- deseq2 %>% reshape::rename(c(X = 'Gene.stable.ID')) 
deseq2$select <- factor(deseq2$select,levels = c('nn','ny','yn','yy'),
                        labels = c('p >= 0.05, FC <  2', 'p < 0.05, FC < 2', 'p >= 0.05, FC >= 2', 'p < 0.05, FC >= 2')) 
volcano_plot_pvalue2 <- ggplot(deseq2, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = select), alpha = 0.6,show.legend =TRUE) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value')
volcano_plot_abundance2 <- ggplot(deseq2, aes(log2FoldChange, 100 * baseMean / sum(baseMean))) +
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position=c(0.15,0.9),legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  labs(x = 'log2 Fold Change', y = 'Abundance (%)')
png('./htseq-count/DESeq2/volcano_plot.png', width = 3000, height = 1600, res = 300, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(volcano_plot_pvalue2, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(volcano_plot_abundance2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
ggsave('./DEG_DESeq2/volcano_plot_pvalue2.png', volcano_plot_pvalue2, width = 6, height = 5)
ggsave('./htseq-count/DESeq2/volcano_plot_abundance2.png', volcano_plot_abundance2, width = 6, height = 5)
##csv file : DEG - padj<0.001 & abs(log2FoldChange) > 2 , order by padj 
diff_gene_deseq2 <-subset(deseq2, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file = "./DEG_DESeq2/DESeq2_htseqcounts.csv")
data2 <- read.csv("./DEG_DESeq2/DESeq2_htseqcounts.csv")
data2 <- rename(data2,c(X = "genes"))
data2 <- as.data.frame(data2[order(data2$padj),]) #padj order
write.csv(data2,file = "./DEG_DESeq2/DESeq2_htseqcounts.csv",row.names = F)
-----------------------------------------------------------------------------------------
#pheatmap
rld <- vst(dds, blind=F)
#rld <- rlogTransformation(dds,blind = F) #느림 
write.csv(assay(rld),file = "./DEG_DESeq2/counts.csv")
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),30)
mat  <- assay(rld)[ topVarGene, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[,c("group","sizeFactor")]) 
pheatmap(mat, annotation_col = anno)
#clusterprofile -GO
rownames(data2) <- data2[,1]
data3 <- data2[,-1]
gene <- rownames(diff_gene_deseq2)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
png('./DEG_DESeq2/GO.png')
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
dev.off()

-------------------------------------------------------------------------------
#KEGG
EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
png('./DEG_DESeq2/KEGG.png')
dotplot(EGG)
dev.off()
test <- data.frame(EGG)
browseKEGG(EGG, 'hsa04015') #ex:apoptosis PATH_id search

##############################################################################
