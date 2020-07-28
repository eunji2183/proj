#DEG - featurecounts/htseqcount : edgeR/DESeq2
rm(list = ls())
working_dir="/home/eunji/R/project/200717_microdust/"
setwd(working_dir)
options(stringsAsFactors = F) #no character2factor
save(list = ls(),file="microdust.RData")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("vsn")
BiocManager::install("hexbin")
BiocManager::install("apeglm")
BiocManager::install("clusterProfiler")
BiocManager::install("limma")
BiocManager::install("VennDiagram")
devtools::install_github("ToledoEM/msigdf")
BiocManager::install("qusage")
BiocManager::install("GSA")
BiocManager::install("ReactomePA")
BiocManager::install("ggupset")
library(DESeq2)
library(edgeR)
library(dplyr)
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
library(AnnotationDbi)
library(limma)
library(VennDiagram)
library(msigdf)
library(org.Hs.eg.db)
library(phantasus) 
library(rtracklayer)
library(qusage)
library(GSA) #GSA.read.gmt 
library(enrichplot)
library(ReactomePA)
library(stringr)
library(ggupset)

################################################################################
#featurecount-rawdata 
setwd("/home/eunji/R/project/200717_microdust/")
path <- "./data/featurecount"
filenames <- dir(path)
filepath <- sapply(filenames,function(x){
  paste(path,x,sep = '/')})
data <- lapply(filepath,function(x){
  read.table(x,header = T,sep = "\t",stringsAsFactors = F)})
names(data) <- str_sub(names(data),1,-5)

for(i in 1:length(data)){
  data[[i]] <- data[[i]][,c(1,7)] %>%
    reshape::rename(c(Geneid="gene_id")) 
  names(data[[i]])[2]=names(data)[[i]]
}

count <- Reduce(function(x,y)merge(x,y,by="gene_id"),data)
count <- count[,c(1,2,3,4,6,5,7)]

gtf <- rtracklayer::import("/home/eunji/project/0_sh/ref/rna/Homo_sapiens.GRCh38.96.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()


save(count,file = "./data/count.RData")
save(ID,file = "./data/ID.RData")

group <- factor(c(rep("con",3),rep("upm",3)),levels=c("con","upm"))

#=======================================================
#DESeq2 
rownames(count) <- count[,1]
FDcount <- count[,-1]
coldata <- data.frame(row.names = colnames(FDcount),group)
dds <- DESeqDataSetFromMatrix(countData = FDcount,colData = coldata,design = ~group)
dds <- DESeq(dds)
res <- results(dds,contrast = c("group","upm","con")) 
resultsNames(dds) 
summary(res)
table(res$padj<0.05)
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
res <- data.frame(gene_id=rownames(res),res)
rownames(res) <- NULL
rank <- merge(ID,res,by="gene_id") %>%
  dplyr::select(gene_name,log2FoldChange)
write.table(rank,file = "./result/FDrank.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1] <- 'gene_id'
resdata <- merge(ID,resdata,by="gene_id")
resdata <- na.omit(resdata)

write.csv(resdata,file="./result/DESeq2/DESeq2_featurecounts.csv")




gene_des <- read.csv("./data/id.txt",sep = "\t" ,header = T,stringsAsFactors = F,fill = T)
gene_des <- rename(gene_des,c(Gene.stable.ID = "Name",Gene.description = "DESCRIPTION"))
ID3 <- merge(ID,gene_des,by="Name")[,c(1,2,4)]
rawdata <- merge(x=count,y=ID3,by="Name")
rawdata <- rawdata[,c(8,9,5:7,2:4)]
rownames(rawdata) <- make.names(rawdata[,1],unique = TRUE)
count <- rawdata[,-1]
rm(CON,UPM,YJ_con,YJ_upm,p5_con,p5_upm,gtf,gene_des)
group <- factor(c(rep("upm",3),rep("con",3)),levels=c("upm","con"))
type=c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
##############################################################################################
#prepare GSEA file
rownames(rawdata) <- NULL
group <- paste(group,collapse = " ")
group <- c(paste(c(3+3,2,1),collapse = " "), "# upm con",group)
write.table(file = "./data/FDGSEA.gct",rawdata,sep = "\t",col.names = T,row.names = F,quote = F)
write.table(file = "group.cls",group,col.names = F,row.names = F,quote = F)


#################################################################################
#featurecount - edgeR
FEcount <- count[,1:6]
FEcount <- FEcount[rowSums(cpm(FEcount) > 1) >=2,]#remove low expression gene
dim(FEcount)
genes=rownames(FEcount)
dgelist <- DGEList(counts = FEcount, genes = genes, group = group)
dgelist_norm <- calcNormFactors(dgelist,method = 'TMM') #TMM normalization
plotMDS(dgelist_norm,col = rep(c('red', 'blue','black'), each = 2), dim = c(1, 2)) #plotMDS
design <- model.matrix(~group)
dge <- estimateDisp(dgelist_norm,design,robust = TRUE)
plotBCV(dge)
cpm=cpm(dgelist)
lcpm=cpm(dgelist,log=TRUE)
et <- exactTest(dge)
tTags <- topTags(et,n=nrow(dgelist$counts))
tTags <- as.data.frame(tTags)
tTag <- merge(x=tTags,y=ID2,by="genes")
tTag <- unique(tTag) #remove duplicate 
tTag <- tTag[c(order(tTag$FDR)),]
rownames(tTag) <- NULL
write.csv(tTag,"./featurecounts/edgeR/edgeR_featurecounts.csv")
de <- decideTestsDGE(et,adjust.method = 'fdr',p=.05)
plotMD(et, status = de, values = c(1, -1), col = c('red', 'blue'))
abline(h = c(-1, 1), col = 'gray', lty = 2)
for (i in 1:nrow(tTag)) {
  if (abs(tTag[i,'logFC']) >= 1) tTag[i,'select_change'] <- 'y' else tTag[i,'select_change'] <- 'n'
  if (tTag[i,'FDR'] %in% NA | abs(tTag[i,'FDR']) >= 0.05) tTag[i,'select_pvalue'] <- 'n' else tTag[i,'select_pvalue'] <- 'y'
  tTag[i,'select'] <- paste(tTag[i,'select_change'], tTag[i,'select_pvalue'], sep = '')
}
tTag$select <- factor(tTag$select,levels = c('nn','ny','yn','yy'),
                      labels = c('p >= 0.05, FC <  2', 'p < 0.05, FC < 2', 'p >= 0.05, FC >= 2', 'p < 0.05, FC >= 2')) 
volcano_plot_pvalue <- ggplot(tTag, aes(logFC, -log(FDR, 10))) +
  geom_point(aes(color = select), alpha = 0.6,show.legend = TRUE) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value')
volcano_plot_abundance <- ggplot(tTag, aes(logFC, 100 * logCPM / sum(logCPM))) +
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position=c(0.15,0.9),legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  labs(x = 'log2 Fold Change', y = 'Abundance (%)')
png('./featurecounts/edgeR/volcano_plot.png', width = 3000, height = 1600, res = 300, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(volcano_plot_pvalue, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(volcano_plot_abundance, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
ggsave('./featurecounts/edgeR/volcano_plot_pvalue.png', volcano_plot_pvalue, width = 6, height = 5)
ggsave('./featurecounts/edgeR/volcano_plot_abundance.png', volcano_plot_abundance, width = 6, height = 5)
-----------------------------------------------------------------------------------------------------------------
#GO : featurecounts- edgeR  
diff_gene_edgeR <-subset(tTag, FDR < 0.05 & abs(logFC) > 1)
rownames(diff_gene_edgeR) <- make.names(diff_gene_edgeR[,1],unique = TRUE)
diff_gene_edgeR<- diff_gene_edgeR[,-1]
gene <- rownames(diff_gene_edgeR)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
png('./featurecounts/edgeR/GO.png')
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
dev.off()

##KEGG 
EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
png('./featurecounts/edgeR/KEGG.png')
dotplot(EGG)
dev.off()
test <- data.frame(EGG)
browseKEGG(EGG, 'hsa04015') #ex:apoptosis PATH_id search

##################################################################################
#featurecount - DESeq2
FDcount <- count[,2:7]
coldata <- data.frame(row.names = colnames(FDcount),group,type)
FDcount <- as.data.frame(FDcount)
dds <- DESeqDataSetFromMatrix(countData = FDcount,colData = coldata,design = ~group)
keep <- rowSums(counts(dds)) >= 20 # count first filtering 
dds <- dds[keep,]
dds <- DESeq(dds) #estimating size factors , dispersions, model and testing
res <- results(dds,contrast = c("group","upm","con")) 
resultsNames(dds) 
summary(res)
table(res$padj<0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file="./featurecounts/con_vs_upm/DESeq2/DESeq2_featurecounts.csv")
FDranked <- resdata[,c(1,3)]
colnames(FDranked) <- FDranked[1,]
FDranked <- FDranked[-1,]
write.table(FDranked,"./data/FDranked.rnk",sep = "\t",row.names = F,quote = F)

normcount <- resdata[,c(1,8:13)]
normcount <- rename(normcount,c(Row.names="Name"))
gene_des <- ID3[,2:3] 
gene_des <- rename(gene_des,c(gene_name = "Name"))
FDGSEA <- merge(gene_des,normcount,by="Name") %>% unique()
write.table(file = "./data/FDGSEA.gct",FDGSEA,sep = "\t",col.names = T,row.names = F,quote = F)

#PCA 
vsd <- vst(dds, blind=TRUE)
png('./featurecounts/DESeq2/SDplot.png')
meanSdPlot(assay(vsd))
dev.off()
pcaData <- plotPCA(vsd,intgroup=c("group","type"),returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))
png('./featurecounts/DESeq2/PCA.png')
ggplot(pcaData,aes(PC1,PC2,color=group.1,shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()
-----------------------------------------------------------------------------------------
##ggplot2 : volcano plot 
deseq2 <- as.data.frame(res)
write.csv(deseq2,file="./featurecounts/DESeq2/deseq2.csv")
deseq2 <- read.csv("./featurecounts/DESeq2/deseq2.csv")
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
  geom_point(aes(color = select), alpha = 0.6,show.legend = TRUE) +
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
png('./featurecounts/DESeq2/volcano_plot.png', width = 3000, height = 1600, res = 300, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(volcano_plot_pvalue2, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(volcano_plot_abundance2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
ggsave('./featurecounts/DESeq2/volcano_plot_pvalue.png', volcano_plot_pvalue2, width = 6, height = 5)
ggsave('./featurecounts/DESeq2/volcano_plot_abundance.png', volcano_plot_abundance2, width = 6, height = 5)
##csv file : DEG - padj<0.001 & abs(log2FoldChange) > 2 , order by padj 
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file = "./featurecounts/DESeq2/DESeq2_featurecounts.csv")
data2 <- read.csv("./featurecounts/DESeq2/DESeq2_featurecounts.csv")
data2 <- rename(data2,c(X = "genes"))
result <- merge(data2,ID2,by="genes")
result <- as.data.frame(result[order(result$padj),])#padj order
result <- unique(result)
write.csv(result,file = "./featurecounts/DESeq2/DESeq2_featurecounts.csv",row.names = F)
----------------------------------------------------------------------------------------------
#plotMA
plotMA(res,ylim=c(-20,20))
gene <- rownames(res)
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#pheatmap
rld <- rlogTransformation(dds,blind = F)
write.csv(assay(rld),file = "./featurecounts/DESeq2/counts.csv")
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),30)
mat  <- assay(rld)[ topVarGene, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[,c("group","sizeFactor")]) 
png('./featurecounts/DESeq2/pheatmap.png')
pheatmap(mat, annotation_col = anno)
dev.off()

#clusterprofile -GO
gene <- rownames(diff_gene_deseq2)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
png('./featurecounts/DESeq2/GO.png')
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
dev.off()


#KEGG
EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
png('./featurecounts/DESeq2/KEGG.png')
dotplot(EGG)
dev.off()
test <- data.frame(EGG)
browseKEGG(EGG, 'hsa04015') #ex:apoptosis PATH_id search

-------------------------------------------------------------------------------------
#GSEA 
deseq2 <- as.data.frame(res)
deseq2 <- rename(deseq2,c(log2FoldChange = "logFC"))

#gene-prerank for GSEA
grank <- data.frame(SYMBOL=rownames(deseq2),logFC=deseq2[,2])
grank <- na.omit(grank[order(grank$logFC,decreasing = T),])
grank <- grank %>% dplyr::mutate(rank = rank(logFC, ties.method = "random"))
FDgene <- grank[,1:2]
write.table(FDgene,file = "./data/FDgrank.rnk",sep = "\t",row.names = F,quote = F,col.names = F)
rownames(grank) <- grank[,1]
grank <- grank[,-1]
gene_list <- grank$logFC
names(gene_list) <- as.character(rownames(grank))


#ENTREZID-prerank for GSEA 
nrDEG <- deseq2[,c(2,6)]
eg <- bitr(geneID = rownames(nrDEG),fromType = "SYMBOL",
           toType = "ENTREZID",OrgDb = org.Hs.eg.db)

nrDEG$ENTREZID <- eg[match(rownames(nrDEG),eg$SYMBOL),2]
nrDEG2 <- nrDEG[,c(1,3)]
rownames(nrDEG2) <- NULL
nrDEG2 <- na.omit(nrDEG2[order(nrDEG2$logFC,decreasing = T),])
nrDEG2 <- nrDEG2 %>% dplyr::mutate(rank = rank(logFC, ties.method = "random"))
nrDEG2 <- nrDEG2[,c(2,1,3)]
rownames(nrDEG2) <- nrDEG2[,1]
nrDEG2 <- nrDEG2[,-1]
gene_list <- nrDEG2$logFC
names(gene_list) <- as.character(rownames(nrDEG2))
FDrank <- data.frame(ENTREZID=rownames(nrDEG2),logFC=nrDEG2[,1])
write.table(FDrank,file = "./data/FDrank.rnk",sep = "\t",row.names = F,quote = F,col.names = F)


nrDEG <- na.omit(nrDEG[order(nrDEG$logFC,decreasing = T),])  #sort by FC
nrDEG2 <- nrDEG %>% dplyr::mutate(rank = rank(logFC, ties.method = "random"))
nrDEG <- data.frame(symbol=rownames(nrDEG),nrDEG)
rownames(nrDEG) <- NULL
nrDEG3 <- merge(nrDEG,nrDEG2,by="ENTREZID")
nrDEG3 <- nrDEG3[,c(2,1,3,4,7)]
nrDEG3 <- rename(nrDEG3,c(logFC.x = "logFC",padj.x = "padj"))
rownames(nrDEG3) <- nrDEG3[,1]
nrDEG3 <- nrDEG3[,-1]
nrDEG <- na.omit(nrDEG3[order(nrDEG3$rank,decreasing = T),])
gene_list <- nrDEG$logFC
names(gene_list) <- as.character(rownames(nrDEG))

ridgeplot(gseGO.res,5)
gseaplot2(gseGO.res,1:5)
gseaplot2(H_GSEA, 1, title = "", color = "green", base_size = 11,
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = T,
          ES_geom = "line")

S <- GSA.read.gmt("data/msigdb.v7.1.symbols.gmt")
len_vec=c()
len_vec[1]=3
for(i in 1:length(S$genesets)){len_vec[i] <- c(length(S$genesets[[i]]))}
pathway_vec <- unlist(Vectorize(rep.int)(S$geneset.names,len_vec),use.names = F)
S7_1 <- as.data.frame(cbind(geneset=pathway_vec,symbol=unlist(S$genesets,use.names = F)))
S_GSEA <- GSEA(gene_list,exponent = 1,TERM2GENE = S7_1,TERM2NAME = NA,nPerm = 1000,minGSSize = 15,maxGSSize = 500,pvalueCutoff = 0.25,pAdjustMethod = "BH",verbose = T,seed = F)
S_results <- as.data.frame(S_GSEA@result)
write.csv(S_results,file = "./result/gsea.csv")

gene <- names(gene_list)
egmt <- enricher(gene,TERM2GENE = S7_1) #hypergeometric test
write.csv(egmt,file = "./result/egmt.csv")

gseaplot2(S_GSEA,geneSetID=1,title=S_GSEA$Description[1],pvalue_table=T)
barplot(egmt,showCategory=20)
dotplot(egmt,showCategory=30)
egmt2 <-simplify(egmt)
cnetplot(egmt2,foldChange=gene_list)
upsetplot(egmt)
heatplot(egmt,foldChange=gene_list)
emapplot(egmt)

H <- msigdf.human %>%
  dplyr::filter(category_code == "hallmark") %>% dplyr::select(geneset,symbol) %>% as.data.frame
H_GSEA <- GSEA(gene_list,TERM2GENE = H,nPerm = 1000,minGSSize = 15,maxGSSize = 500,pvalueCutoff = 0.05)
H_results <- as.data.frame(H_GSEA@result)

write.csv(H_results,"./featurecounts/con_vs_upm/DESeq2/H_result.csv")
write.csv(c1_results,"./featurecounts/con_vs_upm/DESeq2/c1_result.csv")



##################################################################################################################
#htseq-count : rawdata
rm(list = ls())
working_dir="/home/eunji/project/rna/project/5.DEG/"
setwd(working_dir)
count <- read.table("./htseq-count/count.txt",stringsAsFactors=FALSE)
count <- rename(count,c(V1 = "Gene.stable.ID",V2 = "CON",
                        V3 = "UPM", V4 = "p5_con", V5 = "p5_upm",
                        V6 = "YJ_con", V7 = "YJ_upm"))
#ID : http://asia.ensembl.org/index.html 
ID <- read.table("./featurecounts/ID.txt",sep=",",header=TRUE,stringsAsFactors=FALSE,fill = TRUE,quote="")
rawdata <- merge(x=count,y=ID,by="Gene.stable.ID")
ID <- rename(ID,c(Gene.name = "genes"))
ID2 <- ID[,2:3] 
rawdata <- rawdata[,c(8,2:7,9)]
rownames(rawdata) <- make.names(rawdata[,1],unique = TRUE)
count <- rawdata[,-1]
rm(rawdata)
group <- factor(c(rep("con",3),rep("upm",3)),levels=c("con","upm"))
type=c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
###########################################################################################
#htseq-count :edgeR 
HEcount <- count[,c(1,5,2,6)]
HEcount <- HEcount[rowSums(cpm(HEcount) > 1) >=2,]#remove low expression gene
dim(HEcount)
genes=rownames(HEcount)
dgelist <- DGEList(counts = HEcount, genes = genes, group = group)
dgelist_norm <- calcNormFactors(dgelist,method = 'TMM') #TMM normalization
plotMDS(dgelist_norm,col = rep(c('red', 'blue','black'), each = 2), dim = c(1, 2)) #plotMDS
design <- model.matrix(~group)
dge <- estimateDisp(dgelist_norm,design,robust = TRUE)
plotBCV(dge)
cpm=cpm(dgelist)
lcpm=cpm(dgelist,log=TRUE)
et <- exactTest(dge)
tTags <- topTags(et,n=nrow(dgelist$counts))
tTags <- as.data.frame(tTags)
tTag <- merge(x=tTags,y=ID2,by="genes")
tTag <- unique(tTag) #remove duplicate 
tTag <- tTag[c(order(tTag$FDR)),]
rownames(tTag) <- NULL
write.csv(tTag,"./htseq-count/edgeR/edgeR_htseqcounts.csv")
de <- decideTestsDGE(et,adjust.method = 'fdr',p=.05)
plotMD(et, status = de, values = c(1, -1), col = c('red', 'blue'))
abline(h = c(-1, 1), col = 'gray', lty = 2)
for (i in 1:nrow(tTag)) {
  if (abs(tTag[i,'logFC']) >= 1) tTag[i,'select_change'] <- 'y' else tTag[i,'select_change'] <- 'n'
  if (tTag[i,'FDR'] %in% NA | abs(tTag[i,'FDR']) >= 0.05) tTag[i,'select_pvalue'] <- 'n' else tTag[i,'select_pvalue'] <- 'y'
  tTag[i,'select'] <- paste(tTag[i,'select_change'], tTag[i,'select_pvalue'], sep = '')
}
tTag$select <- factor(tTag$select,levels = c('nn','ny','yn','yy'),
                      labels = c('p >= 0.05, FC <  2', 'p < 0.05, FC < 2', 'p >= 0.05, FC >= 2', 'p < 0.05, FC >= 2')) 
volcano_plot_pvalue <- ggplot(tTag, aes(logFC, -log(FDR, 10))) +
  geom_point(aes(color = select), alpha = 0.6,show.legend = TRUE) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value')
volcano_plot_abundance <- ggplot(tTag, aes(logFC, 100 * logCPM / sum(logCPM))) +
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position=c(0.15,0.9),legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  labs(x = 'log2 Fold Change', y = 'Abundance (%)')
png('./htseq-count/edgeR/volcano_plot.png', width = 3000, height = 1600, res = 300, units = 'px')
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
print(volcano_plot_pvalue, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(volcano_plot_abundance, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
ggsave('./htseq-count/edgeR/volcano_plot_pvalue.png', volcano_plot_pvalue, width = 6, height = 5)
ggsave('./htseq-count/edgeR/volcano_plot_abundance.png', volcano_plot_abundance, width = 6, height = 5)
-----------------------------------------------------------------------------------------------------------------
#GO : htseq-count- edgeR  
diff_gene_edgeR <-subset(tTag, FDR < 0.05 & abs(logFC) > 1)
rownames(diff_gene_edgeR) <- make.names(diff_gene_edgeR[,1],unique = TRUE)
diff_gene_edgeR<- diff_gene_edgeR[,-1]
gene <- rownames(diff_gene_edgeR)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
png('./htseq-count/edgeR/GO.png')
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
dev.off()
##KEGG 
EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
png('./htseq-count/edgeR/KEGG.png')
dotplot(EGG)
dev.off()
test <- data.frame(EGG)
browseKEGG(EGG, 'hsa04015') #ex:apoptosis PATH_id search

###########################################################################################################
#htseq-count : DESeq2
HDcount <- count[,c(1,5,2,6)]
coldata <- data.frame(row.names = colnames(HDcount),group,type)
HDcount <- as.data.frame(HDcount)
dds <- DESeqDataSetFromMatrix(countData = HDcount,colData = coldata,design = ~group)
keep <- rowSums(counts(dds)) >= 20 # count first filtering 
dds <- dds[keep,]
dds <- DESeq(dds) #estimating size factors , dispersions, model and testing
res <- results(dds,contrast = c("group","con","upm")) 
resultsNames(dds) 
summary(res)
vsd <- vst(dds, blind=TRUE)
png('./htseq-count/DESeq2/SDplot.png')
meanSdPlot(assay(vsd))
dev.off()
pcaData <- plotPCA(vsd,intgroup=c("group","type"),returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))
png('./htseq-count/DESeq2/PCA.png')
ggplot(pcaData,aes(PC1,PC2,color=group.1,shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()
table(res$padj<0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
##ggplot2 : volcano plot 
deseq2 <- as.data.frame(res)
write.csv(deseq2,file="./htseq-count/DESeq2/deseq2.csv")
deseq2 <- read.csv("./htseq-count/DESeq2/deseq2.csv")
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
ggsave('./htseq-count/DESeq2/volcano_plot_pvalue.png', volcano_plot_pvalue2, width = 6, height = 5)
ggsave('./htseq-count/DESeq2/volcano_plot_abundance2.png', volcano_plot_abundance2, width = 6, height = 5)
##csv file : DEG - padj<0.001 & abs(log2FoldChange) > 2 , order by padj 
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file = "./htseq-count/DESeq2/DESeq2_htseqcounts.csv")
data2 <- read.csv("./htseq-count/DESeq2/DESeq2_htseqcounts.csv")
data2 <- rename(data2,c(X = "genes"))
result <- merge(data2,ID2,by="genes")
result <- as.data.frame(result[order(result$padj),]) #padj order
result <- unique(result)
write.csv(result,file = "./htseq-count/DESeq2/DESeq2_htseqcounts.csv",row.names = F)
-----------------------------------------------------------------------------------------
#pheatmap
rld <- rlogTransformation(dds,blind = F)
write.csv(assay(rld),file = "./htseq-count/DESeq2/counts.csv")
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),30)
mat  <- assay(rld)[ topVarGene, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[,c("group","sizeFactor")]) 
pheatmap(mat, annotation_col = anno)
#clusterprofile -GO
gene <- rownames(diff_gene_deseq2)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
png('./htseq-count/DESeq2/GO.png')
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
dev.off()

-------------------------------------------------------------------------------
#KEGG
EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
png('./htseq-count/DESeq2/KEGG.png')
dotplot(EGG)
dev.off()
test <- data.frame(EGG)
browseKEGG(EGG, 'hsa04015') #ex:apoptosis PATH_id search

####################################################################################

