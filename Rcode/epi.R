#epidemiology paper 
#BRCA - DEG - GO & KEGG -survival 

rm(list = ls())
working_dir="/home/eunji/R/project/200625_epi/"
setwd(working_dir)
options(stringsAsFactors = F)
save(list = ls(),file="./data/200625_epi.RData")
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
devtools::install_github("ToledoEM/msigdf")
install.packages("protobuf-compiler")
BiocManager::install("protolite")
BiocManager::install("phantasus")
library(enrichplot)
library(magrittr) # %>%
library(dplyr) #filter
library(data.table) #fread
library(DESeq2)
library(edgeR)
library(ggplot2)
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
library(UCSCXenaTools)
library(estimate)
library(org.Hs.eg.db)
library(msigdf)
library(protolite)
library(phantasus) # write.gct
library(caret)
library(table1)
library(msigdf) #GSEA msigdb
library(stringr)
####################################################################################
#TCGAGTEx : Tumor vs normal gene expression data (BRCA)
#https://toil.xenahubs.net/download/TcgaTargetGtex_gene_expected_count.gz
#log2(expected_count+1), gencode.v23.annotation.gtf.gz
#phenotypes, primary sites : breast 
pheno <- read.table("./data/TcgaTargetGTEX_phenotype.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,fill = TRUE,quote="")
pheno <- pheno %>%
  dplyr::filter(pheno$X_primary_site == 'Breast') %>%
  dplyr::filter(pheno$X_gender == 'Female')
tumor <- pheno[grep("0.$",pheno$sample),]
normal <- pheno[-grep("0.$",pheno$sample),]
tumor_sample <- c("sample",tumor$sample)
normal_sample <- c("sample",normal$sample[1:191])
#gene expression data
system.time(genedata <- fread("/home/eunji/R/project/200625_epi/data/TcgaTargetGtex_gene_expected_count",header = F,stringsAsFactors = F))
system.time(genedata <- as.data.frame(genedata))
system.time(tgenedata <- t(genedata))
colnames(genedata) <- genedata[1,]
genedata <- genedata[-1,]
normal_merge <- genedata %>%
  dplyr::select(all_of(normal_sample))
tumor_merge <- genedata %>%
  dplyr::select(all_of(tumor_sample))
save(normal_merge,file = "./data/normal_merge.RData")
save(tumor_merge,file = "./data/tumor_merge.RData")
rownames(tumor_merge) <- tumor_merge[,1]
tumor_merge2 <- tumor_merge[,-1]
#############################################################################################
#BRCA clinical matrix
getTCGAdata(project = 'BRCA') 
pheno=getTCGAdata(project = 'BRCA',download = TRUE,forceDownload = TRUE)
clinical=XenaPrepare(pheno)
clin <- clinical %>%
  dplyr::select(sampleID,PAM50Call_RNAseq,OS_Time_nature2012,
                days_to_last_followup,ER_Status_nature2012,
                HER2_Final_Status_nature2012,PR_Status_nature2012,
                Vital_Status_nature2012,additional_radiation_therapy,
                vital_status,radiation_therapy,pathologic_stage,
                pathologic_T,pathologic_N,pathologic_M,menopause_status,
                gender,Age_at_Initial_Pathologic_Diagnosis_nature2012,Tumor_nature2012)
load("~/R/project/200625_epi/data/normal_merge.RData")
load("~/R/project/200625_epi/data/tumor_merge.RData")
normal <- t(normal_merge) %>% as.data.frame 
colnames(normal) <- normal[1,]
normal <- normal[-1,]
normal <- data.frame(sampleID=rownames(normal),normal) 
rownames(normal) <- NULL
norm_merge <- merge(clin,normal,by="sampleID")

tumor <- t(tumor_merge) %>% as.data.frame
colnames(tumor) <- tumor[1,]
tumor <- tumor[-1,]
tumor <- data.frame(sampleID=rownames(tumor),tumor)
rownames(tumor) <- NULL
tu_merge <- merge(clin,tumor,by="sampleID")
tuclin <- tu_merge[,c(1:7,9:19)]
tuclin <- tuclin %>%
  dplyr::mutate(age=ifelse(Age_at_Initial_Pathologic_Diagnosis_nature2012< 50 ,"<50",">=50"))

#training set vs validation set (2:1)
index <- createDataPartition(y=tuclin$age,p=0.667,list = F)
train <- tuclin[index,]
valid <- tuclin[-index,]
train <- train %>%
  dplyr::mutate(dataset= 1)
valid <- valid %>%
  dplyr::mutate(dataset= 2)
data <- rbind(train,valid)
temp <- as.matrix(data[,c(2,5:15,18,19)])
temp[is.na(temp)] <- "Unknown"
data[,c(2,5:15,18,19)] <- temp
data <- data %>%
  dplyr::mutate(pathologic_stage=
                  as.character(factor(pathologic_stage,
                                      levels = c("Stage I","Stage IA","Stage IB",
                                                 "Stage II","Stage IIA","Stage IIB", 
                                                 "Stage III","Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                 "Stage IV","Stage X","Unknown"),
                                      labels = c(1,1,1,2,2,2,3,3,3,3,4,0,0))))

################################################################################################
#table 
data$dataset <- 
  factor(data$dataset,
         levels = c(1,2),
         labels = c("Training dataset","Validation dataset"))
data$pathologic_stage <- 
  factor(data$pathologic_stage,
         levels = c(1,2,3,4,0),
         labels = c("I","II","III","IV","unknown"))
data$menopause_status <- 
  factor(data$menopause_status,
         levels = c("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)",
                    "Peri (6-12 months since last menstrual period)",
                    "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)",
                    "Indeterminate (neither Pre or Postmenopausal)","Unknown"),
         labels = c("premenopause","perimenopause","postmenopause",
                    "Indeterminate","Indeterminate"))

table1(~age + pathologic_stage + menopause_status | dataset,data=data)
########################################################
#TCGAGTEx : Tumor vs normal gene expression data (BRCA)
rownames(train) <- train[,1]
train <- train[,-1]
rownames(valid) <- valid[,1]
valid <- valid[,-1]
rawdata <- merge(normal_merge,tumor_merge,by="sample")
rownames(rawdata) <- rawdata[,1]
rawdata <- rawdata[,-1]
write.csv(rawdata,"./data/rawdata.csv")
rawdata <- read.csv("./data/rawdata.csv",stringsAsFactors = F,header = T)
rawdata <- 2^rawdata - 1
rawdata <- ceiling(rawdata)
save(rawdata,file = "./data/rawdata.RData")
id <- rownames(rawdata)
rawdata <- data.frame(id=id,rawdata)
rownames(rawdata) <- NULL
ID <- read.csv("./data/gencode.v23.annotation.gene.probemap",sep = "\t",stringsAsFactors = F,header = T)
ID2 <- ID[,2:6]
rawcount <- merge(ID,rawdata,by="id")
rownames(rawcount) <- NULL
save(rawcount,file = "./data/rawcount.RData")
count <- rawcount[,c(2:1283)]
rownames(count) <- make.names(count[,1],unique = TRUE)
count <- count[,-1]
colnames(count) <- gsub('.','-',colnames(count),fixed = T)
training <- count %>%
  dplyr::select(all_of(rownames(train)))
validation <- count %>%
  dplyr::select(all_of(rownames(valid)))
normal <- count[,5:195]
train_count <- cbind(normal,training)  
ID2 <- rawcount[,2:6]
##############################################################################################################
#train_count - DEG : DESeq2 
Dcount <- train_count[,1:917]
group <- factor(c(rep("normal",191),rep("tumor",726)),levels=c("normal","tumor"))
coldata <- data.frame(row.names = colnames(Dcount),group)
dds <- DESeqDataSetFromMatrix(countData = Dcount,colData = coldata,design = ~group)
keep <- rowSums(counts(dds)) >= 20 # count first filtering 
dds <- dds[keep,]
system.time(dds <- DESeq(dds)) #estimating size factors , dispersions, model and testing
save(dds,file = "./data/dds.RData")
res <- results(dds,contrast = c("group","normal","tumor"),independentFiltering = FALSE) 
summary(res)
DEG <- as.data.frame(res)
DEG <- na.omit(DEG)
nrDEG <- DEG[,c(2,6)]
table(res$padj<0.01)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file="./result/train_DESeq2.csv")

#####################################################################################################################
#GSEA 
eg <- bitr(geneID = rownames(nrDEG),fromType = "SYMBOL",
           toType = "ENTREZID",OrgDb = org.Hs.eg.db)

nrDEG$ENTREZID <- eg[match(rownames(nrDEG),eg$SYMBOL),2]
nrDEG <- na.omit(nrDEG[order(nrDEG$log2FoldChange,decreasing = T),])  #sort by FC
gene_list <- nrDEG$log2FoldChange
names(gene_list) <- as.character(rownames(nrDEG))
c1 <- msigdf.human %>%
  dplyr::filter(category_code == "c1") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c1_GSEA <- GSEA(gene_list,TERM2GENE = c1,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
c1_results <- as.data.frame(c1_GSEA@result)
c2 <- msigdf.human %>% 
  dplyr::filter(category_code == "c2") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c2_GSEA <- GSEA(gene_list,TERM2GENE = c2,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.02)
c2_results <- as.data.frame(c2_GSEA@result)
c3 <- msigdf.human %>%
  dplyr::filter(category_code == "c3") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c3_GSEA <- GSEA(gene_list,TERM2GENE = c3,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
c3_results <- as.data.frame(c3_GSEA@result)
c4 <- msigdf.human %>%
  dplyr::filter(category_code == "c4") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c4_GSEA <- GSEA(gene_list,TERM2GENE = c4,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
c4_results <- as.data.frame(c4_GSEA@result)
c5 <- msigdf.human %>%
  dplyr::filter(category_code == "c5") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c5_GSEA <- GSEA(gene_list,TERM2GENE = c5,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
c5_results <- as.data.frame(c5_GSEA@result)
c6 <- msigdf.human %>%
  dplyr::filter(category_code == "c6") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c6_GSEA <- GSEA(gene_list,TERM2GENE = c6,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
c6_results <- as.data.frame(c6_GSEA@result)
c7 <- msigdf.human %>%
  dplyr::filter(category_code == "c7") %>% dplyr::select(geneset, symbol) %>% as.data.frame
c7_GSEA <- GSEA(gene_list,TERM2GENE = c7,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
c7_results <- as.data.frame(c7_GSEA@result)
H <- msigdf.human %>%
  dplyr::filter(category_code == "hallmark") %>% dplyr::select(geneset,symbol) %>% as.data.frame
H_GSEA <- GSEA(gene_list,TERM2GENE = H,nPerm = 1000,minGSSize = 5,maxGSSize = 1000,pvalueCutoff = 0.05)
H_results <- as.data.frame(H_GSEA@result)
gseaplot2(H_GSEA, 1, title = "", color = "black", base_size = 11,
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = T,
          ES_geom = "line")
H_term <- H %>%
    dplyr::filter(geneset == "HALLMARK_GLYCOLYSIS")
write.csv(nrDEG,"./data/nrDEG.csv")
nrDEG2 <- read.csv("./data/nrDEG.csv")
nrDEG2 <- rename(nrDEG2,c(X="symbol"))
H1 <- merge(H_term,nrDEG2,by="symbol")
#gseKEGG
rownames(up) <- up[,1]
up <- up[,-1]
kup <- up[,c(2,6)]
kegg_up <- bitr(unique(row.names(up)),
                fromType = "SYMBOL",toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
kup$ENTREZID <- eg[match(rownames(kup),kegg_up$SYMBOL),2]
kup <- na.omit(kup[order(kup$log2FoldChange,decreasing = T),])  #sort by FC
kegg_list <- kup$log2FoldChange
names(kegg_list) <- as.character(kup$ENTREZID)
BRCA_KEGG <- gseKEGG(geneList = kegg_list,nPerm = 1000,
                     keyType = 'kegg',organism = 'hsa',
                     pvalueCutoff = 0.5,pAdjustMethod = 'BH') 
rownames(down) <- down[,1]
down <- down[,-1]
kdown <- down[,c(2,6)]
kegg_down <- bitr(unique(row.names(down)),
                fromType = "SYMBOL",toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
kdown$ENTREZID <- eg[match(rownames(kdown),kegg_down$SYMBOL),2]
kdown <- na.omit(kdown[order(kdown$log2FoldChange,decreasing = T),])  #sort by FC
kegg_list2 <- kdown$log2FoldChange
names(kegg_list2) <- as.character(kdown$ENTREZID)
BRCA_KEGG2 <- gseKEGG(geneList = kegg_list2,nPerm = 1000,
                     keyType = 'kegg',organism = 'hsa',
                     pvalueCutoff = 0.1,pAdjustMethod = 'BH') 
BRCA_KEGG2$Description
BRCA_KEGG2$enrichmentScore
kegg_down_results <- as.data.frame(BRCA_KEGG2@result)

#





#PCA 
vsd <- vst(dds, blind=TRUE)
png('./doc/SDplot.png')
meanSdPlot(assay(vsd))
dev.off()
pcaData <- plotPCA(vsd,intgroup=c("group"),returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))
png('./doc/PCA.png')
ggplot(pcaData,aes(PC1,PC2,color=group.1)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()

##ggplot2 : volcano plot 
deseq2 <- as.data.frame(res)
write.csv(deseq2,file="./doc/deseq2.csv")
deseq2 <- read.csv("./doc/deseq2.csv")
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
ggsave('./doc/volcano_plot_pvalue.png', volcano_plot_pvalue2, width = 6, height = 5)

##csv file : DEG - padj<0.001 & abs(log2FoldChange) > 2 , order by padj 
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file = "./doc/DESeq2.csv")
data2 <- read.csv("./doc/DESeq2.csv")
up <- data2 %>%
  dplyr::filter(log2FoldChange > 1)
down <- data2 %>%
  dplyr::filter(log2FoldChange < -1)
data2 <- rename(data2,c(X = "gene"))
results <- data2[order(data2$padj),]
save(results,file = "./data/results.RData")
result <- merge(ID2,data2,by="gene")
result <- as.data.frame(result[order(result$padj),])#padj order
write.csv(result,file = "./doc/DESeq2.csv",row.names = F)

#pheatmap
rld <- vst(dds,blind = F)
write.csv(assay(rld),file = "./doc/counts.csv")
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),30)
mat  <- assay(rld)[ topVarGene, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[,c("group","sizeFactor")]) 
png('./doc/pheatmap.png')
pheatmap(mat, annotation_col = anno)
dev.off()

#clusterprofile -GO
gene <- rownames(diff_gene_deseq2)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- as.vector(gene$ENTREZID)
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
png('./doc/GO.png')
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")
dev.off()


#KEGG
EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
png('./doc/KEGG.png')
dotplot(EGG)
dev.off()
test <- data.frame(EGG)
browseKEGG(EGG, 'hsa04062') #ex:apoptosis PATH_id search
cytokine <- test %>%
  + dplyr::filter(ID == "hsa04060")
cytokine_gene <- cytokine %>%
  + dplyr::select(geneID)
cytokgene <- data.frame(ID=strsplit(cytokine_gene$geneID,"/"))
cytokgene <- rename(cytokgene,c(c..3624....130399....6356....3953....3977....7048....3590....3627...="ENTREZID"))
cytokgene <- merge(gene,cytokgene,by="ENTREZID")
#GSEA
deseq2 <- read.csv("./doc/deseq2.csv",row.names = 1)
deseq2 <- rename(deseq2,c(log2FoldChange = "logFC"))
if(T){
  genelist <- deseq2$logFC
  names(genelist) <- rownames(deseq2)
  genelist_tr <- bitr(names(genelist),
                      fromType = "SYMBOL",
                      toType = c("ENSEMBL","ENTREZID"),
                      OrgDb = org.Hs.eg.db)
  new_list <- data.frame(SYMBOL=names(genelist),logFC = as.numeric(genelist))
  new_list <- merge(new_list,genelist_tr,by="SYMBOL")
  genelist <- new_list$logFC
  names(genelist) <- new_list$SYMBOL
  genelist <- sort(genelist,decreasing = T)}
msigdb <- msigdf.human
msigdb <- as.data.frame(msigdb)
c2 <- msigdb %>% 
  dplyr::filter(category_code == "c2") %>% dplyr::select(geneset, symbol) %>% as.data.frame
de <- names(genelist)[abs(genelist) > 1]
x <- enricher(genelist,TERM2GENE = c2)
y <- GSEA(genelist,TERM2GENE = msigdb)


##############################################################################
#TCGA BRCA methylation K450
meth <- fread("./data/TCGA-BRCA.methylation450.tsv",stringsAsFactors = F,header = T)
ID_meth <- fread("./data/illuminaMethyl450_hg38_GDC",stringsAsFactors = F)
meth <- as.data.frame(meth)
rownames(meth) <- meth[,1]
meth <- meth[,-1]
colnames(meth) <- str_sub(colnames(meth),1,15)
samples <- colnames(meth)
meth2 <- meth %>% 
  dplyr::select(contains(samples))
tu_meth <- subset(tumor_merge2,select=samples)
