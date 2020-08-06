#colon cancer - metastasis - GSEA 
rm(list = ls())
working_dir="/home/eunji/proj/200721_colon_M/"
setwd(working_dir)
options(stringsAsFactors = F)

library(UCSCXenaTools)
library(TCGAbiolinks)
library(SummarizedExperiment) #assay 
library(magrittr)
library(dplyr)
library(stringr) #str_sub
library(DESeq2)
library(edgeR)
library(limma)
library(rtracklayer)
library(ReactomePA)
library(enrichplot)
library(qusage) #read.gmt
library(clusterProfiler)



#======================================================
#data prepare  

#UCSC xena COADREAD clinical matrix 
getTCGAdata(project = 'COADREAD') 
pheno=getTCGAdata(project = 'COADREAD',download = TRUE,forceDownload = TRUE)
clinical=XenaPrepare(pheno) 


#gene count 
cancer <- c("TCGA-COAD","TCGA-READ")

for (i in 1:2) {
  cancer_select <- cancer[i]
  print(cancer)
  setwd("/home/eunji/proj/200721_colon_M/data/count/")
  suppressMessages({
    query <- GDCquery(
      project = cancer_select,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "HTSeq - Counts")  })
  if (is.null(query)){
    print(paste0("No count data for ", cancer_select ))
  } else{
    
    GDCdownload(query, method = "api", 
                files.per.chunk = 150)
    expdat <- GDCprepare(query = query)
    setwd("/home/eunji/proj/200721_colon_M/data/count/")
    #data-prepeocessing ( sprearman > 0.6)
    #expdat <- TCGAanalyze_Preprocessing(object = expdat,
    #cor.cut = 0.6,
    #filename = paste0(cancer_select,".png"),
    #datatype = "HTSeq - Counts")
    count_matrix=assay(expdat)
    write.csv(count_matrix,
              file = paste( cancer_select,"Counts.csv",
                            sep = "-"))}}

#======================================================
#TMN stage(1,2,3) vs stage4

TMN <- clinical %>% 
  dplyr::select(`_PATIENT`,pathologic_stage) %>% distinct()

names(TMN)[1] <- 'sampleID'

TMN <- TMN %>%
  dplyr::filter(!is.na(pathologic_stage)) %>%
  dplyr::mutate(group=ifelse(pathologic_stage %in% c("Stage IV","Stage IVA","Stage IVB"),'Y','N'))

TMN <- dplyr::arrange(TMN,desc(group)) 

TMN <- TMN %>%
  dplyr::select(sampleID,group)

save(TMN,file = "./result/TMN.RData")


#======================================================
#metastasis (no vs yes) 

clin <- clinical %>% 
  dplyr::select(`_PATIENT`,pathologic_M) %>% distinct()
  
names(clin)[1] <- 'sampleID'

clin <- clin %>%
  dplyr::filter(pathologic_M %in% c("M0","M1","M1a","M1b")) %>%
  dplyr::mutate(M_group=ifelse(pathologic_M=='M0','N',"Y")) %>%
  dplyr::select(sampleID,M_group)

clin <- dplyr::arrange(clin,desc(M_group))

save(clin,file = "./result/clin_metastasis.RData")

#=====================================================
# merge metastasis group & count 

COAD <- read.csv("./data/TCGA-COAD.htseq_counts.tsv",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
READ <- read.csv("./data/TCGA-READ.htseq_counts.tsv",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
count <- merge(COAD,READ,by="row.names")
count <- count[-c(1:5),]
rownames(count) <- count[,1]
count <- count[,-1]
colnames(count) <- str_sub(colnames(count),1,12)
rownames(count) <- str_sub(rownames(count),1,15)

count <- as.data.frame(t(count))
count <- data.frame(sampleID=rownames(count),count)
rownames(count) <- NULL


count$sampleID <- gsub('.','-',count$sampleID,fixed = T)

merge <- merge(clin,count,by="sampleID") 
#merge <- merge(TMN,count,by="sampleID")

rownames(merge) <- merge[,1]
merge <- merge[,-1]

merge <- dplyr::arrange(merge,desc(group))


coldata <- data.frame(row.names = rownames(merge),M_group=merge$M_group)
#coldata <- data.frame(row.names = rownames(merge),group=merge$group)

merge <- as.data.frame(t(merge))
merge <- merge[-1,]

# !!! data preprocessing (데이터정제)
merge <- na.omit(merge)
sapply(merge,mode) 
sapply(merge,class)
merge[,1:547] <- sapply(merge[,1:547],as.numeric)
#merge[,1:601] <- sapply(merge[,1:601],as.numeric)
meta <- 2^merge - 1
merge <- ceiling(meta)


save(coldata,file = "./result/coldata.RData")
save(merge,file = "./result/merge.RData")

#====================================================
#DEG : DESeq2 

#ID 
gtf <- rtracklayer::import("/home/eunji/proj/0_sh/ref/rna/Homo_sapiens.GRCh38.96.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()
save(ID,file = "./result/ID.RData")

load("~/proj/200721_colon_M/result/coldata.RData")
load("~/proj/200721_colon_M/result/merge.RData")
names(coldata)[1] <- 'group'
coldata$group <- factor(coldata$group,levels = c("Y","N"))
dds <- DESeqDataSetFromMatrix(countData = as.matrix(merge),colData = coldata,design = ~group)
dds <- DESeq(dds)
res <- results(dds,contrast = c("group","Y","N"))
resultsNames(dds) 
summary(res)  
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
res <- data.frame(gene_id=rownames(res),res)
rownames(res) <- NULL
rank <- merge(ID,res,by="gene_id") %>%
  dplyr::select(gene_name,log2FoldChange)
write.table(rank,file = "./result/FDrank.rnk",sep = "\t",row.names = F,col.names = F,quote = F)
#write.table(rank,file = "./result/FDrank2.rnk",sep = "\t",row.names = F,col.names = F,quote = F)

#control - results(dds,contrast = c("group","Y","N")) 
#> log2 fold change (MLE): group Y vs N  , second is control 

#grp file for GSEA 
grp <- read.table("./result/geneset.grp")
write.table(grp,"./result/grp.grp",quote = F,row.names = F,col.names = F)

#=====================================================
#clusterprofiler - GSEA 

gmt <- read.gmt("./result/gene_sets.gmt")
gmt <- unlist(gmt)
gmt <- as.data.frame(cbind(geneset="genelist",symbol=gmt))
rank <- read.table("./result/FDrank2.rnk",sep = "\t")
rank <- na.omit(rank[order(rank$V2,decreasing = T),])
gene_list <- rank$V2
names(gene_list) <- as.character(rank$V1)
gene <- names(gene_list)

egmt2 <-GSEA(gene_list,TERM2GENE=gmt,nPerm = 1000,minGSSize = 0,maxGSSize = 500,pvalueCutoff=0.5)

gseaplot2(egmt2,geneSetID=1,title = "",pvalue_table=T,color = "black")
results <- as.data.frame(egmt2@result)
save(results,file = "./result/gsea_result.RData")
