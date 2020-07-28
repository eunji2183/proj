#cachexia‐inducing factors in human cancers (12 types)
#TCGA RNA-seq data & GTEx - pancancer(12) - DEG - survival
#https://doi.org/10.1002/jcsm.12565 

rm(list = ls())
working_dir="/home/eunji/R/project/200629_cac/"
setwd(working_dir)
options(stringsAsFactors = F)
save(list = ls(),file="./data/cac.RData")
BiocManager::install("GenomicDataCommons")
BiocManager::install("CePa")
install.packages("enrichR")
install.packages("numDeriv")
library(dplyr)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(stringr)
library(GenomicDataCommons)
library(CePa)
library(readr)
library(tidyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichR)
library(msigdf)
library(numDeriv)
library(reshape2)
library(UCSCXenaTools)
library(rtracklayer) #gtf file 
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate",repos = rforge,dependencies = TRUE)
library(estimate)
library(ABSOLUTE)


#=======================================================
#data : TCGA & GTEx & secretome(HPA) database


#TCGA 12 cancer types RNA-seq HTSeq - Counts data , 3 times

#cancer  <- TCGAbiolinks:::getGDCprojects()$project_id
#cancer <- str_subset(cancer, "TCGA") #33 cancer types
#cancer <- sort(cancer)

rm(list = ls())
cancer <- c("TCGA-PAAD","TCGA-LAML","TCGA-READ","TCGA-COAD")
cancer <- c("TCGA-STAD","TCGA-LIHC","TCGA-PRAD","TCGA-BRCA")
cancer <- c("TCGA-LUAD","TCGA-LUSC","TCGA-HNSC","TCGA-ESCA")
for (i in length(cancer)) {
  cancer_select <- cancer[i]
  print(cancer)
  setwd("/home/eunji/R/project/200629_cac/data/TCGA/")
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
    expdat <- GDCprepare(query = query, save = TRUE,
                         save.filename = paste0(cancer_select,".rda"))
    setwd("/home/eunji/R/project/200629_cac/data/TCGA/")
    #data-prepeocessing ( sprearman > 0.6)
    #expdat <- TCGAanalyze_Preprocessing(object = expdat,
                                        #cor.cut = 0.6,
                                        #filename = paste0(cancer_select,".png"),
                                        #datatype = "HTSeq - Counts")
    count_matrix=assay(expdat)
    write.csv(count_matrix,
              file = paste( cancer_select,"Counts.csv",
                            sep = "-"))}}


rm(list = ls())

#TCGA clinical data : GDC

cancer <- c("PAAD","LAML","READ","COAD")
cancer <- c("STAD","LIHC","PRAD","BRCA")
cancer <- c("LUAD","LUSC","HNSC","ESCA")

for (i in 1:length(cancer)){
  cancer_select <- cancer[i]
  print(cancer_select)
  setwd("/home/eunji/R/project/200629_cac/data/TCGA_clinical_GDC")
  query <- GDCquery(project = cancer_select,
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query,clinical.info = "patient")
  write.csv(clinical,file = paste0(cancer_select,'.csv'))
}

rm(clinical,query,cancer_select,i,cancer)


#TCGA clinical data : UCSC xena 

cancer <- c("PAAD","LAML","READ","COAD")
cancer <- c("STAD","LIHC","PRAD","BRCA")
cancer <- c("LUAD","LUSC","HNSC","ESCA")

for (i in 1:length(cancer)){
  cancer_select <- cancer[i]
  print(cancer_select)
  setwd("/home/eunji/R/project/200629_cac/data/TCGA_clinical_UCSC/")
  getTCGAdata(project = cancer_select)
  pheno=getTCGAdata(project = cancer_select,download = TRUE,forceDownload = TRUE)
  clinical=XenaPrepare(pheno) 
  write.csv(clinical,file = paste0(cancer_select,'.csv'))
}

rm(clinical,pheno)



#GTEx 

GTEx <- read.gct("./data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
gtex <- as.data.frame(GTEx)
rownames(gtex) <- make.names(str_sub(rownames(gtex),1,15),unique = T)
set.seed(1)
gtex <- gtex[,sample(ncol(gtex),100,replace = F)]
gtex <- rbind(gtex,group=rep("N",100))
gtex <- data.frame(gene_id=rownames(gtex),gtex)
rownames(gtex) <- NULL
save(gtex,file = "./data/gtex.RData")


#secretome (HPA cite secretome genes)

signalp <- read.csv("./data/secretome/protein_class_SignalP.tsv",sep = "\t",header = T,stringsAsFactors = F)
spoctopus <- read.csv("./data/secretome/protein_class_SPOCTOPUS.tsv",sep = "\t",header = T,stringsAsFactors = F)
phobius <- read.csv("./data/secretome/protein_class_Phobius.tsv",sep = "\t",header = T,stringsAsFactors = F)
secretome <- rbind(signalp,spoctopus,phobius) %>% unique()
names(secretome)[1] <- 'gene_name' 
names(secretome)[3] <- 'gene_id'
secretome <- secretome[,c(1,3:11,49)]
save(secretome,file ="./data/secretome.RData")

#paper- secretome genes 
paper_secret <- readxl::read_excel("./data/secretome/jcsm12565-sup-0011-sm_data-sd1-sd2.xlsx",
                                   sheet = "Data S1",
                                   range = "A3:D2936", #or skip=2
                                   col_names = T,
                                   na="NA")

names(paper_secret)[1] <- 'gene_name'
names(paper_secret)[3] <- 'gene_id'
paper_secret <- paper_secret[,c(1,3,4)]

save(paper_secret,file = "./data/paper_secret.RData")


rm(query,expdat,count_matrix,GTEx,phobius,signalp,spoctopus)


#======================================================
# merge GTEx & TCGA 
working_dir="/home/eunji/R/project/200629_cac/"
setwd(working_dir)
options(stringsAsFactors = F)
TGC_dir <- c("./data/TCGA/")
TGC_filelist <- str_subset(list.files(TGC_dir),"csv")

for(i in 1:length(TGC_filelist)) {
  setwd("/home/eunji/R/project/200629_cac/data/TCGA/")
  a <- read.csv(TGC_filelist[i],header = T,sep = ",",stringsAsFactors = F,row.names=1)
  rownames(a) <- make.names(str_sub(rownames(a),1,15),unique = T)
  colnames(a) <- str_sub(colnames(a),1,15)
  colnames(a) <- gsub('.','-',colnames(a),fixed = T)
  a <- as.data.frame(t(a))
  a$group <- as.numeric(as.character(substring(rownames(a),14,15)))
  a$group <- ifelse(a$group < 10 , 'T','N')
  print(table(a$group))
  a <- as.data.frame(t(a))
  a <- data.frame(gene_id=rownames(a),a)
  a <- merge(a,gtex,by="gene_id")
  rownames(a) <- a[,1]
  a <- a[,-1]
  a <- as.data.frame(t(a))
  tumor <- a %>% 
    dplyr::filter(group == "T")
  normal <- a %>%
    dplyr::filter(group == "N")
  merge <- rbind(tumor,normal)
  coldata <- merge %>%
    dplyr::select(group)
  coldata <- data.frame(row.names = rownames(coldata),coldata)
  tumor <- as.data.frame(t(tumor))
  normal <- as.data.frame(t(normal))
  a <- cbind(tumor,normal)
  a <- a[c(1:(length(rownames(a)) -1)),]
  type <- str_sub(TGC_filelist[i],6,9)
  setwd("/home/eunji/R/project/200629_cac/result/coldata/")
  write.csv(coldata,file = paste0(type,'.csv'))
  setwd("/home/eunji/R/project/200629_cac/result/merge/")
  write.csv(a,file = paste0(type,'.csv'))}

rm(a,coldata,tumor,normal,merge)

#========================================================
# secretome genes DEG : DESeq2 

working_dir="/home/eunji/R/project/200629_cac/"
setwd(working_dir)
options(stringsAsFactors = F)
merge_dir <- c("./result/merge/")
merge_filelist <- str_subset(list.files(merge_dir),"csv")

for (i in 1:length(merge_filelist)) {
  setwd("/home/eunji/R/project/200629_cac/result/merge/")
  a <- read.csv(merge_filelist[i],header = T,sep = ",",stringsAsFactors = F)
  names(a)[1] <- 'gene_id'
  a <- merge(paper_secret,a,by="gene_id")
  count <- data.frame(row.names = a$gene_id,a[,-c(1:3)])
  setwd("/home/eunji/R/project/200629_cac/result/coldata/")
  coldata <- read.csv(merge_filelist[i],header = T,sep = ",",stringsAsFactors = F,row.names = 1)
  coldata$group <- factor(coldata$group,levels = c("T","N"))
  dds <- DESeqDataSetFromMatrix(countData = count,colData = coldata,design = ~group)
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("group","T","N"))
  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
  names(resdata)[1] <- 'gene_id'
  a <- merge(paper_secret,resdata,by="gene_id")
  a <- a[order(a$padj),]
  type <- str_sub(merge_filelist[i],1,4)
  setwd("/home/eunji/R/project/200629_cac/result/DEG/")
  write.csv(a,file = paste0(type,'.csv'))}


rm(a,coldata,count,dds,res,resdata)

#======================================================
#pheatmap - padj < 0.01 , FC > 2 
#secreted protein-coding genes(LogFC)-12 cancer 

setwd("/home/eunji/R/project/200629_cac/")
path <- "./result/DEG"
filenames <- dir(path)
filepath <- sapply(filenames,function(x){
  paste(path,x,sep = '/')})
data <- lapply(filepath,function(x){
  read.csv(x,header = T,sep = ",",stringsAsFactors = F,row.names = 1)})
names(data) <- str_sub(names(data),end = 4)
save(data,file = "./data/data.RData")

for (i in 1:length(data)) {
  data[[i]] <- data[[i]] %>%
    dplyr::select(gene_name,log2FoldChange,padj) %>%
    dplyr::filter(padj < 0.01) %>%
    dplyr::select(gene_name,log2FoldChange) %>%
    reshape::rename(c(log2FoldChange = names(data)[[i]]))}


FC <- Reduce(function(x,y)merge(x,y,by="gene_name"),data)
rownames(FC) <- FC[,1]
FC <- FC[,-1]

p_FC <- pheatmap(FC,scale = "row", clustering_distance_row = "correlation",
                 show_colnames = T,show_rownames = F,
                 cluster_cols = T,cluster_rows = T,
                 color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                 border_color = 'white',
                 display_numbers = F,main = "",key=T)
                 #legend = T,legend_breaks = c(min(FC),0,5,max(FC)),
                 #legend_labels = c(min(FC),"0","5","log2FC"))

save(data,file = "./data/FC_list.RData")
save(FC,file = "./data/FC.RData")


rm(p_FC)

#=======================================================
#PCA 

pca <- prcomp(t(FC),scale=T)
pca.data <- data.frame(type=rownames(pca$x),X=pca$x[,1],Y=pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
ggplot(pca.data,aes(X,Y,color=type,label=type)) +
  geom_text(size=3,vjust=0,nudge_y = 0.1) +
  geom_point() +
  xlab(paste0("PC1: ",pca.var.per[1],"% variance")) +
  ylab(paste0("PC2: ",pca.var.per[2],"% variance")) +
  coord_fixed()

rm(pca,pca.data)
#=======================================================
#barplot

load("./data/data.RData")

for (i in 1:length(data)){
  data[[i]] <- data[[i]] %>%
    dplyr::select(gene_name,log2FoldChange,padj) %>%
    dplyr::filter(padj < 0.01) %>%
    dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
    dplyr::select(gene_name,log2FoldChange) %>%
    dplyr::mutate(regulate=ifelse(log2FoldChange > 1,'up','down'))
  data[[i]] <- data[[i]] %>%
    dplyr::count(regulate) %>%
    reshape::rename(c(n = names(data)[[i]]))}

up_down <- Reduce(function(x,y)merge(x,y,by="regulate"),data) 
rownames(up_down) <- up_down[,1]
up_down <- up_down[,-1]
up_down <- as.matrix(up_down)

barplot(up_down, axes = T, main = "",
        horiz = T,beside = F, #horiz (y axes)
        legend.text = T,
        args.legend=list(x="topright",cex=0.5,text.width=70))



#=======================================================
#KEGG & GO pathway analysis

en_sym<-select(org.Hs.eg.db,keys = keys(org.Hs.eg.db),columns = c('ENSEMBL','SYMBOL','ENTREZID'))
names(en_sym)[3] <- 'gene_name'
en_sym <- en_sym[,c(1,3)] %>% unique()
save(en_sym,file = "./data/en_sym.RData")

load("./data/data.RData")

for (i in 1:length(data)){
  data[[i]] <- data[[i]] %>%
    dplyr::select(gene_name,log2FoldChange,padj) %>%
    dplyr::filter(padj < 0.01) %>%
    dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
    dplyr::select(gene_name,log2FoldChange) %>%
    dplyr::mutate(regulate=ifelse(log2FoldChange > 1,'up','down')) %>%
    dplyr::filter(regulate=="up")
  data[[i]] <- merge(en_sym,data[[i]],by="gene_name")[,"ENTREZID"]
  }

up <- data
save(up,file = "./data/up.RData")

load("./data/data.RData") 

for (i in 1:length(data)){
  data[[i]] <- data[[i]] %>%
    dplyr::select(gene_name,log2FoldChange,padj) %>%
    dplyr::filter(padj < 0.01) %>%
    dplyr::filter(log2FoldChange > 1 | log2FoldChange < -1) %>%
    dplyr::select(gene_name,log2FoldChange) %>%
    dplyr::mutate(regulate=ifelse(log2FoldChange > 1,'up','down')) %>%
    dplyr::filter(regulate=="down")
  data[[i]] <- merge(en_sym,data[[i]],by="gene_name")[,"ENTREZID"]}

down <- data
save(down,file = "./data/down.RData")
  

kegg_up <- compareCluster(up,fun = "enrichKEGG",
                          pvalueCutoff=0.05)
kegg_down <- compareCluster(down,fun = "enrichKEGG",
                            pvalueCutoff=0.05)                          
dotplot(kegg_up,showCategory=5,includeAll=F)
dotplot(kegg_down,showCategory=5,includeAll=F)


GO_up <- compareCluster(up,fun = "enrichGO",
                        OrgDb='org.Hs.eg.db',
                        ont='BP',pvalueCutoff=0.05)
GO_down <- compareCluster(down,fun = "enrichGO",
                        OrgDb='org.Hs.eg.db',
                        ont='BP',pvalueCutoff=0.05)
dotplot(GO_up,showCategory=5,includeAll=F)
dotplot(GO_down,showCategory=5,includeAll=F)
ggplot(GO_up,aes(Cluster,Description),showCategory=5,includeAll=F) +
  geom_point(aes(color=p.adjust,size=GeneRatio))

  
#======================================================
#cachexia-inducing factors(CIF) - expression profiles
#dotplot 

load("./data/data.RData")
gene_list = c("CXCL8","IL1B","HGF","TNFSF10","LIF",
              "TGFA","TNFSF11","PDGFB","IL6","CSF3",
              "CCL2","CSF1","IL15","CSF2","TNF",
              "VEGFA","LEP","FGF2","CXCL12","MMP13",
              "IL15","IL10","CD40LG","IFNA1","IL4")

gene <- data.frame(gene_name=gene_list)
save(gene,file = "./data/CIF_25.RData")

gene_count=data.frame()


for(i in 1:length(data)){
  data[[i]] <- merge(gene,data[[i]],by="gene_name") %>%
    dplyr::select(gene_name,log2FoldChange,padj) %>%
    dplyr::mutate(type=names(data)[[i]])
  a <- data[[i]]
  gene_count <- rbind(gene_count,a)}


gene_count$log10 <- -log10(gene_count$padj)
names(gene_count)[5] <- "-log10(padj)" 
count <- gene_count %>%
  dplyr::select(gene_name,type,log2FoldChange,"-log10(padj)")
count$`-log10(padj)`[is.infinite(count$`-log10(padj)`)] <- 225
count <- count %>%
  dplyr::filter(`-log10(padj)` > 1.30103) #padj > 0.05


ggplot(count,aes(type,gene_name),showCategory=25,includeAll=F) +
  geom_point(aes(color=log2FoldChange,
                 size=`-log10(padj)`)) +
  scale_size_area(max_size = 8) +
  scale_color_gradient2(low = "blue",midpoint = 0,mid = "white",high = "red") +
  #scale_color_continuous(low='blue','white',high='red') +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())
        #axis.line = element_line(colour = "black"))

#======================================================
#enrichR

dbs <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018")
msigdb <- as.data.frame(msigdf.human)
GO_BP <- msigdb %>%
  dplyr::filter(msigdb$category_subcode == "bp") %>%
  dplyr::select(geneset, symbol) %>% 
  as.data.frame


load("./data/data.RData")


for(i in 1:length(data)){
  genelist <- data[[i]]$log2FoldChange
  names(genelist) <- data[[i]]$gene_name
  genelist_tr <- bitr(names(genelist),
                      fromType = "SYMBOL",
                      toType = c("ENSEMBL","ENTREZID"),
                      OrgDb = org.Hs.eg.db)
  new_list <- data.frame(SYMBOL=names(genelist),logFC = as.numeric(genelist))
  new_list <- merge(new_list,genelist_tr,by="SYMBOL") %>%
    dplyr::select(SYMBOL,ENTREZID,logFC) %>% unique()
  genelist <- new_list$logFC
  names(genelist) <- new_list$SYMBOL
  data[[i]] <- sort(genelist,decreasing = T)}
  
for (i in 1:length(data)){
  data[[i]] <- data[[i]] %>%
    dplyr::filter(log2FoldChange > 0 & padj < 0.05)
  genelist <- data[[i]]$log2FoldChange
  names(genelist) <- data[[i]]$gene_name
  data[[i]] <- sort(genelist,decreasing = T)}

for(i in 1:length(data)){
  data[[i]] <- enrichr(data[[i]],GO_BP)
}



#======================================================
#tumor purity 

#ESTIMATE purity 

setwd("/home/eunji/R/project/200629_cac/")
path <- "./data/TCGA"
filenames <- dir(path)
filelist <- str_subset(filenames,"csv")
filepath <- sapply(filelist,function(x){
  paste(path,x,sep = '/')})
TCGAcount <- lapply(filepath,function(x){
  read.csv(x,header = T,sep = ",",stringsAsFactors = F,row.names = 1)})
names(TCGAcount) <- str_sub(names(TCGAcount),6,9)

save(TCGAcount,file = "./data/TCGAcount.RData")
load("./data/gtf_df.RData")

for(i in 1:length(TCGAcount)) {
  a <- TCGAcount[[i]]
  rownames(a) <- make.names(str_sub(rownames(a),1,15),unique = T)
  colnames(a) <- str_sub(colnames(a),1,15)
  colnames(a) <- gsub('.','-',colnames(a),fixed = T)
  a <- as.data.frame(t(a))
  a$group <- as.numeric(as.character(substring(rownames(a),14,15)))
  a$group <- ifelse(a$group < 10 , 'T','N')
  print(table(a$group))
  a <- a %>%
    dplyr::filter(group == 'T')
  a <- as.data.frame(t(a))
  a <- data.frame(gene_id=rownames(a),a)
  a <- merge(gtf,a,by="gene_id")
  a$gene_id <- NULL
  names(a)[1] <- 'GeneSymbol'
  a$GeneSymbol <- make.names(a$GeneSymbol,unique = T)
  setwd("/home/eunji/R/project/200629_cac/result/estimate/input_f/")
  write.table(a,file = paste0(names(TCGAcount)[[i]],'.txt'),sep = "\t",col.names = T,row.names = F,quote = F)
}


working_dir="/home/eunji/R/project/200629_cac/"
setwd(working_dir)
options(stringsAsFactors = F)
input_dir <- c("./result/estimate/input_f")
input_filelist <- str_subset(list.files(input_dir),"txt")

for (i in 1:length(input_filelist)) {
  type <- str_sub(input_filelist[i],1,4)
  setwd("/home/eunji/R/project/200629_cac/result/estimate/input_f")
  filterCommonGenes(input.f = input_filelist[[i]],
                    output.f = paste0(paste0(type,'.gct')),
                    id='GeneSymbol')}


working_dir="/home/eunji/R/project/200629_cac/"
setwd(working_dir)
options(stringsAsFactors = F)
gct_dir <- c("./result/estimate/input_f")
gct_filelist <- str_subset(list.files(gct_dir),"gct")

for (i in 1:length(gct_filelist)){
  type <- str_sub(gct_filelist[i],1,4)
  setwd("/home/eunji/R/project/200629_cac/result/estimate/input_f")
  estimateScore(input.ds = gct_filelist[i],
                output.ds = paste0(paste0(type,'_score.gct')),
                platform = "affymetrix")}
 

working_dir="/home/eunji/R/project/200629_cac/"
setwd(working_dir)
options(stringsAsFactors = F)
score_dir <- c("./result/estimate/input_f")
score_filelist <- str_subset(list.files(score_dir),"_score.gct")  

for(i in 1:length(score_filelist)){
  setwd("/home/eunji/R/project/200629_cac/result/estimate/input_f/")
  type <- str_sub(score_filelist[i],1,4)
  scores <- read.table(score_filelist[[i]],skip = 2,header = T)
  rownames(scores) <- scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  scores <- as.data.frame(scores)
  scores <- data.frame(sampleID=rownames(scores),scores)
  setwd("/home/eunji/R/project/200629_cac/result/estimate/output_score/")
  write.table(scores,file =paste0(type,'.csv'),sep = ",",col.names = T,row.names = F,quote = F) 
}



#ABSOLUTE purity 

absolute <- read.table("./data/absolute_purity.txt",sep = "\t",fill = T,header = T,quote = "") %>%
  dplyr::select(array,purity) %>%
  reshape::rename(c(array='sampleID',purity='ABSOLUTE_Purity'))

load("./data/data.RData")
load("~/R/project/200629_cac/data/CIF_25.RData")


for(i in 1:length(data)) {
  data[[i]] <- merge(gene,data[[i]],by="gene_name")
  a <- data[[i]][,-c(2:9)] %>% unique()
  rownames(a) <- a[,1]
  a <- a[,-1]
  a <- as.data.frame(t(a))
  rownames(a) <- gsub('.','-',rownames(a),fixed = T)
  a <- log2(a+1)
  a <- data.frame(sampleID=rownames(a),a)
  data[[i]] <- merge(absolute,a,by="sampleID") %>%
    dplyr::filter(!is.na(ABSOLUTE_Purity)) %>%
    dplyr::mutate(Purity_class=ifelse(ABSOLUTE_Purity > median(ABSOLUTE_Purity),'high','low'))}


#DNA hypermethylation mode purity 

setwd("/home/eunji/R/project/200629_cac/")
path <- "./data/methylation_purity"
filenames <- dir(path)
filepath <- sapply(filenames,function(x){
  paste(path,x,sep = '/')})
methylation_purity <- lapply(filepath,function(x){
  read.csv(x,header = T,sep = "\t",stringsAsFactors = F)})
names(methylation_purity) <- str_sub(names(methylation_purity),end = 4)

data$READ <- NULL  

for(i in 1:length(methylation_purity)){
  methylation_purity[[i]]$SampleName <- str_sub(methylation_purity[[i]]$SampleName,1,15) 
  names(methylation_purity[[i]])[1] <- 'sampleID'
  methylation_purity[[i]] <- methylation_purity[[i]] %>%
    dplyr::select(sampleID,Purity_InfiniumPurify)
}

for (i in 1:11){
  data[[i]] <- merge(methylation_purity[[i]],data[[i]],by="sampleID")
}

for (i in 1:length(data)){
  a <- data[[i]] %>%
    dplyr::select(sampleID,
                  ABSOLUTE_Purity,Purity_class)
  b <- data[[i]] %>%
    dplyr::select(-c(Purity_InfiniumPurify,
                     ABSOLUTE_Purity,Purity_class))
  a <- unique(a)
  b <- unique(b)
  a <- as.factor(a$Purity_class)
  rownames(a) <- a[,1]
  a <- a[,-1]
  rownames(b) <- b[,1]
  b <- b[,-1]
}

pheatmap(BRCA2,annotation_col = BRCA1,scale = "row", clustering_distance_row = "correlation",
                 show_colnames = T,show_rownames = F,
                 cluster_cols = T,cluster_rows = T,
                 color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                 border_color = 'white',
                 display_numbers = F,main = "",key=T)

#====================================================
#gene expression in normal vs tumor 

setwd("/home/eunji/R/project/200629_cac/data/")
#gtf <- rtracklayer::import("/home/eunji/project/0_sh/ref/rna/Homo_sapiens.GRCh38.96.gtf")
#gtf <- as.data.frame(gtf)
#save(gtf,file = "./gtf.RData")

load("./gtf.RData")
gtf <- subset(gtf,select=c("gene_id","gene_name")) %>% distinct()
save(gtf,file = "./gtf_df.RData") 

#gtf <- subset(gtf,select=c("gene_id","gene_name","gene_biotype"))
#gtf_protein_coding <- subset(gtf,gtf$gene_biotype == 'protein_coding')
#gtf_protein_coding <- gtf_protein_coding %>% distinct()

setwd("/home/eunji/R/project/200629_cac/data/TCGA/")
rm(list = ls())
load("/home/eunji/R/project/200629_cac/data/gtf_df.RData")
filelist <- dir(pattern = "csv$")

gene <- 'IL6'

for(i in 1:length(filelist)){
  setwd("/home/eunji/R/project/200629_cac/data/TCGA/")
  a <- read.csv(filelist[i],sep = ",",header = T)
  names(a)[1] <- 'gene_id'
  a <- merge(gtf,a,by='gene_id')
  a$gene_id <- NULL
  a <- subset(a,a$gene_name == 'IL6')
  a <- a %>%
    remove_rownames() %>%
    column_to_rownames(var = 'gene_name')
  a <- as.data.frame(t(a))
  print(dim(a))
  a$group <- as.numeric(as.character(substring(rownames(a),14,15)))
  a$group <- ifelse(a$group < 10,'T','N')
  print(table(a$group))
  a$type <- filelist[i]
  type1 <- strsplit(a$type[1],"-")
  type1 <- type1[[1]][2]
  a$type <- type1
  setwd("/home/eunji/R/project/200629_cac/result/expr/")
  write.csv(a,file = paste0('IL6_',type1,'.csv'))}








#=====================================================
#survival 

setwd("/home/eunji/R/project/200629_cac/data/TCGA_clinical_GDC/")

rm(list = ls())

filelist <- dir(pattern = "csv$")

for(i in 1:length(filelist)){
  if(i==1){
    a <- read.csv(filelist[i],sep = ",",header = T,row.names = 1)
    a <- a %>%
      dplyr::select('bcr_patient_barcode',
                    'vital_status',
                    'days_to_last_followup',
                    'days_to_death') %>%
      dplyr::rename(Barcode = 'bcr_patient_barcode',
                    OS = 'vital_status',
                    Time1 = 'days_to_death',
                    Time2 = 'days_to_last_followup')
    a[,'Time1'][is.na(a[,'Time1'])] <- 0
    a[,'Time2'][is.na(a[,'Time2'])] <- 0
    a$OS.Time <- a$Time2 + a$Time1
    dt <- a
    group <- as.character(filelist[i])
    dt$group <- group
    dt$Time2 <- NULL
    dt$Time1 <- NULL
  }
  else{
    a <- read.csv(filelist[i],sep = ",",header = T,row.names = 1)
    a <- a %>%
      dplyr::select('bcr_patient_barcode',
                    'vital_status',
                    'days_to_last_followup',
                    'days_to_death') %>%
      dplyr::rename(Barcode = 'bcr_patient_barcode',
                    OS = 'vital_status',
                    Time1 = 'days_to_death',
                    Time2 = 'days_to_last_followup')
    a[,'Time1'][is.na(a[,'Time1'])] <- 0
    a[,'Time2'][is.na(a[,'Time2'])] <- 0
    a$OS.Time <- a$Time2 + a$Time1
    dt1 <- a
    group <- as.character(filelist[i])
    dt1$group <- group
    dt1$Time2 <- NULL
    dt1$Time1 <- NULL
    dt <- rbind(dt,dt1)}}
dt$group <- substr(dt$group,start = 6,stop = nchar(dt$group)-4)
write.csv(dt,file = 'TCGA_survival.csv')
  
