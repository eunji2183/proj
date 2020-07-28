#DEG : microdust (con p5_con YJ_con upm p5_upm YJ_upm)

rm(list = ls())
working_dir="/home/eunji/R/project/200717_microdust/"
setwd(working_dir)
options(stringsAsFactors = F) #no character2factor

library(stringr) #str_sub
library(magrittr) # %>%
library(dplyr) # select,filter,distinct
library(DESeq2)
library(enrichR) #enrichr
library(clusterProfiler)
library(org.Hs.eg.db)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate",repos = rforge,dependencies = TRUE)
library(estimate)
library(msigdf)
#===================================================
#data : count-featurecount , ID-Homo_sapiens.GRCh38.96.gtf

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

resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1] <- 'gene_id'
resdata <- merge(ID,resdata,by="gene_id")
resdata <- resdata[order(resdata$padj,decreasing = F),]
resdata <- na.omit(resdata)

#gene_description

des <- read.table("./data/des.txt",sep = "\t",fill = T,quote = "",header = T)
names(des)[1] <- 'gene_name'
names(des)[3] <- 'gene_id'
des <- des %>%
  dplyr::select(gene_id,Gene.description)

resdata <- merge(x=des,y=resdata,by="gene_id",all.y=T)
write.csv(resdata,file="./result/con_vs_upm/DESeq2_con_vs_upm.csv")

#GSEA  cls gct  file 
group <- paste(group,collapse = " ")
group <- c(paste(c(3+3,2,1),collapse = " "), "# con upm",group)
write.table(file = "./result/group.cls",group,col.names = F,row.names = F,quote = F)

gct <- resdata %>%
  dplyr::select(gene_name,Gene.description,CON,p5_con,YJ_con,UPM,p5_upm,YJ_UPM)
names(gct)[1] <- 'NAME'
names(gct)[2] <- 'DESCRIPTION'
write.table(file = "./result/GSEA/FDGSEA.gct",gct,sep = "\t",col.names = T,row.names = F,quote = F)


#GSEA rank file (prerank)
rnk <- resdata %>%
  dplyr::select(gene_name,log2FoldChange)

write.table(rnk,file = "./result/GSEA/FDrank.rnk",sep = "\t",row.names = F,col.names = F,quote = F)


#=======================================================
#enrichr

GO_BP <- msigdb %>%
  dplyr::filter(msigdb$category_subcode == "bp") %>%
  dplyr::select(geneset, symbol) %>% 
  as.data.frame


dbs <- listEnrichrDbs()
head(dbs)
dbs <- c("KEGG_2019_Human")
GO_BP <- read.table("./data/GO_Biological_Process_2018.txt",sep = "\t",fill = T) 
library(knitr)
kable(head(dbs[c(1:3,5:6),]))
symbol <- data.frame(SYMBOL=resdata$gene_name) %>% dplyr::distinct()
symbol <- toupper(symbol$SYMBOL)
df <- bitr(symbol$SYMBOL,fromType = "SYMBOL",
           toType = c("ENSEMBL","ENTREZID"),
           OrgDb = org.Hs.eg.db)
enrichr <- enrichr(gene,dbs)

write.table(symbol,file="./result/GSEA/genelist.txt",row.names = F,col.names = F,quote = F)
