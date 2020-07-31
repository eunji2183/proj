data <- function(path){
  suppressMessages(library(stringr))
  suppressMessages(library(dplyr))
  suppressMessages(library(magrittr))
  filenames <- dir(path,recursive = T,pattern = ".txt",full.names = T)
  data <- lapply(filenames,function(x){
    read.table(x,header = T,sep = "\t",stringsAsFactors = F)})
  for(i in 1:length(data)){
    data[[i]] <- data[[i]][,c(1,7)] %>%
      reshape::rename(c(Geneid="gene_id"))}
  count <- Reduce(function(x,y)merge(x,y,by="gene_id"),data)
  rownames(count) <- count[,1]
  count <- count[,-1]
  name <- strsplit(filenames,split = "/")
  name <- sapply(name,function(x){str_subset(x,pattern = ".txt")})
  colnames(count) <- name
  colnames(count) <- str_sub(colnames(count),1,-5)
  return(count)  
}

DESeq2.1 <- function(count){
  suppressMessages(library(DESeq2))
  suppressMessages(library(dplyr))
  coldata <- data.frame(row.names = colnames(count),group)
  dds <- DESeqDataSetFromMatrix(countData = count,colData = coldata,design = ~group)
  dds <- DESeq(dds)
  res <- results(dds) 
  resultsNames(dds) 
  res <- as.data.frame(res)
  res <- res %>% dplyr::filter(!is.na(log2FoldChange))
  resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
  names(resdata)[1] <- 'gene_id'
  resdata <- merge(ID,resdata,by="gene_id")
  resdata <- resdata[order(resdata$pvalue),]
  resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 1,
                                  ifelse(resdata$log2FoldChange > 1 , 'Up','Down'),'not'))
  
  resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 2.5,resdata$gene_name,NA)
  return(resdata)
}

DESeq2.2 <- function(count){
  suppressMessages(library(TCC))
  suppressMessages(library(DESeq2))
  suppressMessages(library(dplyr))
  tcc <- new("TCC",count,group)
  tccNF <- calcNormFactors(tcc, method='tmm', test.method='deseq2', 
                           iteration = 3, FDR = 0.05, floorPDEG = 0.05)
  tccDE <- estimateDE(tccNF, test.method = 'deseq2', FDR = 0.5)
  tccRES <- getResult(tccDE, sort = T)
  tccRES <- merge(ID,tccRES,by="gene_id")
  tccRES <- dplyr::arrange(tccRES,rank)
  return(tccRES)
}



edgeR <- function(count){
  suppressMessages(library(TCC))
  suppressMessages(library(edgeR))
  suppressMessages(library(dplyr))
  tcc <- new("TCC",count,group)
  tccNF <- calcNormFactors(tcc, method='tmm', test.method='edger', 
                           iteration = 3, FDR = 0.05, floorPDEG = 0.05)
  tccDE <- estimateDE(tccNF, test.method = 'edger', FDR = 0.5)
  tccRES <- getResult(tccDE, sort = T)
  tccRES <- merge(ID,tccRES,by="gene_id")
  tccRES <- dplyr::arrange(tccRES,rank)
  return(tccRES)
}
