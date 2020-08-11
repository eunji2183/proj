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

DESeq2.1 <- function(count,group,ID){
  suppressMessages(library(DESeq2))
  suppressMessages(library(dplyr))
  suppressMessages(library(magrittr))
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
                                  ifelse(resdata$log2FoldChange > 1 , 'UP','DOWN'),'NOT'))
  
  resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 2.5,resdata$gene_name,NA)
  return(resdata)
}

DESeq2.2 <- function(count,group,ID){
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

edgeR.1 <- function(count,group,ID){
  gene_id <- rownames(count)
  dgelist <- DGEList(counts = count,genes = gene_id,group = group)
  dgelist_norm <- calcNormFactors(dgelist,method='TMM')
  design <- model.matrix(~group)
  dge <- estimateDisp(dgelist_norm,design,robust = T)
  cpm=cpm(dgelist)
  lcpm=cpm(dgelist,log = T)
  et <- exactTest(dge)
  tTags <- topTags(et,n=nrow(dgelist$counts))
  tTags <- as.data.frame(tTags)
  names(tTags)[1] <- 'gene_id'
  rownames(tTags) <- NULL
  tTag <- merge(ID,tTags,by="gene_id")
  tTag <- tTag %>%
    dplyr::filter(logFC != 0)
  tTag <- tTag[order(tTag$PValue),]
  return(tTag)
}

#count=count,group=group(coldata),ID=ID
edgeR.2 <- function(count,group,ID){
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

#TPM FPKM normalization 
countToTpm <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))}

countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )}

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}

countToEffCounts <- function(counts, len, effLen){
  counts * (len / effLen)}                  
                  
#df=resdata(dataframe) , gene=gene(character), raw=count
heatmap <- function(df,gene,raw){
  suppressMessages(library(magrittr))
  suppressMessages(library(dplyr))
  suppressMessages(library(pheatmap))
  heat <- df %>%
    dplyr::select(gene_name,colnames(raw))
  heat = heat[-which(duplicated(heat$gene_name)),]
  rownames(heat) <- heat[,1]
  heat <- heat[,-1]
  heat <- heat %>%
    dplyr::filter(rownames(heat) %in% gene)
  p <- pheatmap(heat,scale = "row", clustering_distance_row = "correlation",
                show_colnames = T,show_rownames = T,
                cluster_cols = T,cluster_rows = T,
                color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                border_color = 'white',
                display_numbers = F,main = "",key=T,fontsize_row = 10)
  return(p)
}

#res=resdata , raw=count ,gene=gene(character),ID=ID

datatable <- function(res,raw,gene,ID){
  suppressMessages(library(magrittr))
  suppressMessages(library(dplyr))
  heat <- res %>%
    dplyr::select(gene_name,colnames(raw),pvalue)
  heat = heat[-which(duplicated(heat$gene_name)),]
  heat <- heat %>%
    dplyr::filter(gene_name %in% gene)
  for(i in 1:length(colnames(raw))){
    names(raw)[i] <- paste0('Raw_',names(raw)[i])
    }
  rawdata <- data.frame(gene_id=rownames(raw),raw)
  rownames(rawdata) <- NULL
  rawdata <- merge(ID,rawdata,by="gene_id")
  rawdata$gene_id <- NULL
  data <- merge(rawdata,heat,by="gene_name")
  data <- data[order(data$pvalue),]
  return(data)
}                  
                  
#GO&KEGG : res=resdata,p=0.05(pvalue cutoff),FC=1(logFC cutoff)
GO_KEGG <- function(res,p,FC){
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(magrittr))
  suppressMessages(library(dplyr))
  DE <- res %>%
    dplyr::filter(pvalue < p & abs(log2FoldChange) > FC)
  gene <- DE$gene_name
  gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  de <- gene$ENTREZID
  go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
  EGG <- enrichKEGG(gene = de,organism = 'hsa',pvalueCutoff = 0.05)
  return(list(go,EGG))
}
