#metabric CNA 
metacna <- read.table("./public/brca_metabric/data_CNA.txt",sep = "\t",header = T)
chr <- read.table("./public/chr.txt",sep = "\t")
names(chr) <- c("CHR","Hugo_Symbol")
metacna2 <- merge(chr,metacna,by="Hugo_Symbol")
metacna2 <- metacna2[,-3]
metacna2$CHR <- as.numeric(metacna2$CHR)
kmean <- metacna2[,-2]
kmean <- kmean[!duplicated(kmean$Hugo_Symbol),]
rownames(kmean) <- kmean[,1]
kmean <- kmean[,-1]

HC <- as.data.frame(t(kmean))
HCdist <- dist(HC, diag=TRUE)
hclust <- hclust(HCdist)

# Plot the result
plot(hclust)
plot(hclust, cex=0.2)
library(dendextend)
hclustdend <- as.dendrogram(hclust) # create dendrogram object
nleaves(hclustdend) #“leaves” (= # of genes we clustered) 
nnodes(hclustdend) #nodes (= # of leaves + number of internal joins)
clusters <- cutree(hclustdend, k=4)
table(clusters)
plot(color_branches(hclustdend, k=4),leaflab="none")
clusters.df <- data.frame(SAMPLE = names(clusters), cluster = clusters)
clusters.df["YAL022C",]



heat <- as.matrix(as.numeric(HC))

heatmap(heat)

genegroup_DEG <- function(cancer,N){
  for(i in 2:length(cancer)){
    a <- colgene %>%
      dplyr::filter(TYPE == cancer[i])
    rownames(a) <- a[,1]
    a <- a[,-c(1:3)]
    a <- a[,-c(1:length(colnames(a))/2)]
    if(N=1){
      b <- filter(a, a[,1] != "middle")
      c <- data.frame(row.names = rownames(b),group=b[,1])}
    else if (N=2){
      b <- filter(a, a[,2] != "middle")
      c <- data.frame(row.names = rownames(b),group=b[,2])}
    else{
      b <- filter(a, a[,1] != "middle",a[,2] != "middle")
      b$group <- ifelse(b[,1]==b[,2],b[,1],NA)
      c <- b%>%
        dplyr::filter(group != "NA")%>%
        dplyr::select(group)}
    c$group <- as.factor(c$group)
    c$group <- c[order(c$group),]
    d <- merge(c,expr,by="row.names")
    table(d$group)
    rownames(d) <- d[,1]
    d <- d[,-c(1:2)]
    d <- as.data.frame(t(d))
    e <- as.matrix(sapply(d,as.numeric))
    e[is.na(e)] <- 0
    row.names(e) <- rownames(d)
    f <- DESeqDataSetFromMatrix(countData = e,colData = c,design = ~group)
    f <- DESeq(f)
    h <- results(f) 
    resultsNames(f) 
    h <- as.data.frame(h)
    h <- h %>% dplyr::filter(!is.na(log2FoldChange))
    h <- h[order(h$pvalue),]
    h <- merge(as.data.frame(h),as.data.frame(counts(f,normalized=TRUE)),by="row.names",sort=FALSE)
    l[[i]] <- h
    names(l)[[i]] <- cancer[i]
    lcol[[i]] <- c
    names(lcol)[[i]] <- cancer[i]
  }
  m <- list(l,lcol)
  names(m)[[1]] <- "res"
  names(m)[[2]] <- "coldata"
  return(m)
}


