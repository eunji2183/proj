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
