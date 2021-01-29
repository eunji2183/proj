count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
rownames(count) <- count[,1]
count <- count[,-1]
countData <- count[apply(count, 1, sum) > 1 , ]
names(countData)[13] <-'Y26_R' 
names(countData)[14] <- 'Y5_NR'
countData <- countData[,c(1:3,7:9,14,4:6,10:12,13)]
geneLength <- as.data.frame(fread("./hg38_gene_length.txt",fill = T,header = F))
names(geneLength)[1] <- 'gene_id'
names(geneLength)[2] <- 'gene_length'
