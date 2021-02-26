###  2021/02/26  ###
#protein oding gene transcriptome 

gtf <- rtracklayer::import("/proj2/ref/Homo_sapiens.GRCh38.102.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name,gene_biotype) %>% distinct()
coding_ID <- ID %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(gene_id,gene_name)
setwd("/home/eunji/proj/2D3D/")
count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
geneLength <- as.data.frame(fread("./hg38_gene_length.txt",fill = T,header = F))
names(geneLength)[1] <- 'gene_id'
names(geneLength)[2] <- 'gene_length'

#TPM normalization 
rownames(count) <- count[,1]
count <- count[,-1]
geneLength <- geneLength[order(geneLength$gene_id),]
rownames(geneLength) <- NULL

tpm <- function(counts,lengths){
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk)/1e6
  rpk/coef
} 

tpms <- apply(count,2,function(x)tpm(x,geneLength$gene_length))
tpms <- as.data.frame(tpms)
tpms <- data.frame(gene_id=rownames(tpms),tpms)
rownames(tpms) <- NULL
coding <- merge(coding_ID,tpms,by="gene_id")
coding <- coding[,-1]
coding <- coding[!duplicated(coding$gene_name),]
rownames(coding) <- coding[,1]
coding <- coding[,-1]
coding <- coding[,c(1:3,7:9,14,4:6,10:12,13)]
names(coding) <- rownames(Allcol)
merge <- coding

#FC for 2D 3D patient (R / NR) 
for (i in 1:nrow(merge)) {
  a <- merge[i,8:10]
  b <- merge[i,1:3]
  c <- mean(as.matrix(a))
  d <- mean(as.matrix(b))
  merge$FC2D[i] <- c/d 
}

for (i in 1:nrow(merge)) {
  a <- merge[i,11:13]
  b <- merge[i,4:6]
  c <- mean(as.matrix(a))
  d <- mean(as.matrix(b))
  merge$FC3D[i] <- c/d 
}

for (i in 1:nrow(merge)) {
  a <- merge[i,14]
  b <- merge[i,7]
  c <- mean(as.matrix(a))
  d <- mean(as.matrix(b))
  merge$FCT[i] <- c/d 
}

merge$log2FC2D = apply( merge, 1, function(t) {log2( mean(t[8:10])/ mean(t[1:3]))})
merge$log2FC3D = apply( merge, 1, function(t) {log2( mean(t[11:13])/ mean(t[4:6]))})
merge$log2FCT = apply( merge, 1, function(t) {log2( mean(t[14])/ mean(t[7]))})


#p value
tpmdf <- merge[,c(1:14)]
for(i in 1:nrow(tpmdf)){
  x=as.numeric(tpmdf[i, 8:10])
  y=as.numeric(tpmdf[i, 1:3])
  tpmdf$p.2D[i]=t.test(x, y,alternative = "two.sided")$p.value}

for(i in 1:nrow(tpmdf)){
  x=as.numeric(tpmdf[i, 11:13])
  y=as.numeric(tpmdf[i, 4:6])
  tpmdf$p.3D[i]=t.test(x, y,alternative = "two.sided")$p.value}

two <- merge %>%
  dplyr::select(`2D_NR1`,`2D_NR2`,`2D_NR3`,
                `2D_R1`,`2D_R2`,`2D_R3`,FC2D,log2FC2D)

two <- data.frame(gene_name=rownames(two),two)
rownames(two) <- NULL

two_coding_GSEA <- two %>%
  dplyr::filter(FC2D != Inf) %>%
  dplyr::filter(!is.nan(FC2D)) %>%
  dplyr::select(gene_name,FC2D) 

three <- merge %>%
  dplyr::select(`3D_NR1`,`3D_NR2`,`3D_NR3`,
                `3D_R1`,`3D_R2`,`3D_R3`,FC3D,log2FC3D)
three <- data.frame(gene_name=rownames(three),three)
rownames(three) <- NULL 

p <- merge %>%
  dplyr::select(Y5_NR,Y26_R,FCT,log2FCT)
p<- data.frame(gene_name=rownames(p),p)
rownames(p) <- NULL 

