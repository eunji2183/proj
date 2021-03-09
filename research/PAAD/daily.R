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
three_coding_GSEA <- three%>%
  dplyr::filter(FC3D != Inf) %>%
  dplyr::filter(!is.nan(FC3D)) %>%
  dplyr::select(gene_name,FC3D)
              
p <- merge %>%
  dplyr::select(Y5_NR,Y26_R,FCT,log2FCT)
p<- data.frame(gene_name=rownames(p),p)
rownames(p) <- NULL 
p_coding_GSEA <- p %>%
  dplyr::filter(FCT != Inf) %>%
  dplyr::filter(!is.nan(FCT)) %>%
  dplyr::select(gene_name,FCT)

              
### 2021/03/02 ###
#protein-coding > PCA 
setwd("/home/eunji/proj/2D3D/")
count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
#GTF
gtf <- rtracklayer::import("/proj2/ref/Homo_sapiens.GRCh38.102.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name,gene_biotype) %>% distinct()
save(ID,file="./ID.RData")

coding_ID <- ID %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(gene_id,gene_name)

save(coding_ID,file = "./protein_codingID.RData")

coding_merge <- merge(coding_ID,count,by="gene_id")

coding_merge <- coding_merge[!duplicated(coding_merge$gene_name),]
coding_merge<- coding_merge[,-1]
rownames(coding_merge) <- coding_merge[,1]
coding_merge<- coding_merge[,-1]

coding_merge <- coding_merge[,c(1:3,7:9,14,4:6,10:12,13)]
names(coding_merge) <- rownames(Allcol)

#PCA
codingPCA <- as.matrix(sapply(coding_merge,as.numeric))
codingPCA[is.na(codingPCA)] <- 0
row.names(codingPCA) <- rownames(coding_merge)
dds <- DESeqDataSetFromMatrix(countData = codingPCA,colData = Allcol,design = ~RESPONSE)
dds <- DESeq(dds)
vsd <- vst(dds, blind=TRUE)
vsd <- vst(dds)
pcaData <- plotPCA(vsd,intgroup="RESPONSE",returnData=TRUE)
percentVar <- round(100*attr(pcaData,"percentVar"))
ggplot(pcaData,aes(PC1,PC2,color=group,shape=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  geom_label_repel(aes(label=name), 
                   fontface="bold", color="grey50", box.padding=unit(0.35, "lines"), 
                   point.padding=unit(0.2, "lines"), segment.colour = "grey50",size=02)
#GSEA result pheatmap 
fcmerge <- read.table("/home/eunji/gsea_home/output/feb24/keggcoding/merge2.tsv",
                      sep = "\t",header = T,row.names = 1,fill = T)

names(fcmerge) <- c('2D','2D_p','3D','3D_p','Patient','Patient_p')
p <- fcmerge %>%
  dplyr::select(`2D`,`3D`,Patient)

p1 <- p %>%
  dplyr::filter(`3D` < 0 & Patient < 0) %>%
  dplyr::filter(`2D` > 0)
p2 <- p %>%
  dplyr::filter(`3D` > 0 & Patient > 0) %>%
  dplyr::filter(`2D` < 0)  

p <- rbind(p1,p2)

paletteLength <- 30
myColor <- colorRampPalette(c('navy','white','red'))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(p), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(p)/paletteLength, max(p), length.out=floor(paletteLength/2)))

pheatmap(p, cluster_cols = T,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=12,fontsize_row=11, cellwidth = 30, cellheight = 15,
         color=myColor, breaks=myBreaks)

fcmerge_kegg <- read.table("/home/eunji/gsea_home/output/feb24/hallmark_coding_all/merge2.tsv",
                           sep = "\t",header = T,row.names = 1,fill = T)

names(fcmerge_kegg) <- c('2D','2D_p','3D','3D_p','Patient','Patient_p')
k <- fcmerge_kegg %>%
  dplyr::select(`2D`,`3D`,Patient)

k1 <- k %>%
  dplyr::filter(`3D` < 0 & Patient < 0) %>%
  dplyr::filter(`2D` > 0)
k2 <- k %>%
  dplyr::filter(`3D` > 0 & Patient > 0) %>%
  dplyr::filter(`2D` < 0)  

k <- rbind(k1,k2)

all <- rbind(p,k)

all <- all %>%
  dplyr::filter(`2D` != "---")

all[,1] <- as.numeric(all[,1])
all[,2] <- as.numeric(all[,2])
paletteLength <- 30
myColor <- colorRampPalette(c('navy','white','red'))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(all), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(all)/paletteLength, max(all), length.out=floor(paletteLength/2)))

pheatmap(all, cluster_cols = T,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=12,fontsize_row=8, cellwidth = 20, cellheight = 8,
         color=myColor, breaks=myBreaks)
       
### 2021/03/04 ###
#PCA sample 분리하는데 공헌많이한 gene 

setwd("/home/eunji/proj/2D3D/")
count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
#GTF
gtf <- rtracklayer::import("/proj2/ref/ensembl102/Homo_sapiens.GRCh38.102.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name,gene_biotype) %>% distinct()
save(ID,file="./ID.RData")

count <- merge(ID,count,by="gene_id")
count <- count[,-c(1,3)]
count <- count[!duplicated(count$gene_name),]
rownames(count) <- count[,1]
count <- count[,-1]

install.packages("useful")
library(useful)
kable(corner(count,r=15,c=8), booktabs=T, caption="Gene expression matrix")
count <- count[,c(1:3,7:9,14,4:6,10:12,13)]
names(count) <- rownames(Allcol)
count2 <- as.matrix(sapply(count,as.numeric))
count2[is.na(count2)] <- 0
row.names(count2) <- rownames(count)
dds <- DESeqDataSetFromMatrix(countData = count2,colData = Allcol,design = ~RESPONSE)
dds <- DESeq(dds)
deseqnorm <- as.data.frame(counts(dds,normalized=TRUE))

save(deseqnorm,file = "./DESeq_normcount.RData")

#TPM normalization 
count <- as.data.frame(fread("./HTseq.txt",fill = T,header = T))
geneLength <- as.data.frame(fread("./hg38_gene_length.txt",fill = T,header = F))
names(geneLength)[1] <- 'gene_id'
names(geneLength)[2] <- 'gene_length'

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
merge <- merge(ID,tpms,by="gene_id")
merge <- merge[!duplicated(merge$gene_name),]
merge <- merge[,-c(1,3)]
rownames(merge) <- merge[,1]
merge <- merge[,-1]
merge <- merge[,c(1:3,7:9,14,4:6,10:12,13)]
names(merge) <- rownames(Allcol)

save(merge,file = "./TPM_normcount.RData")

### 2021/03/09 ###               
setwd("/home/eunji/proj/2D3D/")

#TCGA PAAD_clinical , GDAC-mRNAseq raw-count 
PAAD<- as.data.frame(fread("./HYPOXIA/PAAD.clin.merged.txt",fill = T,header = F,stringsAsFactors = F))
PAAD <- as.data.frame(t(PAAD))
colnames(PAAD) <- PAAD[1,]
PAAD <- PAAD[-1,]
PAAD <- PAAD[,c(11:14,18,22,25,216:219,228,229,299,300,321,323,324,326,328:330,341,366,391,393,417,425:427,437)]

save(PAAD,file = "./HYPOXIA/clinical.RData")

PAAD <- PAAD %>%
  dplyr::filter(patient.drugs.drug.drug_name == "gemcitabine")

PAAD$patient.drugs.drug.measure_of_response <- ifelse(is.na(PAAD$patient.drugs.drug.measure_of_response),0,PAAD$patient.drugs.drug.measure_of_response) 
PAAD <- PAAD %>%
  dplyr::filter(patient.drugs.drug.measure_of_response != 0)
PAAD <- PAAD %>%
  dplyr::filter(patient.drugs.drug.measure_of_response == "clinical progressive disease" | patient.drugs.drug.measure_of_response == "complete response")
PAAD$response <- ifelse(PAAD$patient.drugs.drug.measure_of_response == "clinical progressive disease",'NO','YES')
response <- PAAD[,c(5,32)]
response <- data.frame(toupper(response$patient.bcr_patient_barcode),response$response)
names(response) <- c('SAMPLE','RESPONSE')

rownames(response) <- response[,1]
response <- data.frame(row.names = response$SAMPLE,RESPONSE=response$RESPONSE)

#response & hypoxia score 관계 
names(hypoxia)[1] <- 'SAMPLE'
reshy <- merge(response,hypoxia,by="SAMPLE")
names(reshy)[4] <- 'Winter_Hypoxia_score'
names(reshy)[5] <- 'Ragnum_hypoxia_score'
names(reshy)[6] <- "West_hypoxia_score"
  
  
  
install.packages("ggpubr")
install.packages("devtools")
library(ggpubr)
library(devtools)
p <- ggboxplot(reshy, x = "RESPONSE", y = "West_hypoxia_score",
               color = "RESPONSE", palette = "jco",
               bxp.errorbar = T,bxp.errorbar.width = 0.2,
               add = "jitter")+
  labs(title = 'HYPOXIA',
       caption = 'Data source: TCGA-PAAD',
       x='RESPONSE') 
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test",label.x = 0.7,label.y = 23)              
              
