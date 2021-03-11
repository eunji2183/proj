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
              
setwd("/home/eunji/miniconda3/tmp/2D3D")
library(data.table)
library(stringr)
library(ggpubr)
library(devtools)
library(DESeq2)
library(tibble)
library(GSVA)
library(limma)


#TCGA PAAD_clinical , GDAC-mRNAseq raw-count 
PAAD<- as.data.frame(fread("./HYPOXIA/PAAD.clin.merged.txt",fill = T,header = F,stringsAsFactors = F))
rownames(PAAD) <- PAAD[,1]
PAAD <- PAAD[,-1]
PAAD <- as.data.frame(t(PAAD))
PAAD <- PAAD[,c(11:14,18,22,25,216:219,228,229,299,300,321,323,324,326,328:330,341,366,391,393,417,425:427,437)]

save(PAAD,file = "./HYPOXIA/clinical.RData")

PAAD <- PAAD %>%
  dplyr::filter(patient.drugs.drug.drug_name == "gemcitabine")

PAAD <- PAAD %>%
  dplyr::filter(patient.drugs.drug.measure_of_response != "NA")
PAAD <- PAAD %>%
  dplyr::filter(patient.drugs.drug.measure_of_response == "clinical progressive disease" | patient.drugs.drug.measure_of_response == "complete response")
PAAD$response <- ifelse(PAAD$patient.drugs.drug.measure_of_response == "clinical progressive disease",'NO','YES')
response <- PAAD[,c(5,32)]
response <- data.frame(toupper(response$patient.bcr_patient_barcode),response$response)
names(response) <- c('SAMPLE','RESPONSE')

rownames(response) <- response[,1]
response <- data.frame(row.names = response$SAMPLE,RESPONSE=response$RESPONSE)



#response & hypoxia score 관계 
#TCGA hypoxia score 

hypoxia <- as.data.frame(fread("./HYPOXIA/TCGA_HYPOXIA_SCORE.txt",fill = T,header = T))
hypoxia <- hypoxia %>%
  dplyr::filter(tumour_type == "PAAD")
hypoxia$patient_id <- gsub('.','-',hypoxia$patient_id,fixed = T)
names(hypoxia)[1] <- 'SAMPLE'


reshy <- merge(response,hypoxia,by="SAMPLE")
names(reshy)[4] <- 'Winter_Hypoxia_score'
names(reshy)[5] <- 'Ragnum_hypoxia_score'
names(reshy)[6] <- "West_hypoxia_score"
  
  
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


#hypoxia 와 glutaminolysis 관계

count <- as.data.frame(fread("./HYPOXIA/count.txt",fill = T,header = T))
count <- count[!duplicated(count$`HYBRIDIZATION R`),]
rownames(count) <- count[,1]
count <- count[,-1]
count <- round(count)
tcount <- as.data.frame(t(count))
tcount$group <- as.numeric(as.character(substring(rownames(tcount),14,15)))
tcount$group <- ifelse(tcount$group < 5 , 'T','N')
tcount <- tcount %>%
  dplyr::filter(group == "T")
tcount <- tcount[,c(1:(length(colnames(tcount)) -1))]
rownames(tcount) <- str_sub(rownames(tcount),1,12)
count <- as.data.frame(t(tcount))
count <- count[-1,]
tcount <- as.data.frame(t(count))
tcount <- data.frame(SAMPLE=rownames(tcount),tcount)
rownames(tcount) <- NULL


tcount <- merge(hypoxia,tcount,by="SAMPLE")

hycol <- tcount[,1:11]
rownames(hycol)  <- hycol[,1]
hycol <- hycol[,-1]


tcount <- tcount[,-c(2:11)]
rownames(tcount) <- tcount[,1]
tcount <- tcount[,-1]
count <- as.data.frame(t(tcount))



rownames(tcount) <- tcount[,1]
tcount <- tcount[,-1]
tcount <- tcount[order(tcount$RESPONSE,decreasing = T),]

count <- as.data.frame(t(tcount))
count <- count[-1,]

glucol <- data.frame(row.names = rownames(tcount),RESPONSE=tcount[,1])
glucol <- glucol[order(glucol$RESPONSE,decreasing = T),]
glucol <- data.frame(row.names = colnames(count),RESPONSE=glucol)
glucol$RESPONSE <- as.factor(glucol$RESPONSE)

count2 <- as.matrix(sapply(count,as.numeric))
count2[is.na(count2)] <- 0
row.names(count2) <- rownames(count)
dds <- DESeqDataSetFromMatrix(countData = count2,colData = hycol,design = ~Ragnum_hypoxia_score_pan_cancer)
dds <- DESeq(dds)
normglu <- as.data.frame(counts(dds,normalized=TRUE))

glu <- read.table(file = "./2D3D/geneset/glutaminolysis.txt",sep = "\t",header = F)
glu <- glu$V1
glucount <- subset(normglu, (rownames(normglu) %in% glu))


p <- glucount

p <- apply(p,1,function(x){(x-mean(x))/sd(x)}) 

p <- merge(annocol,p,by="row.names")

p <- p %>% dplyr::filter(Row.names != "TCGA-2J-AABI")
p <- p %>% dplyr::filter(Row.names != "TCGA-2J-AABV")
p <- p %>% dplyr::filter(Row.names != "TCGA-FB-AAQ1")
p <- p %>% dplyr::filter(Row.names != "TCGA-IB-7654")
p <- p %>% dplyr::filter(Row.names != "TCGA-US-A779")
p <- p %>% dplyr::filter(Row.names != "TCGA-F2-6880")
p <- p %>% dplyr::filter(Row.names != "TCGA-HZ-7918")

s1 <- ggscatter(p, x = "hypoxia_score", y ="PPAT", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "hypoxia_score", ylab = "PPAT")
s2 <- ggscatter(p, x = "hypoxia_score", y ="GMPS", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "hypoxia_score", ylab = "GMPS")

s3 <- ggscatter(p, x = "hypoxia_score", y = "SLC1A5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "hypoxia_score", ylab = "SLC1A5")

par(mfrow=c(2,2))

grid.arrange(s1, s2, ncol=2)

ggscatter(p, x = "hypoxia_score", y = "SLC7A2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "hypoxia_score", ylab = "SLC7A2")


cor.test(tnorm$HDAC8, tnorm$HDAC6, 
         method = "pearson")


uns <- c("ASNS","CTPS2","GFPT2","GLS2","GLS","GLUL","GOT1","GPT2","GPT","MDH1","PFAS","PHGDH","SLC16A10","SLC6A15","SLC7A2")

p <- p %>% 
  dplyr::select(!(uns))

p <- p %>%
  dplyr::filter(hypoxia_score > 6 | hypoxia_score < -16)

annocol <- data.frame(SAMPLE = rownames(hycol),hypoxia_score=hycol$Ragnum_hypoxia_score_pan_cancer)
annocol <- annocol[order(annocol$hypoxia_score,decreasing = T),]
annocol <- data.frame(row.names = annocol$SAMPLE,hypoxia_score=annocol$hypoxia_score)




p <- data.frame(SAMPLE=rownames(p),as.data.frame(p))
rownames(p) <- NULL
p <- merge(annocol,p,by="SAMPLE")
p <- p[order(p$hypoxia_score,decreasing = T),]


annocol <- data.frame(row.names = p$Row.names,hypoxia_score=p$hypoxia_score)

p <- p[,-2]
rownames(p) <- p[,1]
p <- p[,-1]
p <- t(p)


p[p>5]=5
p[p<-5]=-5

paletteLength <- 5
myColor <- colorRampPalette(c('blue','white','firebrick'))(paletteLength)

myBreaks <- c(seq(min(p),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(p)/paletteLength, max(p), length.out=floor(paletteLength/2)))

ann_colors = list( hypoxia_score = c("white", "#D95F02"))

library(ComplexHeatmap)
ha = HeatmapAnnotation(hypoxia_score = anno_barplot(annocol$hypoxia_score, baseline = 0,
                                                    bar_width = 1,
                                                    gp = gpar(col = "white", fill = "black"), border = F,
                                                    height = unit(3, "cm")))


pheatmap(p, cluster_cols =F ,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=8,fontsize_row=8, cellwidth = 8,cellheight = 11,
         color=myColor,legend = T,breaks = myBreaks,top_annotation=ha,
         show_colnames = F,scale = F)


#response & hypoxia - glutaminolysis gene 
#coldata : response , norm: RESPONSE , glu, hypo, 

hypo <- read.table(file = "./2D3D/geneset/hypoxia.txt",header = F)
hypo <- data.frame(gene_name=hypo$V1,GeneGroup="Hypoxia")
glu <- read.table(file = "./2D3D/geneset/glutaminolysis.txt",sep = "\t",header = F)
glu <- data.frame(gene_name=glu$V1,GeneGroup="Glutaminolysis")
hypoglu <- rbind(hypo,glu)
hypoglu <- data.frame(row.names = hypoglu$gene_name,GeneGroup=hypoglu$GeneGroup)

hgcount <- merge(response,tcount,by="SAMPLE")

rownames(hgcount) <- hgcount[,1]
hgcount <- hgcount[,-1]
hgcount <- hgcount[order(hgcount$RESPONSE,decreasing = T),]

hgcol <- data.frame(row.names = rownames(hgcount),RESPONSE=hgcount[,1])
hgcol$RESPONSE <- as.factor(hgcol$RESPONSE)

hgcount <- as.data.frame(t(hgcount))
hgcount <- hgcount[-1,]

hgcount2 <- as.matrix(sapply(hgcount,as.numeric))
hgcount2[is.na(hgcount2)] <- 0
row.names(hgcount2) <- rownames(hgcount)
dds <- DESeqDataSetFromMatrix(countData = hgcount2,colData = hgcol,design = ~RESPONSE)
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="RESPONSE_YES_vs_NO", type="apeglm")
res <- as.data.frame(resLFC)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),
                 by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$log2FoldChange,decreasing = T),]


#GSEA rnk file 
hgGSEA <- resdata %>%
  dplyr::select(Row.names,log2FoldChange)
write.table(hgGSEA,file = "./hgGSEA.rnk",sep = "\t",row.names = F,col.names = F,quote = F)



normhg <- as.data.frame(counts(dds,normalized=TRUE))

hg <- rownames(hypoglu)
hgcount <- subset(normhg, (rownames(normhg) %in% hg))

h <- hgcount

h <- apply(h,1,function(x){(x-mean(x))/sd(x)}) 
h <- as.data.frame(t(h))

paletteLength <- 30
myColor <- colorRampPalette(c('blue','white','firebrick'))(paletteLength)

myBreaks <- c(seq(min(h),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(h)/paletteLength, max(h), length.out=floor(paletteLength/2)))
library(pheatmap)
pheatmap(h, cluster_cols =T ,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=8,fontsize_row=8, cellwidth = 8,cellheight = 11,
         color=myColor,legend = T,breaks = myBreaks,
         annotation_row = hypoglu,annotation_col = hgcol,
         show_colnames = F)



#ssGSEA 
geneSet <- read.csv2("./2D3D/geneset/hypoglu.gmt",header = F,sep = "\t")
geneSet <- geneSet[,-2]
geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
save(l,file = "./gene_set.Rdata")


dat <- as.matrix(normhg)
ssgsea<- gsva(dat, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ss <- as.data.frame(ssgsea)

ss <- apply(ss,1,function(x){(x-mean(x))/sd(x)})
ss <- ss[,-1]
ss <- ss[-c(38,36),]
ss <- ss[-c(2,6,14,17),]
ss <- ss[-c(14:17),]
ss <- ss[-c(19,21),]
ss <- ss[-14,]
ss <- ss[-23,]
ss <- ss[c(4,1,7,2,3,6,8,10,12,9,11,13,5,26,21,23,25,26,14,16,18,20,15,27,17,19,22,24,28),]



paletteLength <- 100
myColor <- colorRampPalette(c('black','white','yellow'))(paletteLength)

myBreaks <- c(seq(min(ss),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(ss)/paletteLength, max(ss), length.out=floor(paletteLength/2)))


ann_colors = list(
  RESPONSE = c(NO = "black", YES = "firebrick3")
)


pheatmap(ss, cluster_cols =T ,cluster_rows = F,
         clustering_distance_row = "correlation",
         fontsize_col=8,fontsize_row=8, cellwidth = 20,cellheight = 10,
         color=myColor,legend = T,breaks = myBreaks,annotation_row = hgcol,
         show_colnames = T,annotation_colors = ann_colors)


#scatterplot 

sp <- data.frame(row.names=rownames(ss),ss)
sp <- merge(sp,hgcol,by="row.names")


ggplot(sp,aes(HYPOXIA,GLUTAMINOLYSIS ,color=RESPONSE)) + 
  geom_point(size=5) + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0,1), legend.justification=c(0,1)) +
  stat_ellipse(type = "euclid")+
  labs(title="ssGSEA-TCGA",
       x="Hypoxia(NES)", y = "Glutaminolysis(NES)")
