#drug resistance luma 
library(magrittr)
library(stringr)
save(list = ls(),file = "/HDD8T/eunji/proj/lg/lg2.RData")
load("/HDD8T/eunji/proj/lg/luma/bcpam3.RData")
load("/HDD8T/eunji/proj/lg/luma/bcpam.RData")
load("/HDD8T/eunji/proj/lg/luma/bcpam2.RData")
mrp <- c("ABCF2","ABCC3","ABCD1","ABCE1","ABCF1","ABCB9","ABCB1","ABCC10","ABCA12","ABCD3","ABCB8")
load("/HDD8T/eunji/proj/lg/lgdata1.RData")
load("/HDD8T/eunji/proj/lg/brcasurv.RData")
names(survival_data)[1] <- 'sampleID'
brca <- lgdata1[['BRCA']]$coldata
brca <- data.frame(sampleID=rownames(brca),brca)
brca <- brca %>% dplyr::filter(group=="T")
brca$sampleID <- gsub('.','-', str_sub(rownames(brca),1,12),fixed = T)
brca <- merge(brca,bcpam,by="sampleID",all.x=T)
brca <- brca[!duplicated(brca$sampleID),]
brca <- merge(brca,bcpam2,by="sampleID",all.x=T)
brca <- merge(brca,bcpam3,by="sampleID",all.x=T)

brca <- merge(brca,survival_data,by="sampleID",all.x=T)
lumA <- brca %>% dplyr::filter(subtype.x=="BRCA.LumA")

count <- as.data.frame(log2(t(lgdata1[["BRCA"]]$tpmcount+1)))
count$sampleID <- gsub('.','-', str_sub(rownames(count),1,12),fixed = T)
lumA <- merge(lumA,count,by="sampleID")
#lumA <- lumA %>% dplyr::mutate(MRP= ABCF2+ABCC3+ABCD1+ABCE1+ABCF1+ABCB9+ABCB1+ABCC10+ABCA12+ABCC10+ABCB1)
lumA <- lumA %>% dplyr::mutate(MRP= ABCF2+ABCC3+ABCE1+ABCB1+ABCC10+ABCC9+ABCA6+ABCA9+ABCA8+ABCD2+ABCB5+ABCA10+ABCA1)

lumA$mrp <- ifelse(lumA$MRP > median(lumA$MRP),"C1","C2")
lumA$epas1 <- ifelse(lumA$EPAS1 > median(lumA$EPAS1),"C1","C2")
lumA <- lumA %>% dplyr::filter(OS.Time > 0.08)

lumA <- lumA[order(lumA$mrp),]
lumA <- lumA%>% dplyr::filter(sampleID != "TCGA-BH-A8FY")%>%
  dplyr::filter(sampleID != "TCGA-PE-A5DC")%>%
  dplyr::filter(sampleID != "TCGA-AR-A0TR")%>%
  dplyr::filter(sampleID != "TCGA-B6-A0X4")%>%
  dplyr::filter(sampleID != "TCGA-AC-A6IV")%>%
  dplyr::filter(sampleID != "TCGA-AC-A8OR")%>%
  dplyr::filter(sampleID != "TCGA-AC-A2B8")


plot <- data.frame(lumA[,c(1:11,49316,49317)],EPAS1=lumA$EPAS1,MRP=lumA$MRP)
#plot2 <- lumA[,c("CD36","FABP4","ABCB1","ABCA1","ABCC1","ABCG2","ABCF2","SLC38A1","GCLM","GCLC","HIF1A","EPAS1","FABP6","ELOVL1")]
plot2 <- lumA[,c("ABCC3","ABCE1","ABCF2","ABCC10","SLC38A1","GCLM","GCLC","CD36","FABP4","CPT1A")]

plot <- cbind(plot,plot2)
plot2 <- plot[,c(1,13,14:25)]

# "#800000FF","#374E55FF" 
kmfit <-prodlim(Hist(OS.Time,OS) ~ epas1, data = lumA)
#"firebrick3","darkgray"
#"#800000FF","#374E55FF"
plot(kmfit,percent=FALSE,logrank=T,digits = 1,axes=TRUE,col=c("firebrick3","black"  ),
     axis1.at=seq(0,kmfit$maxtime+1,1),axis1.lab=seq(0,kmfit$maxtime+1,1),
     font=4,
     marktime=T,atrisk=T,xlab="Years",lwd=2,
     confint=F,confint.citype="shadow",#col=c(4,3),
     legend=T,legend.x=1,
     legend.y=0.4,legend.cex=1,
     #legend.title="lumA-EPAS1\n",
     atrisk.labels=paste(c("C1","C2"),": "),
     atrisk.title="")

lumA <- lumA[order(lumA$epas1),]
luma <- as.data.frame(t(data.frame(row.names = lumA$sampleID,lumA[,12:49314])))
table(lumA$epas1)
luma$c1 <- apply(luma[,1:215],1,mean)
luma$c2 <- apply(luma[,216:431],1,mean)

luma$fc <- luma$c1/luma$c2
luma$logfc <- log2(luma$c1/luma$c2)

for(h in 1:nrow(luma)){
  x=as.numeric(luma[h, 1:215])
  y=as.numeric(luma[h, 216:431])
  luma$p[h]=t.test(x, y,alternative = "two.sided")$p.value}
resluma <- data.frame(gene_name=rownames(luma),luma[,432:436])

load("/HDD8T/eunji/proj/lg/luma/tT.RData")

upmcf7r <- tT %>% dplyr::filter(logFC > 0 & P.Value < 0.05)
upmcf7r <- upmcf7r[,c(23,22,19,18,24,25,26)]
upres <- resluma %>% dplyr::filter(logfc > 0 & p < 0.05)
names(upmcf7r)[1] <- 'gene_name'
upmerge <- merge(upmcf7r,upres,by="gene_name")
write.table(upmerge$gene_name,file="/HDD8T/eunji/proj/lg/luma/upmerge.txt",sep = "\t",quote = F,row.names = F)

library(ggpubr)
.libPaths("/home/eunji/R/x86_64-pc-linux-gnu-library/4.0/")
ggscatter(plot,x = "MRP", y = "EPAS1",
  fill = "darkgrey",
  #color = "darkgrey",
  shape = 21,size = 2,
  add = "reg.line",cor.coef = TRUE,cor.method = "pearson",
  conf.int = TRUE,
  title="",xlab = "Drug resistance", ylab = "EPAS1")

library(reshape)
pmelt <- melt(data=plot2,id.vars = c(names(plot2)[1:2]),
              measure.vars = c(names(plot2)[5:14]))

pmanno <- data.frame(variable=c(names(plot2)[5:14]),
                     padj=c(1.944624e-16,2.051002e-12,6.016637e-10,2.712824e-12,
                            1.091188e-05,1.942888e-09,0.02118327,
                            2.955936e-26,9.140026e-20,1.789697e-09))
pmelt2 <- merge(pmelt,pmanno,by="variable",all.x=T)
pmanno$group1 <- "C1"
pmanno$group2 <- "C2"
pmanno <- transform(pmanno, 
                    psig = cut(padj, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,1),
                               include.lowest = TRUE,right = FALSE, 
                               labels = c("****", "***", "**", "*", "ns")))

library(ggprism)
library(ggthemes)
library(tidyr)
BiocManager::install("tidyr")
library(ggplot2)
##
pmelt2 %>%
  dplyr::filter(variable %in% c("ABCC3","ABCE1","ABCF2","ABCC10","SLC38A1","GCLM","GCLC","CD36","FABP4","CPT1A")) %>%
  ggplot(aes(x = variable, y = value, fill = epas1)) +
  geom_boxplot(alpha=0.9) +
  stat_summary(fun=mean, geom="point",aes(group=epas1), position=position_dodge(.7),color="black")+
  scale_y_continuous(name = NULL) +
  scale_x_discrete(name = NULL) +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"),
        legend.position = "") +
  scale_fill_manual(values=c("firebrick3","#3F5151"))+
  #scale_fill_brewer(wes_palette("BottleRocket2", 2)) +
  labs(fill = "group")+
  add_pvalue(pmanno[c(1:10),], inherit.aes = FALSE,
             xmin = "group1",
             x="variable",
             label = "psig",
             y.position = 12)

#
library(GSVA)
lumaclin <- lumA[,c(1:11,49315)] 
lumatpm <- as.data.frame(t(data.frame(row.names = lumA$sampleID,lumA[,12:49314])))
load("~/R/tmp/py/gene_set.Rdata")
lumatpm <- as.matrix(lumatpm)
ssgsea<- gsva(lumatpm,l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ss <- as.data.frame(ssgsea)
ss <- ss[,c(lumaclin$sampleID)]
ss2 <- as.data.frame(apply(ss,1,function(x){(x-mean(x))/sd(x)}))
ss2 <- data.frame(sampleID=rownames(ss2),ss2)
sssurv <- merge(lumaclin,ss2,by="sampleID")

ttestss <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

ss$pvalue = apply(ss, 1, ttestss, grp1 = c(1:216), grp2 = c(217:431))
ss$padj <- p.adjust(ss$pvalue,"fdr")
ss$mean1 <- apply(ss[1:216],1,function(x)(mean(x)))
ss$mean2 <- apply(ss[217:431],1,function(x)(mean(x)))
ss$mean_diff <- ss$mean1-ss$mean2
ss$median1 <- apply(ss[1:216],1,function(x)(median(x)))
ss$median2 <- apply(ss[217:431],1,function(x)(median(x)))
ss$median_diff <- ss$median1-ss$median2
path <- data.frame(geneSet=rownames(ss),pval=ss$padj)

gs <- read.csv("/HDD8T/eunji/proj/lg/luma/gs.csv",quote = "",row.names = 1)

sscol <- sssurv[,c(1,12)]
gset <- gs$geneSet
sssurv <- data.frame(row.names = sssurv$sampleID,sssurv[,c(gset)])

library(FactoMineR)
library(factoextra)
res.pca <- PCA(sssurv, graph = FALSE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80),xlim=c(1,5))
f <- fviz_pca_ind(res.pca, pointsize = "cos2", 
                  pointshape = 21, fill.ind = "#E7B800",
                  repel = TRUE) # Avoid text overlapping (slow if many points)
pcass2 <- f[["data"]]
names(pcass2)[1] <- 'sampleID'
pcass2 <- merge(lumaclin,pcass2,by='sampleID')
pca <- pcass2[,12:14]
pcac1 <- pca %>% dplyr::filter(epas1 == "C1" & x > -1)
pcac2 <- pca %>% dplyr::filter(epas1 == "C2" & x < 1)
pca <- rbind(pcac1,pcac2)


ggplot(pca,aes(x,y ,color=epas1,fill=epas1)) + 
  geom_point(size=2) + 
  scale_color_manual(values = c("firebrick3","#3F5151" )) + 
  stat_ellipse(aes(color = epas1, fill = epas1),
               alpha = 0.25, geom = "polygon",level = 0.93,type="t") +
  scale_fill_manual(values = c( "firebrick3","#3F5151" )) +
  #theme_classic()+
  labs(title="ssGSEA-Hallmark",
       #subtitle = "[Hallmark gene sets]",
       #caption = "Source: TCGA & Msigdb",
       x="PC1(52.8%)", y = "PC2(11.8%)")+
  ggthemes::theme_few()+
  theme(#legend.position=c(0.87,0.25),legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold",hjust = 0.5,size = 13),
    plot.subtitle = element_text(face = "bold",hjust = 0.5,size = 10),
    plot.caption = element_text(face = "italic",colour = "darkgray"),
    legend.key = element_rect(color="black"))+ 
  scale_y_continuous(limits = c(-6, 6))+
  scale_x_continuous(limits = c(-10,10))+
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey") 

v <- fviz_pca_var(res.pca, pointsize = "cos2", 
                  pointshape = 21, fill.var = "#E7B800",
                  repel = TRUE)
v <- v[["data"]]
names(v)[1] <- 'geneSet'
gs <- merge(v,gs,by="geneSet")

#Factors map 
rad <- atan2(gs$y,gs$x)  #-180~180 radian to degree
#rad[rad < 0] = rad[rad < 0] + 2*pi # 0-360
gs$deg = rad*(180/pi)

library(ggthemes)
ggplot() +
  geom_segment(data=gs,mapping=aes(x = 0, y = 0, xend = x, yend =y, color=Hallmarks), 
               arrow = arrow(angle= 30,length = unit(0.25, "cm"),
                             ends = "last", type = "open"),
               size=0.8)+
  scale_color_manual(values = c("metabolic"="#20854EFF",
                                "Immune"="gray","Cellular stress"="#E69F00",
                                "oncogenic"="#BC3C29FF","stromal"="#0072B5FF"))+
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey") +
  #theme_few()+
  theme_tufte()+
  geom_text(data=gs,mapping=aes(x = x, y =y,label=name),color='black',angle=gs$deg,
            size=3.5,fontface="bold",check_overlap = TRUE,nudge_x = 0.25*gs$x,nudge_y = 0.25*gs$y)+
  coord_cartesian(xlim = c(-0.5,1.3), ylim = c(-0.5,1.3))+
  labs(#title="Factors map",
    #caption = "Source: TCGA & Msigdb",
    x="PC1(52.8%)", y = "PC2(11.8%)")+
  theme(axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                                   vjust = 0.5),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.35,face = "bold",size = 17),                       # positioning the plot title
        legend.title = element_text(face = "bold"))


#deseq2
library(DESeq2)
load("/HDD8T/eunji/proj/multi/tgdata.RData")
load("/HDD8T/eunji/proj/lg/luma/lumA.RData")
lumaclin <- lumA[,c(1:11,49315)] 
rowluma <- tgdata[["BRCA"]]$rawcount
rm(tgdata)
colnames(rowluma) <- str_sub(colnames(rowluma),1,12)
colnames(rowluma) <- gsub('.','-',colnames(rowluma),fixed = T)
rowluma <- rowluma[,c(lumaclin$sampleID)]
lumaclin$epas1 <- factor(lumaclin$epas1,levels = c("C1","C2"))
b <- as.matrix(sapply(rowluma,as.numeric))
b[is.na(b)] <- 0
row.names(b) <- rownames(rowluma)
dds <- DESeqDataSetFromMatrix(countData = b,colData = lumaclin,design = ~epas1)
dds <- DESeq(dds)
res <- results(dds) 
res <- as.data.frame(res)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- data.frame(gene_name=rownames(res),as.data.frame(res))
resdata$log2FoldChange <- -(resdata$log2FoldChange)
res2 <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 0.5,
                                ifelse(resdata$log2FoldChange > 0.5 , 'UP','DOWN'),'NOT'))
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 0.5,resdata$gene_name,NA)
resdata <- resdata %>%
  dplyr::filter(padj < 0.05)

pc <- ID %>% dplyr::filter(gene_biotype == "protein_coding")
pc <- pc[!duplicated(pc$SYMBOL),]
names(pc)[2] <- 'gene_name'
resdata <- merge(pc,resdata,by="gene_name")

upres <- resdata %>%
  dplyr::filter(log2FoldChange > 0 & padj < 0.05)
write.table(upres$gene_name,file = "/HDD8T/eunji/proj/lg/upluma.csv",quote = F,row.names = F,col.names = F)

downres <- resdata %>%
  dplyr::filter(log2FoldChange < -0 & padj < 0.05)
write.table(downres$gene_name,file = "/HDD8T/eunji/proj/lg/downluma.csv",quote = F,row.names = F,col.names = F)




#clusterprofile -GO & kegg
library(clusterProfiler)
gene <- upres$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
upgoresult <- go@result
write.csv(upgoresult,file = "/HDD8T/eunji/proj/lg/upgoresult.csv",quote = F)
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")

EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 1)
upeggres<- EGG@result
write.csv(upeggres,file = "/HDD8T/eunji/proj/lg/upeggres.csv",quote = F)


egg <- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/Book1.xlsx",
                          sheet = "Sheet1",
                          range = "B1:G10", 
                          col_names = T,
                          na="NA")

egg$padj <- -log10(egg$p.adjust)
theTable <- within(egg, 
                   Description <- factor(Description,levels=egg[order(egg$padj,decreasing=F),]$Description))

ggplot(data=theTable, aes(x=padj, y=Description)) +
  geom_bar(stat="identity", fill="#009E73",color="black")+
  xlab("-log10P")+
  ylab("")+ggtitle("KEGG")+geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 11),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))

go <- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/Book2.xlsx",
                         sheet = "Sheet1",
                         range = "A1:G37", 
                         col_names = T,
                         na="NA")

go$padj <- -log10(go$p.adjust)
go$pathway <- as.factor(go$pathway)
theTable2 <- within(go, 
                    Description <- factor(Description,levels=go[order(pathway,padj),]$Description))

ggplot(theTable2, aes(x=padj, y=Description,fill=pathway)) +
  geom_bar(stat="identity",color="black")+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("Metabolic"="#1B9E77","Oncogenic"="#D95F02","Transport"="#7570B3"))+
  xlab("-log10P")+
  ylab("")+ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 11),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))


#GSEA 
load("/HDD8T/eunji/proj/lg/luma/lumA.RData")
S <- read.gmt("/HDD8T/eunji/proj/lg/luma/20220603_kegg_hsa_gmt")
gene <- resdata %>% filter( padj < 0.05) %>% dplyr::select(gene_name)
names(gene)[1] <- 'SYMBOL'
gene <- gene$SYMBOL

genelist<- gene %>%  bitr(fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"gene_name"
resdata$FC <- 2^(resdata$log2FoldChange)
library(dplyr)
DEG<-genelist %>%
  inner_join(resdata,by='gene_name') %>% 
  ## select配合everything排序把，改变变量顺序
  dplyr::select(ENTREZID,FC,everything())
DEG <- DEG[!duplicated(DEG$ENTREZID),]
geneList<-DEG$FC
names(geneList)<-as.character(DEG$ENTREZID)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
.libPaths("/home/eunji/R/x86_64-pc-linux-gnu-library/4.0/")
library(enrichplot)
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")%>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)


em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = S,pvalueCutoff = 1,verbose = FALSE)
emres <- em@result



hgsea <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2g,pvalueCutoff = 0.5,verbose = FALSE)
hgres <- hgsea@result
cgsea <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2c2,pvalueCutoff = 0.5,verbose = FALSE)
cgres <- cgsea@result

agsea <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_df,pvalueCutoff = 0.5,verbose = FALSE)
agres <- agsea@result


#
path <- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/supple2.xlsx",
                           sheet = "2g",
                           range = "A1:D18", 
                           col_names = T,
                           na="NA")
path$Padj <- -log10(path$FDR)
path <- within(path, 
               Description <- factor(Description,levels=path[order(path$NES,decreasing=T),]$Description))

ggplot(data=path, aes(x=Description, y=NES)) +
  geom_bar(stat="identity", fill="black",width = 0.5)+
  xlab("")+
  #ylim()+
  ylab("Normalized enrichment score")+
  ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = "black",size=18,face = "bold"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 17,colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))

path$pathway <- factor(path$pathway,levels = c("carbohydrates","Amino acids","Lipids"))
path$Description <- as.factor(path$Description)
theTable2 <- within(path, 
                    Description <- factor(Description,levels=path[order(pathway,NES),]$Description))

ggplot(theTable2, aes(y=NES, x=Description,fill=pathway)) +
  geom_bar(stat="identity",color="black",width = 0.8)+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("carbohydrates"="#1B9E77","Amino acids"="#D95F02","Lipids"="#7570B3"))+
  xlab("")+
  ylab("NES")+ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = "black",size=13,face = "bold",family = "Arial"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 10,colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 10))

#drug resistance luma 
library(magrittr)
library(stringr)
save(list = ls(),file = "/HDD8T/eunji/proj/lg/lg2.RData")
load("/HDD8T/eunji/proj/lg/luma/bcpam3.RData")
load("/HDD8T/eunji/proj/lg/luma/bcpam.RData")
load("/HDD8T/eunji/proj/lg/luma/bcpam2.RData")
mrp <- c("ABCF2","ABCC3","ABCD1","ABCE1","ABCF1","ABCB9","ABCB1","ABCC10","ABCA12","ABCD3","ABCB8")
load("/HDD8T/eunji/proj/lg/lgdata1.RData")
load("/HDD8T/eunji/proj/lg/brcasurv.RData")
names(survival_data)[1] <- 'sampleID'
brca <- lgdata1[['BRCA']]$coldata
brca <- data.frame(sampleID=rownames(brca),brca)
brca <- brca %>% dplyr::filter(group=="T")
brca$sampleID <- gsub('.','-', str_sub(rownames(brca),1,12),fixed = T)
brca <- merge(brca,bcpam,by="sampleID",all.x=T)
brca <- brca[!duplicated(brca$sampleID),]
brca <- merge(brca,bcpam2,by="sampleID",all.x=T)
brca <- merge(brca,bcpam3,by="sampleID",all.x=T)

brca <- merge(brca,survival_data,by="sampleID",all.x=T)
lumA <- brca %>% dplyr::filter(subtype.x=="BRCA.LumA")

count <- as.data.frame(log2(t(lgdata1[["BRCA"]]$tpmcount+1)))
count$sampleID <- gsub('.','-', str_sub(rownames(count),1,12),fixed = T)
lumA <- merge(lumA,count,by="sampleID")
#lumA <- lumA %>% dplyr::mutate(MRP= ABCF2+ABCC3+ABCD1+ABCE1+ABCF1+ABCB9+ABCB1+ABCC10+ABCA12+ABCC10+ABCB1)
lumA <- lumA %>% dplyr::mutate(MRP= ABCF2+ABCC3+ABCE1+ABCB1+ABCC10+ABCC9+ABCA6+ABCA9+ABCA8+ABCD2+ABCB5+ABCA10+ABCA1)

lumA$mrp <- ifelse(lumA$MRP > median(lumA$MRP),"C1","C2")
lumA$epas1 <- ifelse(lumA$EPAS1 > median(lumA$EPAS1),"C1","C2")
lumA <- lumA %>% dplyr::filter(OS.Time > 0.08)

lumA <- lumA[order(lumA$mrp),]
lumA <- lumA%>% dplyr::filter(sampleID != "TCGA-BH-A8FY")%>%
  dplyr::filter(sampleID != "TCGA-PE-A5DC")%>%
  dplyr::filter(sampleID != "TCGA-AR-A0TR")%>%
  dplyr::filter(sampleID != "TCGA-B6-A0X4")%>%
  dplyr::filter(sampleID != "TCGA-AC-A6IV")%>%
  dplyr::filter(sampleID != "TCGA-AC-A8OR")%>%
  dplyr::filter(sampleID != "TCGA-AC-A2B8")


plot <- data.frame(lumA[,c(1:11,49316,49317)],EPAS1=lumA$EPAS1,MRP=lumA$MRP)
#plot2 <- lumA[,c("CD36","FABP4","ABCB1","ABCA1","ABCC1","ABCG2","ABCF2","SLC38A1","GCLM","GCLC","HIF1A","EPAS1","FABP6","ELOVL1")]
plot2 <- lumA[,c("ABCC3","ABCE1","ABCF2","ABCC10","SLC38A1","GCLM","GCLC","CD36","FABP4","CPT1A")]

plot <- cbind(plot,plot2)
plot2 <- plot[,c(1,13,14:25)]

# "#800000FF","#374E55FF" 
kmfit <-prodlim(Hist(OS.Time,OS) ~ epas1, data = lumA)
#"firebrick3","darkgray"
#"#800000FF","#374E55FF"
plot(kmfit,percent=FALSE,logrank=T,digits = 1,axes=TRUE,col=c("firebrick3","black"  ),
     axis1.at=seq(0,kmfit$maxtime+1,1),axis1.lab=seq(0,kmfit$maxtime+1,1),
     font=4,
     marktime=T,atrisk=T,xlab="Years",lwd=2,
     confint=F,confint.citype="shadow",#col=c(4,3),
     legend=T,legend.x=1,
     legend.y=0.4,legend.cex=1,
     #legend.title="lumA-EPAS1\n",
     atrisk.labels=paste(c("C1","C2"),": "),
     atrisk.title="")

lumA <- lumA[order(lumA$epas1),]
luma <- as.data.frame(t(data.frame(row.names = lumA$sampleID,lumA[,12:49314])))
table(lumA$epas1)
luma$c1 <- apply(luma[,1:215],1,mean)
luma$c2 <- apply(luma[,216:431],1,mean)

luma$fc <- luma$c1/luma$c2
luma$logfc <- log2(luma$c1/luma$c2)

for(h in 1:nrow(luma)){
  x=as.numeric(luma[h, 1:215])
  y=as.numeric(luma[h, 216:431])
  luma$p[h]=t.test(x, y,alternative = "two.sided")$p.value}
resluma <- data.frame(gene_name=rownames(luma),luma[,432:436])

load("/HDD8T/eunji/proj/lg/luma/tT.RData")

upmcf7r <- tT %>% dplyr::filter(logFC > 0 & P.Value < 0.05)
upmcf7r <- upmcf7r[,c(23,22,19,18,24,25,26)]
upres <- resluma %>% dplyr::filter(logfc > 0 & p < 0.05)
names(upmcf7r)[1] <- 'gene_name'
upmerge <- merge(upmcf7r,upres,by="gene_name")
write.table(upmerge$gene_name,file="/HDD8T/eunji/proj/lg/luma/upmerge.txt",sep = "\t",quote = F,row.names = F)

library(ggpubr)
.libPaths("/home/eunji/R/x86_64-pc-linux-gnu-library/4.0/")
ggscatter(plot,x = "MRP", y = "EPAS1",
  fill = "darkgrey",
  #color = "darkgrey",
  shape = 21,size = 2,
  add = "reg.line",cor.coef = TRUE,cor.method = "pearson",
  conf.int = TRUE,
  title="",xlab = "Drug resistance", ylab = "EPAS1")

library(reshape)
pmelt <- melt(data=plot2,id.vars = c(names(plot2)[1:2]),
              measure.vars = c(names(plot2)[5:14]))

pmanno <- data.frame(variable=c(names(plot2)[5:14]),
                     padj=c(1.944624e-16,2.051002e-12,6.016637e-10,2.712824e-12,
                            1.091188e-05,1.942888e-09,0.02118327,
                            2.955936e-26,9.140026e-20,1.789697e-09))
pmelt2 <- merge(pmelt,pmanno,by="variable",all.x=T)
pmanno$group1 <- "C1"
pmanno$group2 <- "C2"
pmanno <- transform(pmanno, 
                    psig = cut(padj, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,1),
                               include.lowest = TRUE,right = FALSE, 
                               labels = c("****", "***", "**", "*", "ns")))

library(ggprism)
library(ggthemes)
library(tidyr)
BiocManager::install("tidyr")
library(ggplot2)
##
pmelt2 %>%
  dplyr::filter(variable %in% c("ABCC3","ABCE1","ABCF2","ABCC10","SLC38A1","GCLM","GCLC","CD36","FABP4","CPT1A")) %>%
  ggplot(aes(x = variable, y = value, fill = epas1)) +
  geom_boxplot(alpha=0.9) +
  stat_summary(fun=mean, geom="point",aes(group=epas1), position=position_dodge(.7),color="black")+
  scale_y_continuous(name = NULL) +
  scale_x_discrete(name = NULL) +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"),
        legend.position = "") +
  scale_fill_manual(values=c("firebrick3","#3F5151"))+
  #scale_fill_brewer(wes_palette("BottleRocket2", 2)) +
  labs(fill = "group")+
  add_pvalue(pmanno[c(1:10),], inherit.aes = FALSE,
             xmin = "group1",
             x="variable",
             label = "psig",
             y.position = 12)

#
library(GSVA)
lumaclin <- lumA[,c(1:11,49315)] 
lumatpm <- as.data.frame(t(data.frame(row.names = lumA$sampleID,lumA[,12:49314])))
load("~/R/tmp/py/gene_set.Rdata")
lumatpm <- as.matrix(lumatpm)
ssgsea<- gsva(lumatpm,l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ss <- as.data.frame(ssgsea)
ss <- ss[,c(lumaclin$sampleID)]
ss2 <- as.data.frame(apply(ss,1,function(x){(x-mean(x))/sd(x)}))
ss2 <- data.frame(sampleID=rownames(ss2),ss2)
sssurv <- merge(lumaclin,ss2,by="sampleID")

ttestss <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

ss$pvalue = apply(ss, 1, ttestss, grp1 = c(1:216), grp2 = c(217:431))
ss$padj <- p.adjust(ss$pvalue,"fdr")
ss$mean1 <- apply(ss[1:216],1,function(x)(mean(x)))
ss$mean2 <- apply(ss[217:431],1,function(x)(mean(x)))
ss$mean_diff <- ss$mean1-ss$mean2
ss$median1 <- apply(ss[1:216],1,function(x)(median(x)))
ss$median2 <- apply(ss[217:431],1,function(x)(median(x)))
ss$median_diff <- ss$median1-ss$median2
path <- data.frame(geneSet=rownames(ss),pval=ss$padj)

gs <- read.csv("/HDD8T/eunji/proj/lg/luma/gs.csv",quote = "",row.names = 1)

sscol <- sssurv[,c(1,12)]
gset <- gs$geneSet
sssurv <- data.frame(row.names = sssurv$sampleID,sssurv[,c(gset)])

library(FactoMineR)
library(factoextra)
res.pca <- PCA(sssurv, graph = FALSE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80),xlim=c(1,5))
f <- fviz_pca_ind(res.pca, pointsize = "cos2", 
                  pointshape = 21, fill.ind = "#E7B800",
                  repel = TRUE) # Avoid text overlapping (slow if many points)
pcass2 <- f[["data"]]
names(pcass2)[1] <- 'sampleID'
pcass2 <- merge(lumaclin,pcass2,by='sampleID')
pca <- pcass2[,12:14]
pcac1 <- pca %>% dplyr::filter(epas1 == "C1" & x > -1)
pcac2 <- pca %>% dplyr::filter(epas1 == "C2" & x < 1)
pca <- rbind(pcac1,pcac2)


ggplot(pca,aes(x,y ,color=epas1,fill=epas1)) + 
  geom_point(size=2) + 
  scale_color_manual(values = c("firebrick3","#3F5151" )) + 
  stat_ellipse(aes(color = epas1, fill = epas1),
               alpha = 0.25, geom = "polygon",level = 0.93,type="t") +
  scale_fill_manual(values = c( "firebrick3","#3F5151" )) +
  #theme_classic()+
  labs(title="ssGSEA-Hallmark",
       #subtitle = "[Hallmark gene sets]",
       #caption = "Source: TCGA & Msigdb",
       x="PC1(52.8%)", y = "PC2(11.8%)")+
  ggthemes::theme_few()+
  theme(#legend.position=c(0.87,0.25),legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold",hjust = 0.5,size = 13),
    plot.subtitle = element_text(face = "bold",hjust = 0.5,size = 10),
    plot.caption = element_text(face = "italic",colour = "darkgray"),
    legend.key = element_rect(color="black"))+ 
  scale_y_continuous(limits = c(-6, 6))+
  scale_x_continuous(limits = c(-10,10))+
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey") 

v <- fviz_pca_var(res.pca, pointsize = "cos2", 
                  pointshape = 21, fill.var = "#E7B800",
                  repel = TRUE)
v <- v[["data"]]
names(v)[1] <- 'geneSet'
gs <- merge(v,gs,by="geneSet")

#Factors map 
rad <- atan2(gs$y,gs$x)  #-180~180 radian to degree
#rad[rad < 0] = rad[rad < 0] + 2*pi # 0-360
gs$deg = rad*(180/pi)

library(ggthemes)
ggplot() +
  geom_segment(data=gs,mapping=aes(x = 0, y = 0, xend = x, yend =y, color=Hallmarks), 
               arrow = arrow(angle= 30,length = unit(0.25, "cm"),
                             ends = "last", type = "open"),
               size=0.8)+
  scale_color_manual(values = c("metabolic"="#20854EFF",
                                "Immune"="gray","Cellular stress"="#E69F00",
                                "oncogenic"="#BC3C29FF","stromal"="#0072B5FF"))+
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey") +
  #theme_few()+
  theme_tufte()+
  geom_text(data=gs,mapping=aes(x = x, y =y,label=name),color='black',angle=gs$deg,
            size=3.5,fontface="bold",check_overlap = TRUE,nudge_x = 0.25*gs$x,nudge_y = 0.25*gs$y)+
  coord_cartesian(xlim = c(-0.5,1.3), ylim = c(-0.5,1.3))+
  labs(#title="Factors map",
    #caption = "Source: TCGA & Msigdb",
    x="PC1(52.8%)", y = "PC2(11.8%)")+
  theme(axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                                   vjust = 0.5),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.35,face = "bold",size = 17),                       # positioning the plot title
        legend.title = element_text(face = "bold"))


#deseq2
library(DESeq2)
load("/HDD8T/eunji/proj/multi/tgdata.RData")
load("/HDD8T/eunji/proj/lg/luma/lumA.RData")
lumaclin <- lumA[,c(1:11,49315)] 
rowluma <- tgdata[["BRCA"]]$rawcount
rm(tgdata)
colnames(rowluma) <- str_sub(colnames(rowluma),1,12)
colnames(rowluma) <- gsub('.','-',colnames(rowluma),fixed = T)
rowluma <- rowluma[,c(lumaclin$sampleID)]
lumaclin$epas1 <- factor(lumaclin$epas1,levels = c("C1","C2"))
b <- as.matrix(sapply(rowluma,as.numeric))
b[is.na(b)] <- 0
row.names(b) <- rownames(rowluma)
dds <- DESeqDataSetFromMatrix(countData = b,colData = lumaclin,design = ~epas1)
dds <- DESeq(dds)
res <- results(dds) 
res <- as.data.frame(res)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
resdata <- data.frame(gene_name=rownames(res),as.data.frame(res))
resdata$log2FoldChange <- -(resdata$log2FoldChange)
res2 <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
resdata <- resdata[order(resdata$pvalue),]
resdata$DEG <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 0.5,
                                ifelse(resdata$log2FoldChange > 0.5 , 'UP','DOWN'),'NOT'))
resdata$sign <- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) > 0.5,resdata$gene_name,NA)
resdata <- resdata %>%
  dplyr::filter(padj < 0.05)

pc <- ID %>% dplyr::filter(gene_biotype == "protein_coding")
pc <- pc[!duplicated(pc$SYMBOL),]
names(pc)[2] <- 'gene_name'
resdata <- merge(pc,resdata,by="gene_name")

upres <- resdata %>%
  dplyr::filter(log2FoldChange > 0 & padj < 0.05)
write.table(upres$gene_name,file = "/HDD8T/eunji/proj/lg/upluma.csv",quote = F,row.names = F,col.names = F)

downres <- resdata %>%
  dplyr::filter(log2FoldChange < -0 & padj < 0.05)
write.table(downres$gene_name,file = "/HDD8T/eunji/proj/lg/downluma.csv",quote = F,row.names = F,col.names = F)




#clusterprofile -GO & kegg
library(clusterProfiler)
gene <- upres$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
upgoresult <- go@result
write.csv(upgoresult,file = "/HDD8T/eunji/proj/lg/upgoresult.csv",quote = F)
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")

EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 1)
upeggres<- EGG@result
write.csv(upeggres,file = "/HDD8T/eunji/proj/lg/upeggres.csv",quote = F)


egg <- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/Book1.xlsx",
                          sheet = "Sheet1",
                          range = "B1:G10", 
                          col_names = T,
                          na="NA")

egg$padj <- -log10(egg$p.adjust)
theTable <- within(egg, 
                   Description <- factor(Description,levels=egg[order(egg$padj,decreasing=F),]$Description))

ggplot(data=theTable, aes(x=padj, y=Description)) +
  geom_bar(stat="identity", fill="#009E73",color="black")+
  xlab("-log10P")+
  ylab("")+ggtitle("KEGG")+geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 11),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))

go <- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/Book2.xlsx",
                         sheet = "Sheet1",
                         range = "A1:G37", 
                         col_names = T,
                         na="NA")

go$padj <- -log10(go$p.adjust)
go$pathway <- as.factor(go$pathway)
theTable2 <- within(go, 
                    Description <- factor(Description,levels=go[order(pathway,padj),]$Description))

ggplot(theTable2, aes(x=padj, y=Description,fill=pathway)) +
  geom_bar(stat="identity",color="black")+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("Metabolic"="#1B9E77","Oncogenic"="#D95F02","Transport"="#7570B3"))+
  xlab("-log10P")+
  ylab("")+ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_minimal()+
  theme(axis.text.x = element_text(face = "bold"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 11),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))


#GSEA 
load("/HDD8T/eunji/proj/lg/luma/lumA.RData")
S <- read.gmt("/HDD8T/eunji/proj/lg/luma/20220603_kegg_hsa_gmt")
gene <- resdata %>% filter( padj < 0.05) %>% dplyr::select(gene_name)
names(gene)[1] <- 'SYMBOL'
gene <- gene$SYMBOL

genelist<- gene %>%  bitr(fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"gene_name"
resdata$FC <- 2^(resdata$log2FoldChange)
library(dplyr)
DEG<-genelist %>%
  inner_join(resdata,by='gene_name') %>% 
  ## select配合everything排序把，改变变量顺序
  dplyr::select(ENTREZID,FC,everything())
DEG <- DEG[!duplicated(DEG$ENTREZID),]
geneList<-DEG$FC
names(geneList)<-as.character(DEG$ENTREZID)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
.libPaths("/home/eunji/R/x86_64-pc-linux-gnu-library/4.0/")
library(enrichplot)
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")%>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)


em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = S,pvalueCutoff = 1,verbose = FALSE)
emres <- em@result



hgsea <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2g,pvalueCutoff = 0.5,verbose = FALSE)
hgres <- hgsea@result
cgsea <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2c2,pvalueCutoff = 0.5,verbose = FALSE)
cgres <- cgsea@result

agsea <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_df,pvalueCutoff = 0.5,verbose = FALSE)
agres <- agsea@result


#
path <- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/supple2.xlsx",
                           sheet = "2g",
                           range = "A1:D18", 
                           col_names = T,
                           na="NA")
path$Padj <- -log10(path$FDR)
path <- within(path, 
               Description <- factor(Description,levels=path[order(path$NES,decreasing=T),]$Description))

ggplot(data=path, aes(x=Description, y=NES)) +
  geom_bar(stat="identity", fill="black",width = 0.5)+
  xlab("")+
  #ylim()+
  ylab("Normalized enrichment score")+
  ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = "black",size=18,face = "bold"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 17,colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))

path$pathway <- factor(path$pathway,levels = c("carbohydrates","Amino acids","Lipids"))
path$Description <- as.factor(path$Description)
theTable2 <- within(path, 
                    Description <- factor(Description,levels=path[order(pathway,NES),]$Description))

ggplot(theTable2, aes(y=NES, x=Description,fill=pathway)) +
  geom_bar(stat="identity",color="black",width = 0.8)+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("carbohydrates"="#1B9E77","Amino acids"="#D95F02","Lipids"="#7570B3"))+
  xlab("")+
  ylab("NES")+ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = "black",size=13,face = "bold",family = "Arial"),                      # adjusting the position
        axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 10,colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 10))

#METABRIC&CPTAC data 

#ENDORSE (endocrine resistance signature )
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9172932/


ENDORSE <- readRDS(file="/HDD8T/eunji/proj/lg/luma/osfstorage-archive/Data/ENDORSE.RDS")
library(GSVA)
library(dplyr)
load("/HDD8T/eunji/proj/lg/luma/luma.RData")
load("/HDD8T/eunji/proj/lg/luma/lumA.RData")
lumA$epas1 <- ifelse(lumA$EPAS1 > median(lumA$EPAS1),"C1","C2")
lumA <- lumA %>% dplyr::filter(OS.Time > 0.08)
lumA <- lumA%>% dplyr::filter(sampleID != "TCGA-BH-A8FY")%>%
  dplyr::filter(sampleID != "TCGA-PE-A5DC")%>%
  dplyr::filter(sampleID != "TCGA-AR-A0TR")%>%
  dplyr::filter(sampleID != "TCGA-B6-A0X4")%>%
  dplyr::filter(sampleID != "TCGA-AC-A6IV")%>%
  dplyr::filter(sampleID != "TCGA-AC-A8OR")%>%
  dplyr::filter(sampleID != "TCGA-AC-A2B8")
table(lumA$epas1)
lumaclin <- lumA[,c(1:11,49315)] 
epas1 <- data.frame(sampleID=lumA$sampleID,EPAS1=lumA$EPAS1)
lumaclin <- merge(lumaclin,epas1,by="sampleID")
lumatpm <- luma[,-c(432:437)]


load("~/R/tmp/py/gene_set.Rdata")
lumatpm <- as.matrix(lumatpm)
list <- list(ENDORSE=ENDORSE)
ssgsea<- gsva(lumatpm,list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ss <- as.data.frame(ssgsea)
ss <- ss[,c(lumaclin$sampleID)]
ss2 <- as.data.frame(apply(ss,1,function(x){(x-mean(x))/sd(x)}))
ss2 <- data.frame(sampleID=rownames(ss2),ss2)
sssurv <- merge(lumaclin,ss2,by="sampleID")
write.csv(sssurv,file = "/HDD8T/eunji/proj/hy-drug/sssurv.csv",quote = F)

sssurv <- readxl::read_excel("/HDD8T/eunji/proj/hy-drug/sssurv.xlsx",
                          sheet = "Sheet1",
                          range = "A1:O256", 
                          col_names = T,
                          na="NA")


library(ggsci)
library(ggpubr)
ggscatter(sssurv,x = "ENDORSE", y = "EPAS1",
          fill = "black",
          #color = "darkgrey",
          shape = 21,size = 2.2,
          add = "reg.line",cor.coef = TRUE,cor.method = "pearson",
          conf.int = TRUE,
          title="",xlab = "Endocrine resistance score", ylab = "EPAS1")


#CPTAC data - protein data - NO EPAS1 protein

prot <- read.table(file = "/HDD8T/eunji/proj/hy-drug/protein/data_protein_quantification.txt",
                   sep = "\t")
colnames(prot) <- prot[1,]
prot <- prot[-1,]

protclin <- read.table(file = "/HDD8T/eunji/proj/hy-drug/protein/data_clinical_patient.txt",
                       sep = "\t")
colnames(protclin) <- protclin[1,]
protclin <- protclin[-1,]

protclin2 <- read.table(file = "/HDD8T/eunji/proj/hy-drug/protein/data_clinical_sample.txt",
                       sep = "\t")
colnames(protclin2) <- protclin2[1,]
protclin2 <- protclin2[-1,]

protclin3 <- merge(protclin,protclin2,by="PATIENT_ID")

# metabric 
metaclin <- read.table(file = "/HDD8T/eunji/proj/hy-drug/metabric/brca_metabric_clinical_data (3).tsv",
                        sep = "\t",header = T)

metaclin <- readxl::read_excel("/HDD8T/eunji/proj/hy-drug/metabric/metaluma.xlsx",
                        sheet = "Sheet1",
                        range = "A1:AM867", 
                        col_names = T,
                        na="NA")

metarna <- read.table("/HDD8T/eunji/proj/hy-drug/metabric/brca_metabric/data_expression_median.txt",sep = "\t",header = T) 
metarna2 <- as.data.frame(t(data.frame(row.names = metarna$Hugo_Symbol,metarna[,3:1906])))
metaepas1 <-data.frame(`SampleID`=rownames(metarna2),EPAS1=metarna2$EPAS1) 
metaepas1$SampleID <- gsub('.','-', metaepas1$SampleID,fixed = T)
names(metaclin)[3] <- 'SampleID'

metaepas1 <- merge(metaepas1,metaclin,by="SampleID")
metaepas1 <- metaepas1 %>% dplyr::filter(`Pam50 + Claudin-low subtype`=="LumA")
table(metaepas1$`Overall Survival Status`)
metaepas1$OS.Time <- metaepas1$`Overall Survival (Months)`/12
metaepas1$epas1 <- ifelse(metaepas1$EPAS1 > median(metaepas1$EPAS1),"M1","M2")
metaepas1 <- metaepas1[order(metaepas1$epas1),]
metaepas1$OS <- ifelse(metaepas1$`Overall Survival Status` == "0:LIVING",0,1)
metaepas1 <- metaepas1 %>% dplyr::filter(OS.Time > 0.08)
metaepas1 <- metaepas1 %>% dplyr::filter(OS.Time < 19)
metaplot <- data.frame(epas1=metaepas1$epas1,OS=metaepas1$OS,OS.Time=metaepas1$OS.Time)

write.csv(metaplot,file = "/HDD8T/eunji/proj/hy-drug/metabric/metaplot.csv")
metaplot <- readxl::read_excel("/HDD8T/eunji/proj/hy-drug/metabric/metaplot.xlsx",
                               sheet = "Sheet1",
                               range = "A1:D256", 
                               col_names = T,
                               na="NA")
table(metaplot$epas1)

library(prodlim)

# "#800000FF","#374E55FF" 
kmfit <-prodlim(Hist(OS.Time,OS) ~ epas1, data = metaplot)
#"firebrick3","darkgray"
#"#800000FF","#374E55FF"
plot(kmfit,percent=FALSE,logrank=T,digits = 1,axes=TRUE,
     font=4,
     marktime=T,atrisk=T,xlab="Years",lwd=2,
     confint=F,confint.citype="shadow",#col=c(4,3),
     legend=T,legend.x=1.7,
     legend.y=0.7,legend.cex=0.9,
     legend.title="METABRIC\nEPAS1 Group\n",
     atrisk.labels=paste(c("M1","M2"),": "),
     atrisk.title="")
plot(kmfit,percent=FALSE,logrank=T,digits = 1,axes=TRUE)

metarna3 <- data.frame(row.names = metarna$Hugo_Symbol,metarna[,3:1906])
metarna3 <- as.matrix(metarna3)
ssmeta<- gsva(metarna3,list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ssmeta2 <- as.data.frame(ssmeta)
colnames(ssmeta2)<- gsub('.','-', colnames(ssmeta2),fixed = T)

metaepas1$number <- rownames(metaepas1)
names(metaplot)[1] <- "number"
meta <- metaepas1[,c(1:40,45)]
meta2<- merge(metaplot,meta,by="number")

#ssmeta2 <- ssmeta2[,c(meta2$SampleID)]
ssmeta2

ssmeta3 <- as.data.frame(apply(ssmeta2,1,function(x){(x-mean(x))/sd(x)}))
ssmeta3 <- data.frame(SampleID=rownames(ssmeta3),ssmeta3)
#ssmeta4 <- merge(ssmeta3,meta2,by="SampleID")
ssmeta4 <- merge(ssmeta3,metaepas1,by="SampleID")
ssmeta5 <- ssmeta4[,c(1:3,42:45)]
table(ssmeta5$epas1)
table(ssmeta5$EPAS1group)
ssmeta5$EPAS1group <- NULL
write.csv(ssmeta5,file = "/HDD8T/eunji/proj/hy-drug/ssmeta5.csv",quote = F)



ssmeta6 <- readxl::read_excel("/HDD8T/eunji/proj/hy-drug/ssmeta6.xlsx",
                             sheet = "Sheet1",
                             range = "A1:H227", 
                             col_names = T,
                             na="NA")

table(ssmeta6$EPAS1group)
table(ssmeta6$epas1)


library(ggsci)
library(ggpubr)
ggscatter(ssmeta6,x = "ENDORSE", y = "EPAS1",
          fill = "black",
          #color = "darkgrey",
          shape = 21,size = 2.2,
          add = "reg.line",cor.coef = TRUE,cor.method = "pearson",
          conf.int = TRUE,
          title="",xlab = "Endocrine resistance score", ylab = "EPAS1")

ggplot(ssmeta6, aes(x=epas1, y=ENDORSE, fill=epas1)) + 
  #geom_violin(trim=FALSE)+
  #geom_boxplot(width=0.1, fill="white")+
  geom_boxplot()+
  labs(title="",x="", y = "ENDORSE")+
  #theme_bw()+
  theme_classic()+
  scale_fill_manual(values = c("#999999", "#E69F00"))+
  stat_compare_means(method = "t.test",label.y = 1)

#cellularity score
table(meta2$Cellularity)
meta3 <- meta2[-4,]
meta3$Cellularity <- factor(meta3$Cellularity,levels=c("High","Moderate","Low"))
ggplot(meta3, aes(x=Cellularity, y=EPAS1,fill=Cellularity)) +
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_compare_means(label.y = 10.9)


#univariate cox proportional hazard regression analysis 
clin <- data[,c(1:16,21683)]
meth <- data[,c(1,17:21682)]
rownames(meth) <- meth[,1]
meth <- meth[,-1]
covariates <- c(colnames(meth))
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, fustat)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
                         res<-c(beta, HR, CI , wald.test, p.value)
                         names(res)<-c("beta", "HR", "CI 95%", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res,stringAsFactors=F)
res <- subset(res, p.value < 0.05 & abs(as.numeric(beta))< 5)
res <- res[order(res$p.value),]
