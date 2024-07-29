#TCGA LumA TME 
#brca-pam50 
bcpam<- readxl::read_excel("/HDD8T/eunji/proj/lg/brcapam50.xlsx",
                           sheet = "Sheet1",
                           range = "A1:B1070", 
                           col_names = F,
                           na="NA")
names(bcpam) <- c('sampleID','subtype')

bcpam2<- readxl::read_excel("/HDD8T/eunji/proj/lg/brcapam50.xlsx",
                            sheet = "Sheet2",
                            range = "A2:U827", 
                            col_names = T,
                            na="NA")
bcpam2 <- bcpam2[,c(1,21)]
names(bcpam2) <- c('sampleID','subtype')

bcpam3<- readxl::read_excel("/HDD8T/eunji/proj/lg/brcapam50.xlsx",
                            sheet = "Sheet3",
                            range = "A2:C730", 
                            col_names = F,
                            na="NA")
bcpam3 <- bcpam3[,c(1,3)]
names(bcpam3) <- c('sampleID','subtype')

load("/HDD8T/eunji/proj/lg/lgdata1.RData")

brca <- lgdata1[['BRCA']]$coldata
brca <- data.frame(sampleID=rownames(brca),brca)
brca <- brca %>% dplyr::filter(group=="T")
brca$sampleID <- gsub('.','-', str_sub(rownames(brca),1,12),fixed = T)
brca <- merge(brca,bcpam,by="sampleID",all.x=T)
brca <- brca[!duplicated(brca$sampleID),]
brca <- merge(brca,bcpam2,by="sampleID",all.x=T)
brca <- merge(brca,bcpam3,by="sampleID",all.x=T)

load("/HDD8T/eunji/proj/lg/brcasurv.RData")
names(survival_data)[1] <- 'sampleID'
brca <- merge(brca,survival_data,by="sampleID",all.x=T)
lumA <- brca %>% dplyr::filter(subtype.x=="BRCA.LumA")

count <- as.data.frame(log2(t(lgdata1[["BRCA"]]$tpmcount+1)))
count$sampleID <- gsub('.','-', str_sub(rownames(count),1,12),fixed = T)
lumA <- merge(lumA,count,by="sampleID")
lumA <- lumA %>% dplyr::filter(OS.Time > 0.08)
lumA$epas1 <- ifelse(lumA$EPAS1 > median(lumA$EPAS1),"C1","C2")

lumA <- lumA%>% dplyr::filter(sampleID != "TCGA-BH-A8FY")%>%
  dplyr::filter(sampleID != "TCGA-PE-A5DC")%>%
  dplyr::filter(sampleID != "TCGA-AR-A0TR")%>%
  dplyr::filter(sampleID != "TCGA-B6-A0X4")%>%
  dplyr::filter(sampleID != "TCGA-AC-A6IV")%>%
  dplyr::filter(sampleID != "TCGA-AC-A8OR")%>%
  dplyr::filter(sampleID != "TCGA-AC-A2B8")

lumA <- lumA[order(lumA$epas1),]

library(EPIC)
tpmcount <- as.data.frame(t(lgdata1[["BRCA"]]$tpmcount))
tpmcount$sampleID <- gsub('.','-', str_sub(rownames(count),1,12),fixed = T)
luma <- lumA[,c(1:11,49315)]
luma <- luma[order(luma$epas1),]
lumatpm <- merge(luma,tpmcount)
lumatpm <- as.data.frame(t(data.frame(row.names=lumatpm$sampleID,lumatpm[,12:49314])))
lumatpm <- lumatpm[-1,]
library(xCell)
xcell <- xCellAnalysis(lumatpm,rnaseq=T)
xcell <- as.data.frame(t(xcell))
xcell$sampleID <- rownames(xcell)
xcell <- merge(luma,xcell,by="sampleID")
xcell <- xcell[order(xcell$epas1),]
xcell2 <- as.data.frame(t(data.frame(row.names = xcell$sampleID,xcell[,13:79])))

ttestss <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

xcell2$pvalue = apply(xcell2, 1, ttestss, grp1 = c(1:216), grp2 = c(217:431))
xcell2$padj <- p.adjust(xcell2$pvalue,"fdr")
xcell2$mean1 <- apply(xcell2[1:216],1,function(x)(mean(x)))
xcell2$mean2 <- apply(xcell2[217:431],1,function(x)(mean(x)))
xcell2$mean_diff <- xcell2$mean1-xcell2$mean2

luma <- as.data.frame(t(data.frame(row.names = lumA$sampleID,lumA[,12:49314])))
luma$c1 <- apply(luma[,1:216],1,mean)
luma$c2 <- apply(luma[,217:431],1,mean)

luma$fc <- luma$c1/luma$c2
luma$logfc <- log2(luma$c1/luma$c2)

for(h in 1:nrow(luma)){
  x=as.numeric(luma[h, 1:216])
  y=as.numeric(luma[h, 217:431])
  luma$p[h]=t.test(x, y,alternative = "two.sided")$p.value}
luma$padj <- p.adjust(luma$p,"fdr")
resluma <- data.frame(gene_name=rownames(luma),luma[,432:437])

save(luma,file = "/HDD8T/eunji/proj/lg/luma/luma.RData")

#TCGA LumA TME
#Treg & exhausted T cell ("ENTPD1","CD4","TGFB1","CAT")
#exhausted T cell markers - PDCD1, CTLA4, HAVCR2, TIGIT, TNFRSF9, CD27, TOX, and ENTPD1
treg <- lumA[,c("sampleID","epas1","PDCD1","CTLA4","HAVCR2","TIGIT",
                "TNFRSF9","CD27","TOX","ENTPD1","CAT",
                "CD3D","CD4","IL2RA","FOXP3")]

library(reshape)
tmelt <- melt(data=treg,id.vars = c(names(treg)[c(1,2)]),
              measure.vars = c(names(treg)[3:15]))

tanno <- data.frame(variable=c(names(treg)[3:15]),
                    padj=c(5.090505e-05,0.0006234139,1.262409e-05,1.280672e-06,
                           1.093812e-05,7.626017e-06,8.784651e-10,2.198157e-24,1.504747e-10,
                           3.355822e-06,1.656437e-08,2.615121e-07,6.562364e-03))
tmelt2 <- merge(tmelt,tanno,by="variable",all.x=T)
tanno$group1 <- "C1"
tanno$group2 <- "C2"
tanno <- transform(tanno, 
                   psig = cut(padj, breaks = c(0, 0.0001, 0.001, 0.01, 0.05,1),
                              include.lowest = TRUE,right = FALSE, 
                              labels = c("****", "***", "**", "*", "ns")))

library(ggpubr)

tmelt2 %>%
  dplyr::filter(variable %in% c("HAVCR2")) %>%
  ggplot(aes(x = variable, y = value, fill = epas1)) +
  geom_boxplot(alpha=0.9) +
  #stat_summary(fun=mean, geom="point",aes(group=subtype), position=position_dodge(.7),color="black")+
  scale_y_continuous(name = NULL) +
  scale_x_discrete(name = NULL,labels=NULL,breaks=NULL) +
  ggtitle("HAVCR2") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"),
        legend.position = "") +
  scale_fill_manual(values=c("firebrick3","#3F5151"))+
  #scale_fill_brewer(wes_palette("BottleRocket2", 2)) +
  labs(fill = "group",x=NULL,y="Expression")+
  add_pvalue(tanno[3,], inherit.aes = FALSE,
             xmin = "group1",
             x="variable",
             label = "psig",
             y.position = 4.5)

library(dplyr)
library(tidyr)
library(tibble)

xcell3 <- xcell2[1:64,]
xcell3 <- xcell3[order(xcell3$padj),]
write.csv(xcell3,file = "/HDD8T/eunji/proj/lg/luma/xcell3-2.csv")

xcell3<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/xcell3-2.xlsx",
                            sheet = "Sheet1",
                            range = "A1:PW65", 
                            col_names = T,
                            na="NA")
xcell5 <- xcell3[order(xcell3$Immune.Celltype),]
xcell5 <- xcell5[1:21,]
xcell6 <- data.frame(row.names = xcell5$Celltype,xcell5[,5:435])
xcell5_2 <- xcell5[,4:435]
xcell5_3 <- xcell5_2 %>% group_by(Immune.Celltype) %>% summarise_at(vars(names(xcell5_2)[2:432]),funs(sum(.,na.rm=TRUE)))
xcell5_3 <- data.frame(row.names = xcell5_3$Immune.Celltype,xcell5_3[,2:432])

xcell6 <- as.data.frame(apply(xcell5_3, 2, function(x) {x/sum(x)}))
xcell6$C1 <- apply(xcell6[1:216],1,function(x)(mean(x)))
xcell6$C2 <- apply(xcell6[217:431],1,function(x)(mean(x)))
xcell7 <- as.data.frame(t(xcell6))
xcell7 <- xcell7[c(432:433),]
rownames(xcell7) <- c('C1','C2')

dd1 <- xcell7 %>% 
  as.data.frame() %>% 
  rownames_to_column("subtype") %>% 
  pivot_longer(cols=2:10,
               names_to= "Immune cell type",
               values_to = "Proportion")


colorblind_pal()(8)
ggplot(dd1,aes(subtype,Proportion,fill = `Immune cell type`)) +
  geom_col(colour = "black", position = "fill") +
  theme_few()+
  theme(plot.title = element_text(size = 13, family = "Arial", face = "bold"),
                    text = element_text(size = 12, family = "Arial"),
                    #axis.title = element_text(face="bold"),
                    axis.text.x=element_text(size = 10,face = "bold"))+
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#606060","#F8AFA8","#56B4E9","#E69F00","#009E73",
                               "#F0E442","#CC79A7","#0072B2","#D55E00"))


tcell <- xcell5[order(xcell5$Tcell),]
tcellres <- xcell5[,c(1:3,436:439)]
tcell2 <- tcell[1:10,]
tcell2 <- data.frame(row.names = tcell2$Tcell,tcell2[,5:435])

tcell2 <- as.data.frame(apply(tcell2, 2, function(x) {x/sum(x)}))
tcell2$C1 <- apply(tcell2[1:216],1,function(x)(mean(x)))
tcell2$C2 <- apply(tcell2[217:431],1,function(x)(mean(x)))

tcell3 <- as.data.frame(t(tcell2))
tcell3 <- tcell3[c(432:433),]

dd2 <- tcell3 %>% 
  as.data.frame() %>% 
  rownames_to_column("subtype") %>% 
  pivot_longer(cols=2:11,
               names_to= "celltype",
               values_to = "Proportion")

names(dd2)[2] <- 'T cell type'
stata_pal("s2color")(15)
ggplot(dd2,aes(subtype,Proportion,fill = `T cell type`)) +
  geom_col(colour = "black", position = "fill") +
  theme_dark()+
  theme(plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        text = element_text(size = 12, family = "Arial"),
        #axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"))+
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(brewer.pal(10,"Set3")))
library(RColorBrewer)

#Treg marker - "CD3D","CD4","IL2RA","FOXP3"

tregs <- tcell[1:10,]
tregs <- data.frame(row.names = tregs$Tcell,tregs[,5:435])
tregs <- as.data.frame(t(tregs))
tregs$subtype <- c(rep("C1",216),rep("C2",215))

ggplot(tregs, aes(x=subtype, y=Tregs,fill=subtype)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c("firebrick3","#3F5151"))+
  labs(x='', y='Percentage')+
  theme_bw()+
  theme(plot.title = element_text(size = 13, family = "Arial", face = "bold"),
        text = element_text(size = 12, family = "Arial"),
        #axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"))+ylim(0,1.5)+
  geom_text(x=1.5, y=1.5, label="P < 0.05",size=4,family = "Arial",face = "bold")


#stromal cells
library(pheatmap)
scell<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/xcell3-2.xlsx",
                            sheet = "Sheet1",
                            range = "A1:PX65", 
                            col_names = T,
                            na="NA")
scell <- scell[order(scell$stromal),]
scell <- scell[1:14,]
scellclin <- scell[,c(1,436:440)]
scell <- data.frame(row.names = scell$stromal,scell[,5:435])
colnames(scell) <- gsub('.','-',colnames(scell),fixed = T)
subtype <- data.frame(row.names = lumA$sampleID,subtype=lumA$epas1)
scell <- scell[,c(rownames(subtype))]
scell2 <- as.data.frame(t(scell))
scell2 <- sapply(scell2, function(df) (df-mean(df))/sd(df))
scell2 <- as.data.frame(scell2)
rownames(scell2)<- colnames(scell)
scell3 <- as.data.frame(t(scell2))
scell3[scell3 > 2] =2
scell3[scell3 < -2] =-2
subtype$subtype <- as.factor(subtype$subtype)

mat_colors <- list(group = c("firebrick3","#3F5151"))
names(mat_colors$group) <- c("C1","C2")
ann_colors = list(
  subtype = c(C1 = "firebrick3", C2 = "#3F5151"))

paletteLength <- 50 #FCFFA4FF
myColor <- colorRampPalette(c("#374E55FF",'black',"#FFCD00FF"))(paletteLength)
myColor <- colorRampPalette(c("navy",'white',"darkred"))(paletteLength)
myBreaks <- c(seq(min(scell3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scell3)/paletteLength, max(scell3), length.out=floor(paletteLength/2)))

pheatmap(scell3,scale="none",color = myColor,breaks = myBreaks,
         fontsize_col = 12, fontsize_row=12,
         show_rownames = T,show_colnames = F,
         cluster_rows =T, cluster_cols = F,
         annotation_col = subtype,annotation_colors = ann_colors,
         #cutree_rows = 2,cutree_cols = 3,
         legend = T)

#estimate 
library(estimate)
load("/HDD8T/eunji/proj/multi/TCGAtcount.RData")
rawluma <- data.frame(sampleID=tcount$Barcode,tcount[,9:49311])
rawluma <- merge(lumA[,c(1,49315)],rawluma,by="sampleID")
raw2 <- as.data.frame(t(data.frame(row.names = rawluma$sampleID,rawluma[,3:49305])))
raw2 <- data.frame(GeneSymbol=rownames(raw2),raw2)
colnames(raw2) <- gsub('.','-',colnames(raw2),fixed = T)
write.table(raw2,file = "/HDD8T/eunji/proj/lg/luma/input.txt",sep = "\t",col.names = T,row.names = F,quote = F)
filterCommonGenes(input.f = "/HDD8T/eunji/proj/lg/luma/input.txt",
                  output.f = "/HDD8T/eunji/proj/lg/luma/output.gct",
                  id='GeneSymbol')

estimateScore(input.ds = "/HDD8T/eunji/proj/lg/luma/output.gct",
              output.ds = "/HDD8T/eunji/proj/lg/luma/estimate_score.gct",
              platform = "affymetrix")

scores <- read.table("/HDD8T/eunji/proj/lg/luma/estimate_score.gct",skip = 2,header = T)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])
scores <- as.data.frame(scores)
scores <- data.frame(sampleID=rownames(scores),scores)
scores$sampleID <- gsub('.','-',scores$sampleID,fixed = T)
scores <- merge(lumA[,c(1,49315)],scores,by="sampleID")

ggplot(scores,aes(x = epas1, y = ImmuneScore, fill = epas1)) +
  #geom_violin()+
  geom_boxplot(alpha=0.9) +
  stat_compare_means(method = "t.test",label.y = 3000)+
  #stat_summary(fun=mean, geom="point",aes(group=subtype), position=position_dodge(.7),color="black")+
  scale_y_continuous(name = NULL) +
  scale_x_discrete(name = NULL,labels=NULL,breaks=NULL) +
  ggtitle("ImmuneScore") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold"),
        legend.position = "") +
  scale_fill_manual(values=c("firebrick3","#3F5151"))+
  #scale_fill_brewer(wes_palette("BottleRocket2", 2)) +
  labs(fill = "group",x=NULL,y="")


#mutation 
library(maftools)
library(TCGAbiolinks)
#https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/mutation.html

library(maftools)

mut <- read.maf(maf = "/HDD8T/eunji/proj/lg/mc3.v0.2.8.PUBLIC.maf")
summary <- getSampleSummary(mut)
summary$Tumor_Sample_Barcode <- str_sub(summary$Tumor_Sample_Barcode,1,12)
summary = summary[-which(duplicated(summary$Tumor_Sample_Barcode)),]
mut@data$Tumor_Sample_Barcode <- str_sub(mut@data$Tumor_Sample_Barcode,1,12)
mut@variants.per.sample$Tumor_Sample_Barcode <- str_sub(mut@variants.per.sample$Tumor_Sample_Barcode,1,12)
mut@variant.type.summary$Tumor_Sample_Barcode <- str_sub(mut@variant.type.summary$Tumor_Sample_Barcode,1,12)
mut@variant.classification.summary$Tumor_Sample_Barcode <- str_sub(mut@variant.classification.summary$Tumor_Sample_Barcode,1,12)
mut@clinical.data$Tumor_Sample_Barcode <- str_sub(mut@clinical.data$Tumor_Sample_Barcode,1,12)
mut@maf.silent$Tumor_Sample_Barcode <- str_sub(mut@maf.silent$Tumor_Sample_Barcode,1,12)

load("/HDD8T/eunji/proj/lg/luma/lumA.RData")

subtype2 <- data.frame(Tumor_Sample_Barcode=lumA$sampleID,subtype=lumA$epas1)
c1 <- subtype2 %>% dplyr::filter(subtype == "C1")
c2 <- subtype2 %>% dplyr::filter(subtype == "C2")
mut2 <- read.maf(maf = "/HDD8T/eunji/proj/lg/TCGA.BRCA.varscan3.maf")

mutc1 <- subsetMaf(maf=mut2,tsb=c(c1$Tumor_Sample_Barcode))
mutc2 <- subsetMaf(maf=mut2,tsb=c(c2$Tumor_Sample_Barcode))
mutc1@clinical.data$subtype="C1"
mutc2@clinical.data$subtype="C2"
mut3 <- merge_mafs(maf=c(mutc1,mutc2),verbose=T)
mutsig <- prepareMutSig(mut3)

oncoplot(maf = mutc1,top = 20)
oncoplot(maf = mutc2,top = 20)

mut2 <- subsetMaf(maf=mut,tsb=c(subtype2$Tumor_Sample_Barcode)) #N=390
plotmafSummary(maf = mut2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
summary <- getSampleSummary(mut3)
genesum1 <- getGeneSummary(mutc1)
genesum2 <- getGeneSummary(mutc2)

coc1<- somaticInteractions(maf = mutc1, top = 20, pvalue = c(0.05, 0.1),showCounts = F,fontSize = 0.6)
coc2 <- somaticInteractions(maf = mutc2, top = 20, pvalue = c(0.05, 0.1),showCounts = F,fontSize = 0.6)
C1C2 <- mafCompare(m1 = mutc1, m2 = mutc2, m1Name = 'C1', m2Name = 'C2', minMut = 5)
print(C1C2)

coOncoplot(m1 = mutc1, m2 = mutc2, m1Name = 'C1', m2Name = 'C2', removeNonMutated = TRUE)

#Changing colors (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = RColorBrewer::brewer.pal(n = 2, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
#We will plot same top ten mutated genes with FAB classification as annotation and using above defined colors.
oncoplot(maf = mut2, top = 10, annotation = subtype2, removeNonMutated = TRUE)
oncoplot(maf=mut2)

#Changing colors for variant classifications (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

#Color coding for FAB classification; try getAnnotations(x = laml) to see available annotations.
fabcolors = c("firebrick3","darkgray")
names(fabcolors) = c("C1","C2")
fabcolors = list(subtype = fabcolors)

oncoplot(maf = mut3, colors = col,clinicalFeatures = 'subtype', sortByAnnotation = TRUE, annotationColor = fabcolors)
summary <- getSampleSummary(mut3)
summary <- summary %>%
  dplyr::mutate(TMB=total/38) %>%
  dplyr::mutate(TMB_group=ifelse(TMB>mean(TMB),'high','low'))
summary$logTMB <- log(summary$TMB)
summary <- merge(summary,subtype2,by="Tumor_Sample_Barcode")
summary <- summary %>% dplyr::filter(TMB > 0)
summary <- summary %>% dplyr::filter(Tumor_Sample_Barcode != "TCGA-PE-A5DE")
summary <- summary %>% dplyr::filter(Tumor_Sample_Barcode != "TCGA-E9-A1RI")


ggplot(summary, aes(x=subtype, y=logTMB, fill=subtype)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="",x="", y = "logTMB")+
  theme_classic()+
  scale_fill_manual(values = c("firebrick3","darkgray"))+
  stat_compare_means(method = "t.test",label.y = 3.4)
  
