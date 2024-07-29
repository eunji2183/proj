#GSE67916 analysis result
library(GEOquery)
library(limma)
library(umap)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(ggthemes)

# load series and platform data from GEO

gset <- getGEO("GSE67916", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "X00X0000001111111X"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("resistant","sensitive"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=54675)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

ex2 <- merge(tT,ex,by="row.names")

#metabolites-protein 
mp<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/meta-network.xlsx",
                        sheet = "Sheet1",
                        range = "A1:F31207", 
                        col_names = T,
                        na="NA")
mptable<- as.data.frame(table(mp$a))
mp <- mp %>% dplyr::filter (!a %in% c("Water","Hydrogen Ion","Adenosine triphosphate",
                                      "ADP","Phosphate","Pyrophosphate","Adenosine monophosphate",
                                      "NAD","NADP","NADPH","Oxygen","NADH","Coenzyme A","ubiquitin",
                                      "Protein lysine","O-phosphoprotamine","protamine"))

mprot<- as.data.frame(table(mp$b))

library(clusterProfiler)
prot <- mprot$Var1
prot <- bitr(prot,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

EGG <- enrichKEGG(gene = prot$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
proteggres<- EGG@result
write.csv(proteggres,file = "/HDD8T/eunji/proj/lg/luma/proteggres.csv",quote = F)

mp2<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/protegg.xlsx",
                         sheet = "Sheet1",
                         range = "A1:J140", 
                         col_names = T,
                         na="NA")
mp2 <- data.frame(ENTREZID=unlist(strsplit(mp2$geneID,split = "/")))
mp2 <- bitr(mp2$ENTREZID,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")

ex2 <- ex2[!duplicated(ex2$Gene.symbol),]
ex2 <- data.frame(SYMBOL=ex2$Gene.symbol,ex2[,10:24])

save(ex2,file = "/HDD8T/eunji/proj/lg/luma/ex2.RData")

ex2 <- merge(mp2,ex2,by="SYMBOL")
ex2 <- data.frame(row.names = ex2$SYMBOL,ex2[,-c(1:2)])
ex2 <- as.data.frame(t(ex2))

#PCA
library(factoextra)
res.pca <- prcomp(ex2, scale = TRUE)


res.pca <- PCA(ex2, graph = FALSE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 35),xlim=c(1,5))
f <- fviz_pca_ind(res.pca, pointsize = "cos2", 
                  pointshape = 21, fill.ind = "#E7B800",
                  repel = TRUE) # Avoid text overlapping (slow if many points)
pcass2 <- f[["data"]]
names(pcass2)[1] <- 'sampleID'
pcass2$group <- c(rep("TamR",8),rep("TamS",7))


colorblind_pal()(8)
wes_palettes$BottleRocket2
pal_jama("default")(7)
pal_nejm("default")(7)

ggplot(pcass2,aes(x,y ,color=group,fill=group)) + 
  geom_point(size=3) + 
  scale_color_manual(values = c("#E18727FF","#6A6599FF" )) +
  #scale_color_manual(values = c("#BC3C29FF","#3C5488FF" )) + 
  stat_ellipse(aes(color = group, fill = group),
               alpha = 0.15, geom = "polygon",level = 0.93,type="t") +
  scale_fill_manual(values = c( "#E18727FF","#6A6599FF" )) +theme_classic()+
  labs(title="",
       #subtitle = "[Hallmark gene sets]",
       #caption = "Source: TCGA & Msigdb",
       x="PC1(30.8%)", y = "PC2(16.9%)")+
  ggthemes::theme_few()+
  theme(#legend.position=c(0.87,0.25),legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold",hjust = 0.5,size = 13),
    plot.subtitle = element_text(face = "bold",hjust = 0.5,size = 10),
    plot.caption = element_text(face = "Arial",colour = "darkgray"),
    legend.key = element_rect(color="black"))

#volcano plot 
tT <- merge(ex,tT,by="row.names")
tT <- tT[!duplicated(tT$Gene.symbol),]
tT$DEG <- as.factor(ifelse(tT$adj.P.Val < 0.05 & abs(tT$logFC) > 0,
                           ifelse(tT$logFC > 0 , 'UP','DOWN'),'NOT'))

tT$sign <- ifelse(tT$adj.P.Val < 0.05 & abs(tT$logFC) > 2.5,tT$Gene.symbol,NA)
tT <- tT[order(tT$adj.P.Val),]

tT$Gene.symbol[str_which(tT$Gene.symbol, "C8orf44-SGK3///SGK3")] = 'SGK3'
tT$Gene.symbol[str_which(tT$Gene.symbol, "LOC100505984///ITGB6")] = 'ITGB6'
tT$Gene.symbol[str_which(tT$Gene.symbol, "LRRD1///CYP51A1")] = 'CYP51A1'
tT$Gene.symbol[str_which(tT$Gene.symbol, "MIR1248///SNORA81///SNORA4///SNORD2///SNORA63///EIF4A2")] = 'EIF4A2'
tT$Gene.symbol[str_which(tT$Gene.symbol, "MIR3620///ARF1")] = 'ARF1'
tT$Gene.symbol[str_which(tT$Gene.symbol, "MIR612///NEAT1")] = 'NEAT1'
tT$Gene.symbol[str_which(tT$Gene.symbol, "NPIPA5///NPIPB5///LOC613037///NPIPB4///NPIPB3")] = 'NPIPA5'
tT$Gene.symbol[str_which(tT$Gene.symbol, "RRN3P1///RRN3P2///RRN3")] = 'RRN3'
tT$Gene.symbol[str_which(tT$Gene.symbol, "CHURC1-FNTB///FNTB")] = 'FNTB'
tT$Gene.symbol[str_which(tT$Gene.symbol, "LOC101928927///SNORA9///SNHG15")] = 'SNHG15'
tT$Gene.symbol[str_which(tT$Gene.symbol, "SNORD50B///SNORD50A")] = 'SNORD50A'
tT$Gene.symbol[str_which(tT$Gene.symbol, "BCL2L2-PABPN1///PABPN1")] = 'PABPN1'
tT$Gene.symbol[str_which(tT$Gene.symbol, "HIST2H2AA4///HIST2H2AA3")] = 'HIST2H2AA3'
tT$Gene.symbol[str_which(tT$Gene.symbol, "ATP6V1G2-DDX39B///SNORD84///DDX39B")] = 'DDX39B'
tT$Gene.symbol[str_which(tT$Gene.symbol, "BCL2L2-PABPN1///PABPN1///BCL2L2")] = 'BCL2L2'

tT$Gene.symbol[str_which(tT$Gene.symbol, "MIR4784///MZT2A///MZT2B")] = 'MZT2A'
tT$Gene.symbol[str_which(tT$Gene.symbol, "NEDD8-MDP1///MDP1")] = 'MDP1'
tT$Gene.symbol[str_which(tT$Gene.symbol, "PRR5-ARHGAP8///ARHGAP8")] = 'ARHGAP8'
tT$Gene.symbol[str_which(tT$Gene.symbol, "LOC101928615///FNDC3B")] = 'FNDC3B'
tT$Gene.symbol[str_which(tT$Gene.symbol, "LOC101930578///TPTE2P2")] = 'TPTE2P2'
tT$Gene.symbol[str_which(tT$Gene.symbol, "MIR6734///ELOVL1")] = 'ELOVL1'

tT <- tT %>% dplyr::filter(ID != "AFFX-M27830_5_at")

save(tT,file = "/HDD8T/eunji/proj/lg/luma/tT.RData")
tT <- tT %>% dplyr::filter(adj.P.Val < 0.1)

write.csv(tT,file = "/HDD8T/eunji/proj/lg/luma/tTres.csv",quote = F)

tT$sign2 <- ifelse(tT$Gene.symbol %in% c("CD36","FABP6","EPAS1","HIF1A",
                                         "FBP1","ABCC3","ABCD1","ABCF2",
                                         "ACSL1","ACADL","ADAD10","CAT",
                                         "SLC1A5","SLC38A1","GCLC"),tT$Gene.symbol,NA)

tT$sign2 <- ifelse(tT$Gene.symbol %in% c("EPAS1"),tT$Gene.symbol,NA)

tT <- tT %>% dplyr::filter(adj.P.Val < 0.1)


library(ggrepel)
ggplot(data = tT, aes(x=logFC,y=-log10(adj.P.Val),col = DEG)) +
  geom_point(alpha=0.3)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  #geom_hline(yintercept=2 ,linetype=4) +
  #geom_vline(xintercept=0 ,linetype=4 ) +
  scale_color_manual(name = "", values = c("#E18727FF","#6A6599FF", "gray"), limits = c("UP", "DOWN", "NOT")) +
  geom_label_repel(aes(label=sign2), face="Arial",fontface="bold", color="grey50", box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"), segment.colour = "grey50")



#pheatmap
gene <- c(#"EPAS1",
          "EPAS1","HIF1A",
          "ACSL1","ACADL","ELOVL1",
          "SLC1A5","SLC38A1","GCLC","GCLM",
          "TGFB1","TIMP1",
          "ABCF2","ABCC3","ABCD1","ABCC1","ABCA5")

annotation_row = data.frame(row.names = gene,Gene.symbol=gene,Number=1:17,
                            GeneGroup=factor(rep(c("Hypoxia","Lipid met","Glutamine Met","EMT","Drug resistant"),
                                                 c(2,4,4,2,5))))

annorow <- data.frame(row.names = gene,
                      GeneGroup=factor(rep(c("Hypoxia","Lipid met","Glutamine Met","EMT","Drug resistant"),
                                           c(2,3,4,2,5))))

ex3 <- data.frame(ID=rownames(ex),ex)
ex3 <- merge(tT,ex3,by="ID")
heat <- ex3 %>% dplyr::filter(Gene.symbol %in% c(gene))
heat <- merge(annotation_row,heat,by="Gene.symbol")
heat <- heat[order(heat$Number,decreasing = F),]
heat2 <- data.frame(row.names = heat$Gene.symbol,heat[,6:20]) 
annotation_col=data.frame(row.names = colnames(heat2),
                          Group=factor(rep(c("TamR","TamS"),
                                           c(8,7))))

write.csv(heat2,file = "/HDD8T/eunji/proj/lg/luma/heat6.csv")
heat2 <- read.csv(file = "/HDD8T/eunji/proj/lg/luma/heat6.csv",row.names = 1)
heat2 <- heat2[-3,]

library(pheatmap)

colorblind_pal()(8)
pal_lancet("lanonc")(9)
pal_npg("nrc")(10)
pal_uchicago("dark")(9)
stata_pal("s1color")(15)
ann_colors = list(
  Group = c(TamR = "#E18727FF", TamS = "#6A6599FF"),
  GeneGroup = c(Hypoxia = "#F0E442", `Lipid met` = "#009E73", `Glutamine Met` = "#0072B2",
                EMT="#D55E00",`Drug resistant`="#925E9FFF"))


pheatmap(heat2,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(13), 
         fontsize_col = 10, fontsize_row=12,
         cellwidth = 10, cellheight = 15,show_rownames = T,show_colnames = T,legend = T,
         border="black",cluster_cols = F,cluster_rows = F,
         annotation_row = annorow,annotation_legend = T,
         annotation_col = annotation_col,annotation_colors = ann_colors)

met<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/figure4-met.xlsx",
                         sheet = "Sheet1",
                         range = "A1:C13", 
                         col_names = T,
                         na="NA")
met$Log10p <- -log10(met$padj) 
met <- met[order(met$Log10p,decreasing = T),]
met <- met[order(met$NES,decreasing = T),]
met$pathway <- factor(met$pathway,levels = c(met[order(met$NES,decreasing = T),]$pathway))

ggplot(met, aes(x = NES, y = pathway, fill = Log10p)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip()+
  scale_fill_gradient(low = "white", high = "#C16622FF")+
  xlab("")+
  ylab("")+ggtitle("")+
  theme_few()+
  theme(axis.text.x = element_text(family = "Arial",size = 9,color = "black",angle = 90,vjust = 0.5, hjust=1,face = "bold"),                      # adjusting the position
        #axis.title.x = element_text(family  = "Arial",size = 12,color = "black"),                   # face the x axit title/label
        axis.text.y = element_text(family  = "Arial",size = 10,color = "black",face = "bold"),                   # face the y axis title/label
        axis.title.y = element_text(family = "Arial",hjust = 0.1,face = "bold",size = 15))


