## GSE67916 analysis result 

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
BiocManager::install("GEOquery")
BiocManager::install("umap")
install.packages("weasanderson")
BiocManager::install("ggthemes")
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
ex3 <- ex2 %>% dplyr::filter(P.Value < 0.1 &logFC > 0)
ex4 <- ex2 %>% dplyr::filter(P.Value < 0.1 &logFC < 0)
ex5 <- rbind(ex3,ex4)

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
    plot.caption = element_text(face = "italic",colour = "darkgray"),
    legend.key = element_rect(color="black"))

#volcano plot 
tT <- merge(ex,tT,by="row.names")
tT <- ex5
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

tT <- tT %>% dplyr::filter(ID != "AFFX-M27830_5_at")

save(tT,file = "/HDD8T/eunji/proj/lg/luma/tT.RData")
tT <- tT %>% dplyr::filter(adj.P.Val < 0.1)

write.csv(tT,file = "/HDD8T/eunji/proj/lg/luma/tTres.csv",quote = F)

tT$sign2 <- ifelse(tT$Gene.symbol %in% c("CD36","EPAS1","HIF1A"),tT$Gene.symbol,NA)



library(ggrepel)
ggplot(data = tT, aes(x=logFC,y=-log10(adj.P.Val),col = DEG)) +
  geom_point(alpha=0.3)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept=2 ,linetype=4) +
  geom_vline(xintercept=0 ,linetype=4 ) +
  scale_color_manual(name = "", values = c("#E18727FF","#6A6599FF", "gray"), limits = c("UP", "DOWN", "NOT")) +
  geom_label_repel(aes(label=sign2), fontface="bold", color="grey50", box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"), segment.colour = "grey50")


gene <- c("EPAS1","VEGFA",
          "CD36","FABP6","ACADL","ACAD10","CPT2","ACAT2",
          "SLC1A5","SLC38A1","GCLC","GCLM","GSS",
          "TGFB1","TIMP1","TIMP2","CDH2","TEAD1",
          "ABCC3","ABCE1","ABCF2","ABCD1")
symbol <- c("200878_at","211527_x_at",
            "206488_s_at","210445_at","206068_s_at","237482_s_at","204263_s_at","209608_s_at",
            "208916_at","224580_at","1555330_at","236140_at","211630_s_at",
            "203085_s_at","201666_at","231579_s_at","203440_at","214600_at",
            "208161_s_at","201872_s_at","207622_s_at","205142_x_at")

annotation_row = data.frame(row.names = gene,Gene.symbol=gene,Number=1:22,
                            GeneGroup=factor(rep(c("Hypoxia","Lipid met","Glutamine Met","EMT","Drug resistant"),
                                                 c(2,6,5,5,4))))

annorow <- data.frame(row.names = gene,
                      GeneGroup=factor(rep(c("Hypoxia","Lipid met","Glutamine Met","EMT","Drug resistant"),
                                           c(2,6,5,5,4))))

#ex3 <- data.frame(ID=rownames(ex),ex)
#ex3 <- merge(tT,ex3,by="ID")
#heat <- ex3 %>% dplyr::filter(Gene.symbol %in% c(gene))
heat <- ex5 %>% dplyr::filter(ID %in% c(symbol))
heat <- merge(annotation_row,heat,by="Gene.symbol")
heat <- heat[order(heat$Number,decreasing = F),]
heat2 <- data.frame(row.names = heat$Gene.symbol,heat[,12:26]) 
rownames(heat2)[c(4,7)] <- c("FABP4","CPT1A")
gene <- c("EPAS1","VEGFA",
          "CD36","FABP4","ACADL","ACAD10","CPT1A","ACAT2",
          "SLC1A5","SLC38A1","GCLC","GCLM","GSS",
          "TGFB1","TIMP1","TIMP2","CDH2","TEAD1",
          "ABCC3","ABCE1","ABCF2","ABCD1")
symbol <- c("200878_at","211527_x_at",
            "206488_s_at","210445_at","206068_s_at","237482_s_at","204263_s_at","209608_s_at",
            "208916_at","224580_at","1555330_at","236140_at","211630_s_at",
            "203085_s_at","201666_at","231579_s_at","203440_at","214600_at",
            "208161_s_at","201872_s_at","207622_s_at","205142_x_at")

annotation_row = data.frame(row.names = gene,Gene.symbol=gene,Number=1:22,
                            GeneGroup=factor(rep(c("Hypoxia","Lipid met","Glutamine Met","EMT","Drug resistant"),
                                                 c(2,6,5,5,4))))

annorow <- data.frame(row.names = gene,
                      GeneGroup=factor(rep(c("Hypoxia","Lipid met","Glutamine Met","EMT","Drug resistant"),
                                           c(2,6,5,5,4))))




annotation_col=data.frame(row.names = colnames(heat2),
                          Group=factor(rep(c("TamR","TamS"),
                                           c(8,7))))


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
         cellwidth = 10, cellheight = 12,show_rownames = T,show_colnames = T,legend = T,
         border="black",cluster_cols = F,cluster_rows = F,
         annotation_row = annorow,annotation_legend = T,
         annotation_col = annotation_col,annotation_colors = ann_colors)

#GSEA
gene <- c("EPAS1","VEGFA",
          "CD36","FABP6","ACADL","ACAD10","CPT2","ACAT2",
          "SLC1A5","SLC38A1","GCLC","GCLM","GSS",
          "TGFB1","TIMP1","TIMP2","CDH2","TEAD1",
          "ABCC3","ABCE1","ABCF2","ABCD1")
symbol <- c("200878_at","211527_x_at",
            "206488_s_at","210445_at","206068_s_at","237482_s_at","204263_s_at","209608_s_at",
            "208916_at","224580_at","1555330_at","236140_at","211630_s_at",
            "203085_s_at","201666_at","231579_s_at","203440_at","214600_at",
            "208161_s_at","201872_s_at","207622_s_at","205142_x_at")

tT3 <- tT %>% dplyr::filter(Row.names %in% c(symbol))
tT4 <- tT %>% dplyr::filter(!Gene.symbol %in% c(gene))
tT4 <- tT4[!duplicated(tT4$Gene.symbol),]
tT5 <- rbind(tT3,tT4)

library(org.Hs.eg.db)
genelist<-tT5$Gene.symbol
genelist<- genelist %>%  bitr(fromType = "SYMBOL",
                              toType = c("ENSEMBL", "ENTREZID"),
                              OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"Gene.symbol"
DEG<-genelist %>%
  inner_join(tT5,by='Gene.symbol') %>% 
  ## select配合everything排序把，改变变量顺序
  dplyr::select(ENTREZID,logFC,everything())
DEG <- DEG[!duplicated(DEG$ENTREZID),]
geneList<-DEG$logFC
names(geneList)<-as.character(DEG$ENTREZID)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
library(enrichplot)
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")%>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_df,pvalueCutoff = 0.05,verbose = FALSE)
em2 <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2c2,pvalueCutoff = 0.05,verbose = FALSE)

gseaplot2(em, geneSetID = 169, title = em$ID[169],color="#E18727FF",pvalue_table=F,ES_geom = "dot")
gseaplot2(em, geneSetID = 1, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

#ssGSEA
ex4 <- data.frame(row.names = tT$Gene.symbol,tT[,2:16])
load("~/R/tmp/py/gene_set.Rdata")
ex4 <- as.matrix(ex4)
library(GSA)
geneSet <-GSA.read.gmt("/HDD8T/eunji/proj/lg/luma/msigdb.v7.5.1.symbols.gmt")
names(geneSet[["genesets"]]) <- geneSet$geneset.names
geneSet <- geneSet$genesets
ssgsea<- gsva(ex4,geneSet,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ss <- as.data.frame(ssgsea)
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

ss$pvalue = apply(ss, 1, ttestss, grp1 = c(1:8), grp2 = c(9:15))
ss$padj <- p.adjust(ss$pvalue,"fdr")
ss$mean1 <- apply(ss[1:8],1,function(x)(mean(x)))
ss$mean2 <- apply(ss[9:15],1,function(x)(mean(x)))
ss$mean_diff <- ss$mean1-ss$mean2
ss <- data.frame(geneset=rownames(ss),ss)
ss <- ss[order(ss$padj),]
sshall <- ss[grep("HALLMARK", ss$geneset, ignore.case = TRUE),]
sshall <- sshall %>% dplyr::filter(padj < 0.05)
write.csv(sshall,file = "/HDD8T/eunji/proj/lg/luma/sshall.csv")

sshall <- read.csv("/HDD8T/eunji/proj/lg/luma/sshall.csv",sep = ",")
names(sshall)[1] <- 'geneSet'
sshall$logFDR <- -log10(sshall$padj) 
table(sshall$Hallmarks)

pal_nejm("default")(7)
ggplot(data=sshall,aes(x = mean_diff, y = logFDR,label=geneSet),col=Hallmarks) +
  geom_point(
    aes(color = Hallmarks, fill = Hallmarks),
    size = 5, alpha = 0.7, 
    shape = 21)+
  theme_clean() +
  scale_color_manual(values = c("Metabolism"="#20854EFF","Immune"="#7876B1FF",
                                "Cellular Stress"="#E69F00",
                                "Oncogenic"="#BC3C29FF","Stromal"="#0072B5FF",
                                "Other"="gray"))+
  scale_fill_manual(values = c("Metabolism"="#20854EFF","Immune"="#7876B1FF",
                                "Cellular Stress"="#E69F00",
                                "Oncogenic"="#BC3C29FF","Stromal"="#0072B5FF",
                                "Other"="gray"))+
  geom_vline(xintercept=0 ,linetype=1 ,color="gray") +
  geom_text_repel(aes(label = sshall$geneSet),size = 3,color="#404040")+
  labs(#title="Factors map",
    #caption = "Source: TCGA & Msigdb",
    x="Mean_diff(TamR-TamS)", y = "log10(FDR)")

#metabolites-based NES 
genelist<-tT$Gene.symbol
genelist<- genelist %>%  bitr(fromType = "SYMBOL",
                              toType = c("ENSEMBL", "ENTREZID"),
                              OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"Gene.symbol"
tT$FC <- 2^(tT$logFC)
DEG<-genelist %>%
  inner_join(tT,by='Gene.symbol') %>% 
  ## select配合everything排序把，改变变量顺序
  dplyr::select(ENTREZID,logFC,everything())
DEG <- DEG[!duplicated(DEG$ENTREZID),]
geneList<-DEG$logFC
names(geneList)<-as.character(DEG$ENTREZID)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
library(enrichplot)
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")%>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_df,pvalueCutoff = 0.05,verbose = FALSE)
emres <- em@result
gseaplot2(em, geneSetID = 169, title = em$ID[169],color="#E18727FF",pvalue_table=F,ES_geom = "dot")

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

#clusterprofile -GO & kegg
library(clusterProfiler)
ex6 <- ex2 %>% dplyr::filter(adj.P.Val < 0.05 &logFC > 0.2)
ex6 <- ex6[!duplicated(ex6$Gene.symbol),]
gene <- upmerge$gene_name
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
load("/HDD8T/eunji/proj/2D3D/ID.RData")
names(ID)[2] <- 'SYMBOL' 
gene2 <- merge(ID,gene,by="SYMBOL")
gene2 <- gene2 %>% dplyr::filter(gene_biotype == "protein_coding")
de <- gene2$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
upgoresult <- go@result
write.csv(upgoresult,file = "/HDD8T/eunji/proj/lg/luma/upgoresult.csv",quote = F)
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")

EGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.5)
upeggres<- EGG@result
write.csv(upeggres,file = "/HDD8T/eunji/proj/lg/luma/upeggres.csv",quote = F)


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

go<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/Book3.xlsx",
                         sheet = "Sheet1",
                         range = "A1:F20", 
                         col_names = T,
                         na="NA")

go$padj <- -log10(go$p.adjust)
go$pathway <- as.factor(go$Pathway)
theTable2 <- within(go, 
                    Description <- factor(Description,levels=go[order(pathway,padj),]$Description))

ggplot(theTable2, aes(x=padj, y=Description,fill=pathway)) +
  geom_bar(stat="identity",color="black")+
  #scale_fill_brewer(palette="Dark2")+
  scale_fill_manual(values = c("Metabolic"="#1B9E77","Oncogenic"="#D95F02","Transport"="#7570B3"))+
  xlab("-log10(p.adjust)")+
  ylab("")+ggtitle("")+
  #geom_vline(xintercept = 1.30,linetype='dashed')+
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold",vjust = 0.5,hjust = 1,family = "arial"),                      # adjusting the position
        axis.title.x = element_text(face = "bold",family = "arial"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 10,family = "arial"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 15))


# GSE67916, GSE128458,  GSE106681, GSE118713, GSE86538, GSE26459

##안되는 애들 - GSE104985, GSE111151, GSE55343,GSE129544,GSE75369

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
library(DESeq2)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE128458", "file=GSE128458_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- "0011XXXXXXXXXXXXXXXXXXXX11XXXXXXXXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("control","TamR"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
normalized_counts <- counts(ds)
ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj),] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))
tT$log2FoldChange <- -tT$log2FoldChange

dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts
dat <- as.data.frame(dat)
dat <- data.frame(GeneID=rownames(dat),dat)
tT3 <- merge(tT,dat,by="GeneID")

count1 <- as.data.frame(t(data.frame(row.names = tT3$Symbol,tT3[,10:15])))

plot2 <- data.frame(sampleID=rownames(count1),EPAS1=count1$EPAS1,
                    response=c("Con","Con","TamR","TamR","TamR","TamR"))

ggplot(data=plot2,aes(x = response, y = EPAS1, fill = response)) +
  geom_violin(trim=F)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold",colour = "black"),
        legend.position = "")+theme_bw()+
  labs(title="",x="", y = "Log10(EPAS1)+1")+
  ylim(0.5,4.5)
  stat_compare_means(method = "t.test", label.y = 4.4)





# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE26459", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "XXXXXX000XXXXXXXXX111XXX"
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
exprs(gset) <- log10(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("control","TamR"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=15000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

count2 <- merge(tT,ex,by="row.names")
count2<- count2[!duplicated(count2$Gene.symbol),]

count2 <- as.data.frame(t(data.frame(row.names = count2$Gene.symbol,count2[,10:15])))

plot2 <- data.frame(sampleID=rownames(count2),EPAS1=count2$EPAS1,
                    response=c("Con","Con","Con","TamR","TamR","TamR"))
plot2$EPAS1 <- plot2$EPAS1-3

ggplot(data=plot2,aes(x = response, y = EPAS1, fill = response)) +
  geom_violin(trim=F)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold",colour = "black"),
        legend.position = "")+theme_bw()+
  labs(title="",x="", y = "Log10(EPAS1)+1")+
  ylim(0.5,5)
stat_compare_means(method = "t.test", label.y = 4.4)

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
library(DESeq2)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE86538", "file=GSE86538_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- "0011XXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("control","TamR"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj)[1:17351],] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))

tT$log2FoldChange <- -tT$log2FoldChange

dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts
dat <- as.data.frame(dat)
dat <- data.frame(GeneID=rownames(dat),dat)
tT3 <- merge(tT,dat,by="GeneID")

count3 <- as.data.frame(t(data.frame(row.names = tT3$Symbol,tT3[,10:13])))

plot3 <- data.frame(sampleID=rownames(count3),EPAS1=count3$EPAS1,
                    response=c("Con","Con","TamR","TamR"))

ggplot(data=plot3,aes(x = response, y = EPAS1, fill = response)) +
  geom_violin(trim=F)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold",colour = "black"),
        legend.position = "")+theme_bw()+
  labs(title="",x="", y = "Log10(EPAS1)+1")+
  ylim(0.5,4.5)
stat_compare_means(method = "t.test", label.y = 4.4)


# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
library(DESeq2)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE118713", "file=GSE118713_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- "000111XXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("control","TamR"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj)[1:17790],] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))

dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts
dat <- as.data.frame(dat)
dat <- data.frame(GeneID=rownames(dat),dat)
tT3 <- merge(tT,dat,by="GeneID")

count4 <- as.data.frame(t(data.frame(row.names = tT3$Symbol,tT3[,10:15])))

plot3 <- data.frame(sampleID=rownames(count4),EPAS1=count4$EPAS1,
                    response=c("Con","Con","Con","TamR","TamR","TamR"))

ggplot(data=plot3,aes(x = response, y = EPAS1, fill = response)) +
  geom_violin(trim=F)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold",colour = "black"),
        legend.position = "")+theme_bw()+
  labs(title="",x="", y = "Log10(EPAS1)+1")+
  ylim(0.5,5)
stat_compare_means(method = "t.test", label.y = 4.4)


# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
library(DESeq2)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE106681", "file=GSE106681_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- "00001111XXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
tbl <- tbl[ ,sel]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("control","TamR"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

ds <- DESeq(ds, test="Wald", sfType="poscount")

# extract results for top genes table
r <- results (ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod ="fdr")

tT <- r[order(r$padj)[1:18200],] 
tT <- merge(as.data.frame(tT), annot, by=0, sort=F)

tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))

dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts
dat <- as.data.frame(dat)
dat <- data.frame(GeneID=rownames(dat),dat)
tT4 <- merge(tT,dat,by="GeneID")

count4 <- as.data.frame(t(data.frame(row.names = tT4$Symbol,tT4[,10:17])))

plot4 <- data.frame(sampleID=rownames(count4),EPAS1=count4$EPAS1,
                    response=c("Con","Con","Con","Con","TamR","TamR","TamR","TamR"))

ggplot(data=plot4,aes(x = response, y = EPAS1, fill = response)) +
  geom_violin(trim=F)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(plot.title = element_text(size = 14, family = "Arial", face = "bold",hjust = 0.5),
        text = element_text(size = 12, family = "Arial"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,face = "bold",colour = "black"),
        legend.position = "")+theme_bw()+
  labs(title="",x="", y = "Log10(EPAS1)+1")+
  ylim(0.5,4.5)
                  
