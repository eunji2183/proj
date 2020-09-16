#excel-pheatmap 
setwd("/home/eunji/R/project/200717_heatmap")

library(pheatmap)

normalized <- readxl::read_excel("./data/200324 heatmap (intensity log and RNA copy number) (1).xlsx",
                                   sheet = "Sheet1",
                                   range = "B2:F24", 
                                   col_names = T)
normalized <- as.data.frame(normalized)
rownames(normalized) <- normalized[,1]
normalized <- normalized[,-1]
pheatmap(normalized,scale = "row", clustering_distance_row = "correlation",
                 show_colnames = T,show_rownames = T,
                 cluster_cols = T,cluster_rows = T,
                 color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                 border_color = 'white',
                 display_numbers = F,main = "",key=T)

raw <- readxl::read_excel("./data/200324 heatmap (intensity log and RNA copy number) (1).xlsx",
                                 sheet = "Sheet1",
                                 range = "B2:j24", 
                                 col_names = T)
raw <- as.data.frame(raw)
rownames(raw) <- raw[,1]
raw <- raw[,-1]
raw <- raw[,5:8]
colnames(raw) <- c("control1","control2","UPM1","UPM2")
pheatmap(raw,scale = "row", clustering_distance_row = "correlation",
         show_colnames = T,show_rownames = T,
         cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         border_color = 'white',
         display_numbers = F,main = "",key=T)
save(list = ls(),file="./data/heatmap.RData")

---------------------------------------------------------------------------
#microdust heatmap 
library(magrittr)
library(dplyr)
BiocManager::install("TCC")

count=read.csv("./featurecounts/count.csv",header = T,row.names = 1)
count <- count[,c(1,2,3,5,6)]
names(count)[2] <- 'Control1'
names(count)[3] <- 'Control2'
names(count)[4] <- 'PM1'
names(count)[5] <- 'PM2'
rownames(count) <- count[,1]
count <- count[,-1]
group <- rep(c("con","upm"),c(2,2))
gtf <- rtracklayer::import("/home/eunji/ref/Homo_sapiens.GRCh38.96.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()
gene <- c("TNF","IL12RB1","CCR6","FPR1",
          "EIF2AK3","ERN1","ATF3","HSP90B1","CALR","XBP1","CYP1A1","CRELD2",
          "NFATC2","OR2A4","CALM2A",
          "F2RL3","NQO1","TXNRD1","GCLM","MEF2C","GRIA4","SIK1",
          "AHR","AHRR")
gene <- data.frame(gene_name=gene)
res <- DESeq2.1(count = count,group = group,ID=ID)
heat <- res[,c(1,2,9,10,11,12)]
heat <- merge(gene,heat,by="gene_name")
heat <- heat[-3,]
heat <- heat[,-2]
heat <- heat[c(19,14,10,7,8,3,13,4,5,6,21,16,9,17,20,11,15,12,18,1,2),]
rownames(heat) <- heat[,1]
heat <- heat[,-1]
DESeq2.2 <- function(count,group,ID){
  suppressMessages(library(TCC))
  suppressMessages(library(DESeq2))
  suppressMessages(library(dplyr))
  tcc <- new("TCC",count,group)
  tccNF <- calcNormFactors(tcc, method='tmm', test.method='deseq2', 
                           iteration = 3, FDR = 0.05, floorPDEG = 0.05)
  tccDE <- estimateDE(tccNF, test.method = 'deseq2', FDR = 0.5)
  tccRES <- getResult(tccDE, sort = T)
  tccRES <- merge(ID,tccRES,by="gene_id")
  tccRES <- dplyr::arrange(tccRES,rank)
  return(tccRES)
}
res2 <- DESeq2.2(count = count,group = group,ID=ID)
library(pheatmap)
annotation_row = data.frame(row.names = rownames(heat),
  GeneGroup=factor(rep(c("Inflammation","ER stress","Ca2+ signaling","Melanogenesis","Aryl hydrocarbon"),
                       c(3,8,1,7,2))))
pheatmap(heat,scale="row",color = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize_col = 12, fontsize_row=12,
         cellwidth = 17, cellheight = 17,show_rownames = T,show_colnames = T,legend = T,
         border="white",cluster_cols = T,cluster_rows = F,
         annotation_row = annotation_row,annotation_legend = T)

