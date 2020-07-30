#nchar('NGGCTCCATGGATGAGAAGAAGAAGGAGTGGTAACAGTCACAGATCCCTACCTGCCTGGCTAAGACCCGTGGCCGTCAAGGACTGGTTCGGGGTGGATTCA')

source('functions.R')
source("libs.R")


fls <- dir("kallisto/output/", recursive = T, pattern = ".tsv", full.names = T)
lapply(fls, function(fl){
  read_delim(fl, "\t", escape_double = FALSE, trim_ws = TRUE)
}) -> FLs

SampleIDs <- gsub(pattern = "/abundance.tsv", replacement = "", fls) %>% gsub(pattern = ".*/", replacement = "")
names(FLs) <- SampleIDs

transcripts_to_genes <- read_delim("kallisto/index/homo_sapiens/transcripts_to_genes.txt", 
                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                   trim_ws = TRUE)

DF <- do.call(cbind, FLs)
dim(DF)
rownames(DF) <- DF$CON_1.target_id
df_cnt <- DF %>% select(contains("est_counts")) 
df_tpm <- DF %>% select(contains("tpm"))
df_cnt <- df_cnt[rowSums(df_cnt) > 0, ]
df_tpm <- df_tpm[rowSums(df_tpm) > 0, ]

dim(df_tpm)
df_tpm$Trans_ID <- rownames(df_tpm)
df_tpm <- merge(df_tpm, transcripts_to_genes[, c(1,3)], by.x = 'Trans_ID', by.y = "X1")
df_tpm_mt <- as.matrix(df_tpm[, -c(1, ncol(df_tpm))])
rownames(df_tpm_mt) <- df_tpm$X3
df_tpm_GN <- rmdup(df_tpm_mt) %>% as.data.frame()
class(df_tpm_GN)



### EdgeR ###
GRP <- c('con', 'con', 'con', 'con', 'upm', 'upm', 'upm', 'upm')
CNT <- df_tpm_GN
tcc <- new("TCC", CNT, GRP)
tccNF <- calcNormFactors(tcc, method='tmm', test.method='edger', iteration = 3, FDR = 0.05, floorPDEG = 0.05)
tccDE <- estimateDE(tccNF, test.method = 'edger', FDR = 0.5)
tccRES <- getResult(tccDE, sort = T)
#View(tccRES)

table(tccRES$estimatedDEG == 1)


ggplot(tccRES, aes(a.value, m.value, col = estimatedDEG)) +
  geom_point(alpha = 0.5) +
  scale_color_gradientn(colours = c("gray", 'red')) +
  theme_bw() +
  ggtitle("DEG(EdgeR) plot")

write.csv(tccRES, "DEGresults/EdgeR_table.csv")



### DESeq2 ###
GRP <- c('con', 'con', 'con', 'con', 'upm', 'upm', 'upm', 'upm')
CNT <- df_tpm_GN
tcc <- new("TCC", CNT, GRP)
tccNF <- calcNormFactors(tcc, method='tmm', test.method='deseq2', iteration = 3, FDR = 0.05, floorPDEG = 0.05)
tccDE <- estimateDE(tccNF, test.method = 'deseq2', FDR = 0.5)
tccRES <- getResult(tccDE, sort = T)
#View(tccRES)

table(tccRES$estimatedDEG == 1)


ggplot(tccRES, aes(a.value, m.value, col = estimatedDEG)) +
  geom_point(alpha = 0.5) +
  scale_color_gradientn(colours = c("gray", 'red')) +
  theme_bw() +
  ggtitle("DEG(DESeq2) plot")

write.csv(tccRES, "DEGresults/DESeq2_table.csv")



#### HEATMAP ####
deseq2 <- read.csv("DEGresults/DESeq2_table.csv", row.names = 1)
edger <- read.csv("DEGresults/EdgeR_table.csv", row.names = 1)

gid <- deseq2 %>% filter(estimatedDEG == 1) %>% select(gene_id) %>% unlist %>% as.character()
png("DEGresults/heatmap_deseq2.png", width = 700, height = 700)
heatmap.2(log10(as.matrix(df_tpm_GN[gid,])+1), trace = 'none', 
          cexCol = 1, col = bluered(100), density.info = 'none', margins = c(10, 10), Colv = F)
dev.off()
heatmap.2(as.matrix(df_tpm_GN[gid,]), trace = 'none', cexCol = 1, col = bluered(100), density.info = 'none', margins = c(10, 10))
as.matrix(df_tpm_GN[gid,])

gid <- edger %>% filter(estimatedDEG == 1) %>% select(gene_id) %>% unlist %>% as.character()
png("DEGresults/heatmap_edger.png", width = 700, height = 700)
heatmap.2(log10(as.matrix(df_tpm_GN[gid,])+1), trace = 'none', 
          cexCol = 1, col = bluered(100), density.info = 'none', margins = c(10, 10), Colv = F)
dev.off()
heatmap.2(as.matrix(df_tpm_GN[gid,]), trace = 'none', cexCol = 1, col = bluered(100), density.info = 'none', margins = c(10, 10))
as.matrix(df_tpm_GN[gid,])


#### table for GSEA ####

head(df_tpm_GN)
colnames(df_tpm_GN) <- c(paste('CON', 1:4, sep = "_"), paste('UPM', 1:4, sep = '_'))
cbind(GeneSymbol=rownames(df_tpm_GN), df_tpm_GN) %>% head
write_delim(cbind(GeneSymbol=rownames(df_tpm_GN), df_tpm_GN), 'TPM_GSEAinput.csv', delim = "\t")


## read gsea results

gseaR <- read.csv("GSEAresults/EJgmts_191201/gseapy.gsea.gene_set.report.csv")
View(gseaR)
