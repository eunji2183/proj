#https://doi.org/10.1042/BSR20194337
#Mining TCGA database for tumor mutation burden 
#and their clinical significance in bladder cancer
#TCGA- BLCA -TMB -survival -DEG - GO & KEGG - PPI network 

rm(list = ls())
setwd("/home/eunji/proj/200723_BLCA/")
BiocManager::install("ComplexHeatmap")
BiocManager::install("NMF")
BiocManager::available("BSgenome")

#somaticsignatures
BiocManager::install("SomaticSignatures")
BiocManager::install("SomaticCancerAlterations")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
BiocManager::install("VariantAnnotation")
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(VariantAnnotation)

library(BSgenome)
library(NMF)
library(stringi)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(maftools) #read.maf
library(stringr)
library(DESeq2)
library(pheatmap)
library(ComplexHeatmap)
library(clusterProfiler)
library(ggplot2)
library(survminer)
library(survival)



#======================================================
# TCGA : rnaseq , clinical data 

cancer_type <- "TCGA-BLCA"

Download_TCGA<-function(cancer_type){
  suppressMessages(library(TCGAbiolinks)) 
  suppressMessages(library(SummarizedExperiment))
  suppressMessages(library(dplyr))
  suppressMessages(library(DT))
  dir= "/home/eunji/proj/200723_BLCA/data/" #should change this before you run
  out_dir=paste0(dir,"/",cancer_type)
  dir.create(out_dir,recursive = T)
  setwd(out_dir)
  #get GDC version information
  gdc_info = getGDCInfo()
  #RNAseq HTSeq-Counts
  counts_query <- GDCquery(project = cancer_type, 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
  GDCdownload(counts_query)
  counts_expdat <- GDCprepare(query =counts_query)
  count_matrix= as.data.frame(assay(counts_expdat))
  write.csv(count_matrix,file = paste(cancer_type,"_Counts",".csv"))
  sample_NT <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("NT"))
  sample_TP <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("TP"))
  write.csv(count_matrix[,sample_NT],file = paste(cancer_type,"_Counts_normal",".csv"))
  write.csv(count_matrix[,sample_TP],file = paste(cancer_type,"_Counts_tumor",".csv"))
  
  ##origanize the clinical data (shuould do some ajustion maybe try next time)
  
  CLC_query <- GDCquery(project = cancer_type, 
                        data.category = "Clinical", 
                        file.type = "xml")
  GDCdownload(CLC_query)
  clinical <- GDCprepare_clinic(CLC_query, clinical.info = "patient")
  clinical_trait <- clinical  %>%
    dplyr::select(bcr_patient_barcode,gender,vital_status,
                  age_at_initial_pathologic_diagnosis,
                  neoplasm_histologic_grade,
                  days_to_death,days_to_last_followup,race_list,
                  person_neoplasm_cancer_status,
                  stage_event_pathologic_stage,             
                  stage_event_tnm_categories  ) %>%
    dplyr::distinct( bcr_patient_barcode, .keep_all = TRUE)
  #organize dead data
  dead_patient <- clinical_trait  %>%
    dplyr::filter(vital_status == 'Dead') %>%
    dplyr::select(-days_to_last_followup) %>%
    reshape::rename(c(bcr_patient_barcode = 'Barcode',
                      gender = 'Gender',
                      vital_status = 'OS',
                      days_to_death='OS.Time',
                      race_list = 'Race',
                      person_neoplasm_cancer_status='cancer_status',
                      age_at_initial_pathologic_diagnosis = 'Age',
                      neoplasm_histologic_grade = 'Grade',
                      stage_event_pathologic_stage = 'Stage',
                      stage_event_tnm_categories = 'TNM' )) %>%
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=OS.Time/365)
  
  #organize the alive data
  alive_patient <- clinical_trait %>%
    dplyr::filter(vital_status == 'Alive') %>%
    dplyr::select(-days_to_death) %>%
    reshape::rename(c(bcr_patient_barcode = 'Barcode',
                      gender = 'Gender',
                      vital_status = 'OS',
                      days_to_last_followup='OS.Time',
                      race_list = 'Race',
                      person_neoplasm_cancer_status='cancer_status',
                      age_at_initial_pathologic_diagnosis = 'Age',
                      neoplasm_histologic_grade = 'Grade',
                      stage_event_pathologic_stage = 'Stage',
                      stage_event_tnm_categories = 'TNM' )) %>%
    dplyr::mutate(OS=ifelse(OS=='Dead',1,0))%>%
    dplyr::mutate(OS.Time=OS.Time/365)
  
  #combine clincial data
  survival_data <- rbind(dead_patient,alive_patient)
  write.csv(survival_data , file = paste(cancer_type,"_clinical",".csv"))
  
  #download mutation data
  cancer <- strsplit(cancer_type,split = "-")
  cancer <- cancer[[1]][2]
  mut <- GDCquery_Maf(tumor = cancer,pipelines = "varscan2") 
  mut <- read.maf(maf = mut) #read MAF file 
  mut.mutsig.corrected=prepareMutSig(maf = mut)
  write.table(mut.mutsig.corrected,file = paste0(cancer_type,".maf"),quote = F,sep = "\t",row.names = F)
}

Download_TCGA(cancer_type)  

#=========================================================
  
#mutation data : TCGA MAF 


mut <- read.maf(maf = "./data/TCGA-BLCA/TCGA-BLCA.maf") #read MAF file 

summary <- getSampleSummary(mut)
gene_summary <- getGeneSummary(mut)

datatable(getSampleSummary(mut),
          filter = 'top',
          options = list(scrollX=T,keys=T,pageLength=5),
          rownames = F)

plotmafSummary(maf = mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,
               titvRaw = F,top = 10)

oncoplot(maf = mut, top = 20, removeNonMutated = TRUE)
oncostrip(maf = mut)
oncostrip(maf = mut,genes=c("TP53","TTN","KMT2D"))

titv = titv(maf = mut, plot = FALSE, useSyn = TRUE)
plotTiTv(res = titv)

somaticInteractions(maf = mut,top = 25,pvalue = c(0.05,1))

lollipopPlot(maf = mut,gene = "TP53",
             #AACol='Protein_Change',
             showMutationRate = T)
rainfallPlot(maf = mut,detectChangePoints = T,pointSize = 0.6)

BLCA.mutload=tcgaCompare(maf = mut,cohortName = "Example-BLCA")

tnm = trinucleotideMatrix(maf=mut,ref_genome = "/home/eunji/proj/0_sh/ref/rna/Homo_sapiens_assembly38.fasta",add = T,useSyn = T)

sign=extractSignatures(mat = tnm,nTry=6,plotBestFitRes = F)

#======================================================
#TMB : total mutation /38 

t_mut <- mut@variants.per.sample
t_mut <- t_mut %>%
  dplyr::mutate(TMB=Variants/38) %>%
  dplyr::mutate(TMB_group=ifelse(TMB>median(TMB),'high','low'))

coldata <- t_mut %>%
  dplyr::select(Tumor_Sample_Barcode,TMB_group) %>% distinct()


#=======================================================
#DEG 

count <- read.table("./data/TCGA-BLCA/TCGA-BLCA-Counts_tumor-.csv",sep = ",",row.names = 1,header = T,stringsAsFactors = F)
count <- as.data.frame(t(count))
rownames(count) <- gsub('.','-',rownames(count),fixed = T)
count <- data.frame(sampleID=rownames(count),count)
rownames(count) <- NULL
count$sampleID <- str_sub(count$sampleID,1,15)
names(coldata)[1] <- 'sampleID'
coldata$sampleID <- str_sub(coldata$sampleID,1,15)

TMB_count <- merge(coldata,count,by="sampleID") %>% 
  distinct() %>% 
  as.data.frame() %>%
  dplyr::arrange(TMB_group)
TMB_count = TMB_count[-which(duplicated(TMB_count$sampleID)),]

rownames(TMB_count) <- TMB_count[,1]
TMB_count <- TMB_count[,-1]

#high(N=203) vs low(N=204)
coldata <- data.frame(row.names=rownames(TMB_count),TMB_group=TMB_count$TMB_group)
coldata$TMB_group <- factor(coldata$TMB_group,levels=c("high","low"))

TMB_count <- as.data.frame(t(TMB_count))
TMB_count <- TMB_count[-1,]

TMB_count <- na.omit(TMB_count)
sapply(TMB_count,mode)
sapply(TMB_count,class)
TMB_count[,1:407] <- sapply(TMB_count[,1:407],as.numeric)

dds <- DESeqDataSetFromMatrix(countData = TMB_count,colData = coldata,design = ~TMB_group)
dds <- DESeq(dds)
res <- results(dds,contrast = c("TMB_group","high","low")) 
resultsNames(dds) 
summary(res)
table(res$padj<0.05)
res <- as.data.frame(res)
res <- res %>% dplyr::filter(!is.na(log2FoldChange))
gtf <- rtracklayer::import("/home/eunji/proj/0_sh/ref/rna/Homo_sapiens.GRCh38.96.gtf")
gtf <- as.data.frame(gtf)
ID <- gtf %>%
  dplyr::select(gene_id,gene_name) %>% distinct()
save(ID,file="./data/ID.RData")
res <- data.frame(gene_id=rownames(res),res)
rownames(res) <- NULL
res <- merge(ID,res,by="gene_id")
res <- res[order(res$padj),]
save(res,file = "./result/res.RData")

#padj<0.05 & |logFC| > 0.5

subset <- res %>%
  dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.5) %>%
  dplyr::select(gene_id,gene_name)


norm_count <- as.data.frame(counts(dds,normalized=TRUE))
norm_count <- data.frame(gene_id=rownames(norm_count),norm_count)
rownames(norm_count) <- NULL

sig <- merge(subset,norm_count,by="gene_id")
sig <- sig %>%
  dplyr::select(-gene_id)
rownames(sig) <- sig[,1]
sig <- sig[,-1]
colnames(sig) <- gsub('.','-',colnames(sig),fixed = T)

#pheatmap 
pheatmap(sig,annotation_col = coldata,scale = "row", clustering_distance_row = "correlation",
         show_colnames = T,show_rownames = T,
         cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         border_color = 'white',
         display_numbers = F,main = "",key=T)
#======================================================
#GO & KEGG 

gene <- rownames(sig)
gene <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
de <- gene$ENTREZID
go <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db",ont = "all")
dotplot(go,split ="ONTOLOGY") +facet_grid(ONTOLOGY~.,scale ="free")

KEGG <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05)
dotplot(KEGG)
test <- data.frame(KEGG)
browseKEGG(KEGG, 'hsa04060')

#======================================================
#survival 
setwd("/home/eunji/proj/200723_BLCA/")
clinical <- read.table("./data/TCGA-BLCA/TCGA-BLCA _clinical .csv",sep = ",",fill = T,header = T,row.names = 1)
TMB <- data.frame(Barcode=rownames(coldata),coldata)
TMB$Barcode <- str_sub(TMB$Barcode,1,12)
rownames(TMB) <- NULL
survival <- merge(clinical,TMB,by="Barcode")
sfit <- survfit(Surv(OS.Time,OS)~TMB_group,data = survival)
ggsurvplot(sfit,pval = TRUE)




#======================================================
#somaticsignatures
BiocManager::install("SomaticSignatures")
BiocManager::install("SomaticCancerAlterations")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
BiocManager::install("VariantAnnotation")
BiocManager::install("GenomicRanges")
library(SomaticSignatures)
library(SomaticCancerAlterations)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(VariantAnnotation)
library(GenomicRanges)

sca_metadata = scaMetadata()
sca_data = unlist(scaLoadDatasets())
sca_data$study = factor(gsub("(.*)_(.*)", "\\1", toupper(names(sca_data))))
sca_data = unname(subset(sca_data, Variant_Type %in% "SNP"))
sca_data = keepSeqlevels(sca_data, hsAutosomes())


