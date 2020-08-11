# Tumor mutational burden is a determinant of 
# immune-mediated survival in breast cancer
# https://doi.org/10.1080/2162402X.2018.1490854

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(rtracklayer) #gtf file 
devtools::install_github("mariodeng/FirebrowseR")
library(FirebrowseR)
BiocManager::install("maftools")
library(maftools)

setwd("/home/eunji/proj/200722_TMB_BRCA/")

#========================================================
#data - TCGA & METABRIC clinical , RNA-seq , WES , CNA 

cancer_type <- "TCGA-BRCA"

Download_TCGA<-function(cancer_type){
  suppressMessages(library(TCGAbiolinks)) 
  suppressMessages(library(SummarizedExperiment))
  suppressMessages(library(dplyr))
  suppressMessages(library(DT))
  suppressMessages(library(maftools))
  dir= "/home/eunji/proj/200722_TMB_BRCA/data/TCGA/" #should change this before you run
  out_dir=paste0(dir,"/",cancer_type)
  dir.create(out_dir,recursive = T)
  setwd(out_dir)
  counts_query <- GDCquery(project = cancer_type, 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
  GDCdownload(counts_query)
  counts_expdat <- GDCprepare(query =counts_query)
  count_matrix= as.data.frame(assay(counts_expdat))
  CLC_query <- GDCquery(project = cancer_type, 
                        data.category = "Clinical", 
                        file.type = "xml")
  GDCdownload(CLC_query)
  write.csv(count_matrix,file = paste(cancer_type,"_Counts",".csv"))
  sample_NT <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("NT"))
  sample_TP <- TCGAquery_SampleTypes(barcode = colnames(count_matrix),typesample = c("TP"))
  write.csv(count_matrix[,sample_NT],file = paste(cancer_type,"_Counts_normal",".csv"))
  write.csv(count_matrix[,sample_TP],file = paste(cancer_type,"_Counts_tumor",".csv"))
  
  #clinical data
  CLC_query <- GDCquery(project = cancer_type, 
                        data.category = "Clinical", 
                        file.type = "xml")
  GDCdownload(CLC_query)
  clinical <- GDCprepare_clinic(CLC_query, clinical.info = "patient")
  write.csv(clinical,file = paste(cancer_type,"_clinical",".csv"))
  
  
  #survival data
  clinical_trait <- clinical  %>%
    dplyr::select(bcr_patient_barcode,gender,vital_status,                            
                  days_to_death,days_to_last_followup,race_list,
                  person_neoplasm_cancer_status,
                  stage_event_pathologic_stage,             
                  stage_event_tnm_categories  ) %>%
    distinct( bcr_patient_barcode, .keep_all = TRUE)
  
  
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
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=OS.Time/365)
  
  #combine clincial data
  survival_data <- rbind(dead_patient,alive_patient)
  write.csv(survival_data , file = paste(cancer_type,"_survival",".csv"))
  
  #download Copy Number Variation data
  CNV_query <- GDCquery(project = cancer_type, 
                        data.category = "Copy Number Variation", 
                        data.type = "Copy Number Segment")
  
  GDCdownload(CNV_query)
  CNV_expdat <- GDCprepare(query = CNV_query)
  CNV_count_matrix=as.data.frame(CNV_expdat)
  write.csv(CNV_count_matrix,file = paste(cancer_type,"_CNV",".csv"))
}

Download_TCGA(cancer_type)

#WES - oncocator (GDAC-oncocator file, 20160128)

library(stringr)
library(dplyr)
library(magrittr)
setwd("/home/eunji/proj/200722_TMB_BRCA/")
path <- "./data/TCGA/oncocator/"
filenames <- dir(path,recursive = T,pattern = ".txt",full.names = T)
data <- lapply(filenames,function(x){
  read.table(x,header = F,sep = "\t",stringsAsFactors = F,skip = 4,fill = T,quote = "")})
onco <- Reduce(function(x,y)rbind(x,y),data)
colnames <- colnames(read.table("./data/TCGA/oncocator/TCGA-A1-A0SB-01.hg19.oncotator.hugo_entrez_remapped.maf.txt",sep = "\t",skip = 3,header = T,fill = T,stringsAsFactors = F))
colnames(onco) <- colnames
save(onco,file = "./data/TCGA/oncocator/onco.RData")
write.table(onco,file = "./data/TCGA/oncocator/onco.maf",row.names = F,sep = "\t",quote = F)

#==========================================================
#data preprocessing 

#1 RNAseq : exclude male & metastasis tissue 

clin <- read.csv("./data/TCGA/TCGA-BRCA/TCGA-BRCA _clinical .csv",fill = T,sep = ",",header = T)

clin <- clin %>%
  dplyr::filter(gender == "FEMALE") 
clin = clin[-which(duplicated(clin$bcr_patient_barcode)),]

count <- read.csv("./data/TCGA/TCGA-BRCA/TCGA-BRCA _Counts .csv",sep = ",",fill = T,stringsAsFactors = F,header = T)
names(count)[1] <- 'gene_id'
load("~/proj/0_sh/ref/rna/hg38_ID.RData")
count <- merge(ID,count,by="gene_id")
count$gene_id <- NULL
count = count[-which(duplicated(count$gene_name)),]
rownames(count) <- count[,1]
count <- count[,-1]
colnames(count) <- str_sub(colnames(count),1,12)
colnames(count) <- gsub('.','-',colnames(count),fixed = T)
tcount <- as.data.frame(t(count))
rownames(tcount) <- gsub('.','-',rownames(tcount),fixed = T)
sample <- readxl::read_excel("./source_data/eSupplement_2_06_03_18.xlsx",
                             sheet = "eWorksheet 1",
                             range = "A4:A1086",
                             col_names = F)
sample$...1 <- str_sub(sample$...1,1,12) 
sample <- unique(sample$...1)
tcount <- tcount[c(rownames(tcount) %in% sample),]
count <- as.data.frame(t(tcount))
#=====================================================================
#maftools

library(maftools)
mut <- read.maf(maf = "./data/TCGA/oncocator/onco.maf")
summary <- getSampleSummary(mut)
summary$Tumor_Sample_Barcode <- str_sub(summary$Tumor_Sample_Barcode,1,12)
summary = summary[-which(duplicated(summary$Tumor_Sample_Barcode)),]
summary <- summary %>%
  dplyr::mutate(TMB=total/38) %>%
  dplyr::mutate(TMB_group=ifelse(TMB>mean(TMB),'high','low')) %>%
  dplyr::filter(TMB < 30) %>%
  dplyr::mutate(Breast_Tumors=1:nrow(summary))
gene_summary <- getGeneSummary(mut)

datatable(getSampleSummary(mut),
          filter = 'top',
          options = list(scrollX=T,keys=T,pageLength=5),
          rownames = F)

plotmafSummary(maf = mut, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE,
               titvRaw = F,top = 10)
oncoplot(maf = mut, top = 20, removeNonMutated = TRUE)
oncostrip(maf = mut)
oncostrip(maf = mut,genes=c("TP53","TTN","KMT2D"))
library(ggplot2)
BiocManager::install("ggridges")
library(ggridges)
ggplot(summary,aes(x=Breast_Tumors,y=TMB,height=TMB,fill=TMB_group))+
  geom_ridgeline() +
  ggtitle("TMB distribution")
  
OncogenicPathways(maf = mut) #pathway
PlotOncogenicPathways(maf = mut,pathways = "WNT")
mutsig <- prepareMutSig(maf=mut) #prepare mutsig file

#=======================================================
#CIBERSORT - tumor infiltrting lymphocytes(TIL)
source("./code/CIBERSORT.R")
LM22 <- read.delim("./data/CIBERSORT/LM22.txt")
exp <- count[rownames(count)%in%LM22$Gene.symbol,]
sig_matrix <- "./data/CIBERSORT/LM22.txt"

#TPM normalization 

gene_length <- read.table("./data/hg38_gene_length.txt",sep = "\t",header = F)
names(gene_length)[1] <- 'gene_id'
names(gene_length)[2] <- 'gene_length'
ID <- merge(ID,gene_length,by="gene_id")
ID$gene_id <- NULL

tpm <- function(counts,lengths){
  rpk <- counts/(lengths/1000)
  coef <- sum(rpk)/1e6
  rpk/coef
}
exp <- data.frame(gene_name=rownames(exp),exp)
rownames(exp) <- NULL
merge <- merge(ID,exp,by="gene_name")
genes <- data.frame(merge[,1:2])
counts <- data.frame(row.names=merge$gene_name,merge[,3:1078])
tpms <- apply(counts,2,function(x)tpm(x,genes$gene_length))
colSums(tpms) #1e+06
tpm_exp <- as.data.frame(tpms)
colnames(tpm_exp) <- gsub('.','-',colnames(tpm_exp),fixed = T)
tpm_exp <- data.frame(gene_name=rownames(tpm_exp),tpm_exp)
rownames(tpm_exp) <- NULL
write.table(tpm_exp,file = "./result/ciber_matrix_tpm.txt",sep = "\t",quote = F)
mixture_file = "./result/ciber_matrix_tpm.txt"
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=10, QN=F)
res <- as.data.frame(res_cibersort)
library(dplyr)
library(tidyr)
library(tibble)
res1 <- res %>%
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols=2:23,
               names_to= "celltype",
               values_to = "Proportion")


library(ggplot2)
ggplot(res1,aes(sample,Proportion,fill = celltype)) + 
  geom_bar(position = "stack",stat = "identity")+
  theme_bw()+
  guides(fill=guide_legend(ncol=1))
