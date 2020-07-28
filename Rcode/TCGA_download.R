#TCGA analysis 

setwd("/home/eunji/proj/TCGA/")
rm(list = ls())

library(TCGAbiolinks)
library(DESeq2)
library(stringr) #str_subset
library(maftools)
library(dplyr)
library(readr)
library(tidyverse)

cancer  <- TCGAbiolinks:::getGDCprojects()$project_id
cancer <- str_subset(cancer, "TCGA") 
cancer_type <- sort(cancer) #33  types&proj
cancer <- strsplit(cancer_type,split = "-")
for(i in 1:length(cancer)){
  cancer[i] <- cancer[[i]][2]
}
cancer <- unlist(cancer) #33 cancer types 

cancer_type <- "TCGA-BRCA"

Download_TCGA<-function(cancer_type){
  suppressMessages(library(TCGAbiolinks)) 
  suppressMessages(library(SummarizedExperiment))
  suppressMessages(library(dplyr))
  suppressMessages(library(DT))
  dir= "/home/eunji/proj/TCGA/data/" #should change this before you run
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
  
  
  #download mutation data
  cancer <- strsplit(cancer_type,split = "-")
  cancer <- cancer[[1]][2]
  mut <- GDCquery_Maf(tumor = cancer,pipelines = "varscan2") 
  mut <- read.maf(maf = mut) #read MAF file 
  mut.mutsig.corrected=prepareMutSig(maf = mut)
  write.table(mut.mutsig.corrected,file = paste0(cancer_type,".maf"),quote = F,sep = "\t",row.names = F)
  
  #download Copy Number Variation data
  CNV_query <- GDCquery(project = cancer_type, 
                        data.category = "Copy Number Variation", 
                        data.type = "Copy Number Segment")
  
  GDCdownload(CNV_query)
  CNV_expdat <- GDCprepare(query = CNV_query)
  CNV_count_matrix=as.data.frame(CNV_expdat)
  write.csv(CNV_count_matrix,file = paste(cancer_type,"_CNV",".csv"))
  
  #download methylation
  meth_query <- GDCquery(project =cancer_type,
                         legacy = TRUE,
                         data.category = "DNA methylation",
                         platform = "Illumina Human Methylation 450")
  GDCdownload(meth_query,method = "api", files.per.chunk = 150)
  meth_expdat <- GDCprepare(query = meth_query)
  meth_count_matrix=assay(meth_expdat)
  write.csv(meth_count_matrix,file = paste(cancer_type,"_methylation2",".csv"))
  ####download miR data
  miR_query <- GDCquery(project = cancer_type, 
                        data.category = "Transcriptome Profiling", 
                        data.type = "miRNA Expression Quantification", 
                        workflow.type = "BCGSC miRNA Profiling")
  GDCdownload(miR_query)
  miR_expdat <- GDCprepare(query = miR_query)
  miR_expdat_matrix=assay(miR_expdat)
  write.csv(miR_expdat_matrix,file = paste(cancer_type,"_miRNAs2",".csv"))
  message(paste0(cancer_type," Download Finished!"))
}

Download_TCGA(cancer_type)
