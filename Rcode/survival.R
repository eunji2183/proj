working_dir="/home/eunji/R/colon/clinical/"
setwd(working_dir)
rm(list = ls())
options(stringsAsFactors = F) #no character2factor
BiocManager::install("RTCGAToolbox")
BiocManager::install("ggplot2")
BiocManager::install("survival")
BiocManager::install("survminer")
BiocManager::install("dplyr")
BiocManager::install("asaur")
install.packages("magrittr")
install.packages("ffbase")
install.packages("UCSCXenaTools",dependencies = TRUE) 
library(magrittr)
library(dplyr)
library(survminer)
library(survival)
library(RTCGAToolbox) 
library(ggplot2)
library(readr)
library(ffbase)
library(data.table)
library(UCSCXenaTools)
library(asaur)
###################################################################################
#UCSC Xena - COADREAD_clinicalMatrix 
getTCGAdata(project = 'COADREAD') 
phenocoad=getTCGAdata(project = 'COADREAD',download = TRUE,forceDownload = TRUE)
coad_clinical=XenaPrepare(phenocoad) 
head(coad_clinical)

#Firehose COADREAD data
getFirehoseDatasets()
getFirehoseRunningDates() 
COADData <- getFirehoseData(dataset="COADREAD",runDate = "20150204",forceDownload = TRUE,clinical = TRUE,Mutation = TRUE,mRNAArray = TRUE)
#gdac-COADREAD.clinical data 
clindata=read.csv("/home/eunji/R/colon/clinical/gdac.broadinstitute.org_COADREAD.Merge_Clinical.Level_1.2016012800.0.0/COADREAD.clin.merged.csv",header = F,stringsAsFactors = F)
tclindata=t(clindata)  #col&row exchange 
dt <- as.data.frame(tclindata)
dt1 <- dt[,c(17,18,5,1:4,6:16,19:3545)] #row exchange 
####################################################################################
#isoform data -log(tpm+0.001)
isodata <- fread("/home/eunji/R/colon/clinical/TcgaTargetGtex_rsem_isoform_tpm")
isodata <- as.data.frame(isodata)
aimp2_dx2 <- filter(isodata,sample == "ENST00000395236.2")
rownames(aimp2_dx2) <- aimp2_dx2[,1]
aimp2<- aimp2_dx2[,-1]
taimp2 <- t(aimp2) # col&row exchange 
taimp2_2 <- as.data.frame(taimp2) 
write.table(taimp2,"/home/eunji/R/colon/taimp2.txt",sep = "\t")
taimp2_3 <- read.csv("/home/eunji/R/colon/clinical/taimp2.csv",header = T,stringsAsFactors = F)
str(taimp2_3)
merge1 <- merge(x=coad_clinical,y=taimp2_3,by='sampleID')
merge2 <- merge1[,c(1,134,2:133)]
write.csv(merge2,"/home/eunji/R/colon/aimp2.csv",sep = ",") 
#######################################################################################
#aimp2_dx2 isoform expression level, tpm values are extracted, log2(x+0.001) transformed 
aimp2_exp<- merge2[,1:2]
tumor <- aimp2_exp[grep("0.$",aimp2_exp$sampleID),]
normal <- aimp2_exp[grep("1.$",aimp2_exp$sampleID),]
tumor1 <- mutate(tumor, group = "tumor(n=383)")
normal1 <- mutate(normal, group = "normal(n=51)")
exp <- rbind(tumor1,normal1)
exp$group <- as.factor(exp$group)
p <- ggboxplot(exp,x='group',y='ENST00000395236.2',fill = 'group',bxp.errorbar = T,bxp.errorbar.width = 0.2)
p+stat_compare_means()
#######################################################################################
#aimp2_dx2 survival 
clinical <- read.csv("/home/eunji/R/colon/clinical/aimp2.csv",header = T,stringsAsFactors = F)
clinical1 <- clinical[grep("0.$",clinical$sampleID),]
clinical_trait <- clinical1 %>%
  dplyr::select(bcr_patient_barcode,ENST00000395236.2,gender,vital_status,
                days_to_death,days_to_last_followup,person_neoplasm_cancer_status,
                pathologic_stage,pathologic_T,pathologic_N,pathologic_M,radiation_therapy,history_of_neoadjuvant_treatment) %>%
  distinct(bcr_patient_barcode,.keep_all=TRUE)
dead_patient <- clinical_trait %>%
  dplyr::filter(vital_status == 'DECEASED') %>%
  dplyr::filter(radiation_therapy == 'NO') %>%
  dplyr::select(-days_to_last_followup) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_death = 'OS.Time',
                    person_neoplasm_cancer_status = 'cancer_status',
                    pathologic_stage = 'stage')) %>%
  mutate(OS=ifelse(OS=='DECEASED',1,0)) %>%
  mutate(OS.Time=OS.Time/365) %>%
  mutate(stage_group=as.character(factor(stage,levels = c("Stage I","Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                          "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                          "Stage IV","Stage IVA","Stage IVB"),
                                         labels = c(1,2,2,2,2,3,3,3,3,4,4,4))))
alive_patient <- clinical_trait %>%
  dplyr::filter(vital_status == 'LIVING') %>%
  dplyr::filter(radiation_therapy == 'NO') %>%
  dplyr::select(-days_to_death) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_last_followup='OS.Time',
                    person_neoplasm_cancer_status='cancer_status',
                    pathologic_stage = 'stage')) %>%
  mutate(OS=ifelse(OS=='DECEASED',1,0))%>%
  mutate(OS.Time=OS.Time/365) %>%
  mutate(stage_group=as.character(factor(stage,levels = c("Stage I","Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                          "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                          "Stage IV","Stage IVA","Stage IVB"),
                                         labels = c(1,2,2,2,2,3,3,3,3,4,4,4))))
survival_data <- rbind2(dead_patient,alive_patient) 
group <- ifelse(survival_data$ENST00000395236.2>median(survival_data$ENST00000395236.2),'aimp2-dx2_high','aimp2-dx2_low')
sfit <- survfit(Surv(OS.Time,OS)~group,data = survival_data)
ggsurvplot(sfit,conf.int = FALSE,pval = TRUE)
cfit <- coxph(Surv(OS.Time,OS)~group,data = survival_data) 
table(survival_data$OS)
write.csv(survival_data,"/home/eunji/R/colon/clinical/survival_data.csv")
#############################################################################
#quantile survival analysis 
prob=c(0.25,0.75)
upper=quantile(survival_data$ENST00000395236.2,prob[2],names = FALSE)
lower=quantile(survival_data$ENST00000395236.2,prob[1],names = FALSE)
surv <- survival_data %>%
  mutate(iso_group =  ifelse(ENST00000395236.2>upper,"upper",
                             ifelse(ENST00000395236.2>lower,"middle","lower"))) %>%
  filter(iso_group == "upper" | iso_group == "lower")
sfit <- survfit(Surv(OS.Time,OS)~iso_group,data = surv)
ggsurvplot(sfit,pval = TRUE)
table(surv$iso_group) 
###############################################################################
# stage specificity survival 
fit <- coxph(Surv(OS.Time,OS)~Gender+stage_group+ENST00000395236.2,data = survival_data)
fit
stage <- survfit(Surv(OS.Time,OS)~iso_group+strata(surv$stage_group),data = surv)
ggsurvplot(stage,pval = TRUE)
stage1 <- surv %>% dplyr::filter(stage_group == 1)
stage2 <- surv %>% dplyr::filter(stage_group == 2)
stage3 <- surv %>% dplyr::filter(stage_group == 3)
stage4 <- surv %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(OS.Time,OS)~iso_group,data = stage1)
stage2fit <- survfit(Surv(OS.Time,OS)~iso_group,data = stage2)
stage3fit <- survfit(Surv(OS.Time,OS)~iso_group,data = stage3)
stage4fit <- survfit(Surv(OS.Time,OS)~iso_group,data = stage4)
ggsurvplot(stage1fit,pval = TRUE,risk.table = TRUE)
surv_object <- Surv(time = gene_surv$OS.Time,event = gene_surv$OS)
fit1 <- survfit(formula = surv_object ~ gene_surv$quantile, data = gene_surv)
ggsurvplot(fit1,data = gene_surv,pval = TRUE,conf.int = TRUE)
surv_object2 <- Surv(time = survival_data$OS.Time,event = survival_data$OS)
fit2 <- survfit(formula = surv_object ~ survival_data$stage ,data = survival_data)
sfit <- survfit(Surv(times, patient.vital_status)~admin.disease_code, data=clin)
summary(sfit, times=seq(0,365*5,365))



################################################################################
##비례위험가정 만족여부 
#kaplan-curve (cross-over 확인)
sfit <- survfit(Surv(OS.Time,OS)~quantile,data = gene_surv)
ggsurvplot(sfit,conf.int = FALSE,pval = TRUE)
#scheonfeld method (coxph)
scheon<- coxph(Surv(OS.Time,OS)~quantile,data = gene_surv)
summary(scheon)
cox.zph(scheon)
plot(cox.zph(scheon))
##non-proportional hazard 의 교정 
iso <- survSplit(Surv(OS.Time,OS)~.,data = gene_surv,cut = c(2.00000000,4.10000000),episode = "tgroup")


#log-log plot - 그룹간에 평행이 나오면 만족 
plot(sfit,fun="cloglog", lty=1:2, col=c("Black", "Grey50"), lwd=2, font=2, ylim=c(-3,2), font.lab=2, main="log-log KM curves by Rx", ylab="log-log survival", xlab="time(year)")
legend("topright",lty=1:2,legend=c("Group1","Group2"),col=c("Black", "Grey50"))
#expected plot 
sfit <- survfit(Surv(OS.Time,OS)~quantile,data = gene_surv)
plot(sfit,lty="dashed", col=c("Black", "Grey50"), lwd=2, font=2, font.lab=2, main="Observed Versus Expected Plots by Clinic", ylab="Survival probability", xlab="Survival time")
par(new=T)   
expected <- coxph(Surv(OS.Time,OS)~quantile,data = gene_surv)
summary(expected)
kmfit.exp <- survfit(expected,data=gene_surv) 
plot(kmfit.exp,lty="solid", col=c("Black", "Grey50"), lwd=2, font=2, font.lab=2)
     
     
    
