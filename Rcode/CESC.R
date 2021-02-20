setwd("/home/eunji/proj/cervical/")
options(stringsAsFactors = F)
BiocManager::install("RTCGAToolbox")
BiocManager::install("ggplot2")
BiocManager::install("survival")
BiocManager::install("survminer")
BiocManager::install("dplyr")
BiocManager::install("asaur")
BiocManager::install("DESeq2")
install.packages("magrittr")
install.packages("ffbase")
install.packages("UCSCXenaTools",dependencies = TRUE) 
install.packages("reshape")
install.packages("ggpubr")
install.packages("survminer")
library(DESeq2)
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
library(reshape)
library(ggpubr)
library(stringr)
library(survminer)

####################################################################
#UCSC Xena - CESC_clinicalMatrix , GDAC-mRNAseq raw-count 
getTCGAdata(project = 'CESC') 
phenocesc=getTCGAdata(project = 'CESC',download = TRUE,forceDownload = TRUE)
CESC_clinical=XenaPrepare(phenocesc) 
clinical1 <- CESC_clinical[grep("0.$",CESC_clinical$sampleID),]
clinical_trait <- clinical1 %>%
  dplyr::select(bcr_patient_barcode,gender,vital_status,
                days_to_death,days_to_last_followup,person_neoplasm_cancer_status,
                clinical_stage,pathologic_T,pathologic_N,pathologic_M,radiation_therapy,history_of_neoadjuvant_treatment,human_papillomavirus_type) %>%
  distinct(bcr_patient_barcode,.keep_all=TRUE)
dead_patient <- clinical_trait %>%
  dplyr::filter(vital_status == 'DECEASED') %>%
  dplyr::filter(!is.na(clinical_stage)) %>%
  dplyr::select(-days_to_last_followup) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_death = 'OS.Time',
                    person_neoplasm_cancer_status = 'cancer_status',
                    clinical_stage = 'stage',
                    human_papillomavirus_type = 'HPV')) %>%
  mutate(OS=ifelse(OS=='DECEASED',1,0)) %>%
  mutate(OS.Time=OS.Time/365) %>%
  mutate(stage_group=as.character(factor(stage,levels = c("Stage I","Stage IA","Stage IA1","Stage IA2",
                                                          "Stage IB","Stage IB1","Stage IB2",
                                                          "Stage II", "Stage IIA","Stage IIA1","Stage IIA2",
                                                          "Stage IIB",
                                                          "Stage III", "Stage IIIA", "Stage IIIB", 
                                                          "Stage IVA","Stage IVB"),
                                         labels = c(1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4))))
alive_patient <- clinical_trait %>%
  dplyr::filter(vital_status == 'LIVING') %>%
  dplyr::filter(!is.na(clinical_stage)) %>%
  dplyr::select(-days_to_death) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_last_followup = 'OS.Time',
                    person_neoplasm_cancer_status = 'cancer_status',
                    clinical_stage = 'stage',
                    human_papillomavirus_type = 'HPV')) %>%
  mutate(OS=ifelse(OS=='DECEASED',1,0)) %>%
  mutate(OS.Time=OS.Time/365) %>%
  mutate(stage_group=as.character(factor(stage,levels = c("Stage I","Stage IA","Stage IA1","Stage IA2",
                                                          "Stage IB","Stage IB1","Stage IB2",
                                                          "Stage II", "Stage IIA","Stage IIA1","Stage IIA2",
                                                          "Stage IIB",
                                                          "Stage III", "Stage IIIA", "Stage IIIB", 
                                                          "Stage IVA","Stage IVB"),
                                         labels = c(1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4))))
survival_data <- rbind2(dead_patient,alive_patient) 
names(survival_data)[1] <- 'patient'
#add PFI(progression-free interval)
PFI <- read.table("./data/CESC_survival.txt",sep = "\t",stringsAsFactors = F,header = T,fill = T)
names(PFI)[2] <- 'patient'
PFI <- PFI[,c(2,9,10)]
survival <- merge(survival_data,PFI,by="patient")
survival <- survival %>%
  mutate(PFI.Time=PFI.time/365)

#count 
count <- as.data.frame(fread("./data/gdac.broadinstitute.org_CESC.mRNAseq_Preprocess.Level_3.2016012800.0.0/rawcount.txt",sep = "\t",
                    header = T,fill = T,stringsAsFactors = F))
count <- count[!duplicated(count$`HYBRIDIZATION R`),]
rownames(count) <- count[,1]
count <- count[,-1]
count <- round(count)
tcount <- as.data.frame(t(count))
tcount$group <- as.numeric(as.character(substring(rownames(tcount),14,15)))
tcount$group <- ifelse(tcount$group < 10 , 'T','N')
tumor <- tcount %>%
  dplyr::filter(group == "T")
normal <- tcount %>%
  dplyr::filter(group == "N")
merge <- rbind(tumor,normal)
coldata <- merge %>%
  dplyr::select(group)
coldata <- data.frame(row.names = rownames(coldata),coldata)
tumor <- as.data.frame(t(tumor))
normal <- as.data.frame(t(normal))
a <- cbind(tumor,normal)
a <- a[c(1:(length(rownames(a)) -1)),]
coldata$group <- factor(coldata$group,levels = c("T","N"))
b <- as.matrix(sapply(a, as.numeric))  
row.names(b) <- rownames(a)
dds <- DESeqDataSetFromMatrix(countData =b,colData = coldata,design = ~group)
dds <- DESeq(dds)
res <- results(dds) 
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
normcount <- resdata[,-c(2:7)]
rownames(normcount) <- normcount[,1]
normcount <- normcount[,-1]
HDAC6_expr <- as.data.frame(t(normcount)) %>%
  dplyr::select(HDAC6)
HDAC6_expr$group <- as.numeric(as.character(substring(rownames(HDAC6_expr),14,15)))
HDAC6_expr$group <- ifelse(HDAC6_expr$group < 10 , 'T','N')
HDAC6_expr$log_HDAC6 <- log(HDAC6_expr$HDAC6)
p <- ggboxplot(HDAC6_expr,x='group',y='log_HDAC6',fill = 'group',bxp.errorbar = T,bxp.errorbar.width = 0.2)
p+stat_compare_means()


#correlation plot (TP53-HDAC6)
tnorm <- as.data.frame(t(normcount))
tnorm$group <- as.numeric(as.character(substring(rownames(tnorm),14,15)))
tnorm$group <- ifelse(tnorm$group < 10 , 'T','N')
tnorm <- tnorm %>%
  dplyr::filter(group == "T")
tnorm <- tnorm[,-20502]
rownames(tnorm) <- make.names(str_sub(rownames(tnorm),1,12),unique = T)
rownames(tnorm) <- gsub('.','-',rownames(tnorm),fixed = T)
ggscatter(tnorm, x = "HDAC8", y = "HDAC6", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "HDAC8", ylab = "HDAC6")
cor.test(tnorm$HDAC8, tnorm$HDAC6, 
                method = "pearson")

#OS survival 
HDAC6 <- tnorm %>%
  dplyr::select(HDAC6)
HDAC6[,2] <- rownames(HDAC6)
names(HDAC6)[2] <- 'patient'
rownames(HDAC6) <- NULL
os <- merge(HDAC6,survival,by="patient")
os <- unique(os)
os <- os %>%
  mutate(zscore = (HDAC6-mean(HDAC6))/sd(HDAC6))
os$log_HDAC6 <- log(os$HDAC6 + 1) #log(norm +1)
os$OS.Time <- ifelse(os$OS.Time > 5, 5, os$OS.Time)
os$PFI.Time <- ifelse(os$PFI.Time > 5, 5, os$PFI.Time)
os$zscore_group <- ifelse(os$zscore > 2,"High","Low" )


prob=c(0.25,0.75)
upper=quantile(os$log_HDAC6,prob[2],names = FALSE)
lower=quantile(os$log_HDAC6,prob[1],names = FALSE)
os$quan_group<- ifelse(os$log_HDAC6>upper,"High",
                     ifelse(os$log_HDAC6>lower,"middle","Low"))

os$HPV <- ifelse(is.na(os$HPV),0,os$HPV) #NA>0 
os$HPV_group <- ifelse(os$HPV == 0,"NO",
                       ifelse(os$HPV == "Not applicable","NO","YES"))

save(tnorm,file = "./data/tnorm.RData")
save(os,file = "./data/os.RData")


#quantile survival analysis 
os_quan <- os %>%
  dplyr::filter(quan_group == "High" | quan_group == "Low")

require("survival")
sfit <- survfit(Surv(OS.Time,OS)~quan_group,data = os_quan)
sfit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = os_quan)
ggsurvplot(sfit,pval = TRUE)
ggsurvplot(
  stage1fit,
  data = stage1,
  size = 1,                 # change line size
  palette =
    c("#FF0000", "#0000FF"),# custom color palettes
  conf.int = F,  # Add confidence interval
  conf.int.style = "step",
  pval = TRUE, # Add p-value
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  legend.title = "quartile",
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
) +
  ggtitle("CESC HPV(-) stageI-HDAC6")

ggsurvplot(
  stage4fit,
  data = stage4,
  size = 1,       # change line size
  conf.int = F,  # Add confidence interval
  pval = TRUE, # Add p-value
  font.x = c(10), 
  font.y = c(10),
  font.tickslab = c(10),
  break.time.by = 1,
  xlim = c(0,5),
  xlab = "Time in years",
  ylab = "Overall survival",
  surv.scale="percent",
  palette = c("#FF0000", "#0000FF"),
  legend  = c(0.15,0.4),
  legend.title = "zscore",
  legend.labs = c("High", "Low"),
  conf.int.style = "step",
  risk.table = T,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.20, # Useful to change when you have multiple groups
  risk.table.title="Numbers of risk",
  ggtheme = theme_bw(),# Change ggplot2 theme
  tables.theme = theme_cleantable()
) +
  ggtitle("CESC stageIII-HDAC6")


ggsurvplot(
  stage1fit,
  data = stage1,
  size = 1,       # change line size
  conf.int = F,  # Add confidence interval
  pval = TRUE, # Add p-value
  font.x = c(10), 
  font.y = c(10),
  font.tickslab = c(10),
  break.time.by = 1,
  xlim = c(0,5),
  xlab = "Time in years",
  ylab = "Progression-free survival",
  surv.scale="percent",
  palette = c("#FC4E07","#2E9FDF"),
  legend  = c(0.15,0.4),
  legend.title = "zscore",
  legend.labs = c("High", "Low"),
  conf.int.style = "step",
  risk.table = T,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.height = 0.20, # Useful to change when you have multiple groups
  risk.table.title="Numbers of risk",
  ggtheme = theme_bw(),# Change ggplot2 theme
  tables.theme = theme_cleantable()
) +
  ggtitle("CESC stageI-HDAC6")

                       
                      
#stage -OS
fit <- coxph(Surv(OS.Time,OS)~stage_group+HDAC6,data = os_quan)
stage <- survfit(Surv(OS.Time,OS)~quan_group+strata(os_quan$stage_group),data = os_quan)
ggsurvplot(stage,pval = TRUE)
stage1 <- os_quan %>% dplyr::filter(stage_group == 1)
stage2 <- os_quan %>% dplyr::filter(stage_group == 2)
stage3 <- os_quan %>% dplyr::filter(stage_group == 3)
stage4 <- os_quan %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage1)
stage2fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage2)
stage3fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage3)
stage4fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage4)
ggsurvplot(stage1fit,pval = TRUE,risk.table = TRUE)
ggsurvplot(
  stage3fit,
  data = stage3,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = T,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
)

#stage - PFI
fit <- coxph(Surv(PFI.Time,PFI)~stage_group+HDAC6,data = os_quan)
stage <- survfit(Surv(PFI.Time,PFI)~quan_group+strata(os_quan$stage_group),data = os_quan)
stage1 <- os_quan %>% dplyr::filter(stage_group == 1)
stage2 <- os_quan %>% dplyr::filter(stage_group == 2)
stage3 <- os_quan %>% dplyr::filter(stage_group == 3)
stage4 <- os_quan %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage1)
stage2fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage2)
stage3fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage3)
stage4fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage4)


##zscore -survival 
require("survival")
sfit <- survfit(Surv(OS.Time,OS)~zscore_group,data = os)
sfit <- survfit(Surv(PFI.Time,PFI)~zscore_group,data = os)
ggsurvplot(
  sfit,
  data = os_zscore,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
)

#stage -OS
fit <- coxph(Surv(OS.Time,OS)~stage_group+HDAC6,data = os)
stage <- survfit(Surv(OS.Time,OS)~zscore_group+strata(os$stage_group),data = os)
ggsurvplot(stage,pval = TRUE)
stage1 <- os %>% dplyr::filter(stage_group == 1)
stage2 <- os %>% dplyr::filter(stage_group == 2)
stage3 <- os %>% dplyr::filter(stage_group == 3)
stage4 <- os %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(OS.Time,OS)~zscore_group,data = stage1)
stage2fit <- survfit(Surv(OS.Time,OS)~zscore_group,data = stage2)
stage3fit <- survfit(Surv(OS.Time,OS)~zscore_group,data = stage3)
stage4fit <- survfit(Surv(OS.Time,OS)~zscore_group,data = stage4)
ggsurvplot(stage1fit,pval = TRUE,risk.table = TRUE)

#stage -PFI

fit <- coxph(Surv(PFI.Time,PFI)~stage_group+HDAC6,data = os_zscore)
stage <- survfit(Surv(PFI.Time,PFI)~zscore_group+strata(os_zscore$stage_group),data = os_zscore)
ggsurvplot(stage,pval = TRUE)
stage1 <- os_zscore %>% dplyr::filter(stage_group == 1)
stage2 <- os_zscore %>% dplyr::filter(stage_group == 2)
stage3 <- os_zscore %>% dplyr::filter(stage_group == 3)
stage4 <- os_zscore %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(PFI.Time,PFI)~zscore_group,data = stage1)
stage2fit <- survfit(Surv(PFI.Time,PFI)~zscore_group,data = stage2)
stage3fit <- survfit(Surv(PFI.Time,PFI)~zscore_group,data = stage3)
stage4fit <- survfit(Surv(PFI.Time,PFI)~zscore_group,data = stage4)
ggsurvplot(stage1fit,pval = TRUE,risk.table = TRUE)
ggsurvplot(
  stage1fit,
  data = stage1,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
)

#HPV- quartile  
os_HPV <- os %>%
  dplyr::filter(HPV_group == "YES")
no_HPV <- os %>%
  dplyr::filter(HPV_group == "NO")
quan_HPV <- os_HPV %>%
  dplyr::filter(quan_group == "High" | quan_group == "Low")
quan_no <- no_HPV %>%
  dplyr::filter(quan_group == "High" | quan_group == "Low")

require("survival")
sfit <- survfit(Surv(OS.Time,OS)~quan_group,data = quan_no)
sfit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = quan_no)
ggsurvplot(sfit,pval = TRUE)
ggsurvplot(
  sfit,
  data = quan_HPV,
  size = 1,                 # change line size
  palette =
    c("#FF0000", "#0000FF"),# custom color palettes
  conf.int = F,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
)

#stage -OS
fit <- coxph(Surv(OS.Time,OS)~stage_group+HDAC6,data = os_quan)
stage <- survfit(Surv(OS.Time,OS)~quan_group+strata(quan_no$stage_group),data = quan_no)
ggsurvplot(stage,pval = TRUE)
stage1 <- quan_no %>% dplyr::filter(stage_group == 1)
stage2 <- quan_no %>% dplyr::filter(stage_group == 2)
stage3 <- quan_no %>% dplyr::filter(stage_group == 3)
stage4 <- quan_no %>% dplyr::filter(stage_group == 4)
stage1fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage1)
stage2fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage2)
stage3fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage3)
stage4fit <- survfit(Surv(OS.Time,OS)~quan_group,data = stage4)
ggsurvplot(
  stage1fit,
  data = stage1,
  size = 1,                 # change line size
  palette =
    c("#FF0000", "#0000FF"),# custom color palettes
  conf.int = F,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
)

#stage - PFI
fit <- coxph(Surv(PFI.Time,PFI)~stage_group+HDAC6,data = os_quan)
stage <- survfit(Surv(PFI.Time,PFI)~quan_group+strata(quan_no$stage_group),data = quan_no)

stage1fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage1)
stage2fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage2)
stage3fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage3)
stage4fit <- survfit(Surv(PFI.Time,PFI)~quan_group,data = stage4)


#except stage3 
os_124 <- os %>%
  dplyr::filter(stage_group != 3) %>%
  dplyr::filter(quan_group == "High" | quan_group == "Low")
os_12 <- os %>%
  dplyr::filter(stage_group == 1 | stage_group == 2) %>%
  dplyr::filter(quan_group == "High" | quan_group == "Low")
require("survival")
sfit <- survfit(Surv(OS.Time,OS)~zscore_group,data = os_124)
sfit <- survfit(Surv(PFI.Time,PFI)~zscore_group,data = os_124)
ggsurvplot(sfit,pval = TRUE)
ggsurvplot(
  sfit,
  data = os_124,
  size = 1,                 # change line size
  palette =
    c("#FF0000", "#0000FF"),# custom color palettes
  conf.int = F,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("High", "Low"),    # Change legend labels
  risk.table.height = 0.18, # Useful to change when you have multiple groups
  ggtheme = theme_classic()      # Change ggplot2 theme
)

