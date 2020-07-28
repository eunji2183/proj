# ovarian serous cystadenocarcinoma , DNA methylation 
#https://pubmed.ncbi.nlm.nih.gov/30446011/
#univariate cox - multivariate cox - survival- ROC 
rm(list = ls())
working_dir="/home/eunji/R/project/200612_OSC/"
setwd(working_dir)
options(stringsAsFactors = F)
save(list = ls(),file="./data/OSC.RData")
BiocManager::install("pROC")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("missMethyl")
BiocManager::install("minfiData")
BiocManager::install("Gviz")
BiocManager::install("DMRcate")
install.packages("table1")
install.packages("caret")
BiocManager::install("maftools")
library(reshape) #rename function
library(UCSCXenaTools) #XenaPrepare, getTCGAdata
library(magrittr)  # %>% function 
library(dplyr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(pROC)
library(table1)
library(tidyr)
library(caret)
library(maftools)
library(boot) #example
library(survival)
library(survminer)
library(dplyr)
library(DT)
#########################################################################################
## TCGA UCSC xena : clinical data 
OV <- getTCGAdata(project = 'OV',download = T,Methylation = T,MethylationType = c("27K")) 
OV=XenaPrepare(OV) 
meth <- as.data.frame(OV$HumanMethylation27.gz)
clin <- as.data.frame(OV$OV_clinicalMatrix)
tumor <- clin[grep("0.$",clin$sampleID),]
tumor <- tumor %>%
  dplyr::filter(histological_type == "Serous Cystadenocarcinoma")%>%
  dplyr::select(sampleID,age_at_initial_pathologic_diagnosis,
                additional_pharmaceutical_therapy,
                additional_radiation_therapy,days_to_death,
                days_to_last_followup,gender,vital_status,clinical_stage,
                neoplasm_histologic_grade,tumor_residual_disease,
                anatomic_neoplasm_subdivision)%>%
  dplyr::filter(!is.na(days_to_last_followup) & !is.na(vital_status))
tumor <- tumor[-grep("G[1BX]| ",tumor$neoplasm_histologic_grade),] 
tumor$clinical_stage[is.na(tumor$clinical_stage)] <- "Unknown"
tumor$neoplasm_histologic_grade[is.na(tumor$neoplasm_histologic_grade)] <- "Others"
tumor$tumor_residual_disease[is.na(tumor$tumor_residual_disease)] <- "Unknown"
tumor$anatomic_neoplasm_subdivision[is.na(tumor$anatomic_neoplasm_subdivision)] <- "Unknown"
tumor <- tumor %>%
  dplyr::mutate(age=ifelse(age_at_initial_pathologic_diagnosis< 60 ,"<60",">=60")) %>%
  dplyr::mutate(stage=as.character(factor(clinical_stage,levels = c("Stage IA","Stage IB","Stage IC","Stage IIA","Stage IIB", "Stage IIC",
                                                                    "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                    "Stage IV","Unknown"),
                                          labels = c(1,1,1,2,2,2,3,3,3,4,0)))) %>%
    dplyr::mutate(fustat=ifelse(vital_status=='DECEASED',1,0)) %>%
    dplyr::mutate(time=days_to_last_followup/30) 
#########################################################################################################################
#methylation data (27K)
meth <- as.data.frame(OV$HumanMethylation27.gz)
rownames(meth) <- meth[,1]
meth <- meth[,-1]
meth <- na.omit(meth)
tmeth <- as.data.frame(t(meth))
meth <- data.frame(sampleID=rownames(tmeth),tmeth)
rownames(meth) <- NULL
merge <- merge(tumor,meth,by="sampleID")
######################################################################################
#training set vs validation set (2:1)
index <- createDataPartition(y=merge$age,p=0.667,list = F)
train <- merge[index,]
valid <- merge[-index,]
train <- train %>%
  dplyr::mutate(dataset= 1)
valid <- valid %>%
  dplyr::mutate(dataset= 2)
data <- rbind(train,valid)

################################################################################################
#table 
data$dataset <- 
  factor(data$dataset,
         levels = c(1,2),
         labels = c("Training dataset","Validation dataset"))
data$stage <- 
  factor(data$stage,
         levels = c(1,2,3,4,0),
         labels = c("I","II","III","IV","Unknown"))
data$tumor_residual_disease <- 
  factor(data$tumor_residual_disease,
         levels = c("No Macroscopic disease","1-10 mm","11-20 mm",">20 mm","Unknown"))
data$neoplasm_histologic_grade <- 
  factor(data$neoplasm_histologic_grade,
         levels = c("G2","G3","G4","Others"))
table1(~age + stage + neoplasm_histologic_grade +
        tumor_residual_disease +
         anatomic_neoplasm_subdivision | dataset,data=data)
####################################################################################
#univariate cox proportional hazard regression analysis 
options(scipen = 999) # exponential to number 
clin <- data[,c(1:16,21683)]
meth <- data[,c(1,17:21682)]
rownames(meth) <- meth[,1]
meth <- meth[,-1]
covariates <- c(colnames(meth))
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, fustat)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
                         res<-c(beta, HR, CI , wald.test, p.value)
                         names(res)<-c("beta", "HR", "CI 95%", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                         })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res,stringAsFactors=F)
res <- subset(res, p.value < 0.05 & abs(as.numeric(beta))< 5)
res <- res[order(res$p.value),]
save(res,file = "./result/uni_res.RData")

##########################################################################
# meth sites  survival (>1)
cpg=c('cg04907664','cg13652336','cg25123470') 
splots <- lapply(cpg, function(g){
  data$gene=ifelse(data[,g]>median(data[,g]),'high','low')
  sfit1=survfit(Surv(time, fustat)~gene, data=data)
  ggsurvplot(sfit1,pval =TRUE, data = data, risk.table = TRUE)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 3, nrow = 1, risk.table.height = 0.2)

######################################################################################
#multivariate cox hazard regression analysis 
#risk score fomula 
cpg=c(rownames(res))
train_clin <- train[,1:16]
train_meth <- data.frame(train_clin,train[,c(cpg)])
fmla <- as.formula(paste("Surv(time, fustat) ~", paste0(cpg, collapse = "+")))
cox <- coxph(fmla,data = train_meth)
coef <- cox$coefficients
rownames(train_meth) <- train_meth[,1]
train_meth <- train_meth[,-1]
train_meth2 <- data.frame(train_meth[,c(cpg)])
riskscore <- as.matrix(train_meth2) %*% coef
risk_data <- data.frame(train_clin,train_meth2,riskscore=riskscore)
rownames(risk_data) <- NULL

valid_clin <- valid[,1:16]
valid_meth <- data.frame(valid_clin,valid[,c(cpg)])
fmla <- as.formula(paste("Surv(time, fustat) ~", paste0(cpg, collapse = "+")))
cox <- coxph(fmla,data = valid_meth)
coef <- cox$coefficients
rownames(valid_meth) <- valid_meth[,1]
valid_meth <- valid_meth[,-1]
valid_meth2 <- data.frame(valid_meth[,c(cpg)])
riskscore1 <- as.matrix(valid_meth2) %*% coef
risk_data1 <- data.frame(valid_clin,valid_meth2,riskscore=riskscore1)
rownames(risk_data1) <- NULL

#######################################################################################
# risk score survival (training & validation)
#training dataset
risk_data$risk_group <- ifelse(risk_data$riskscore > median(risk_data$riskscore),'high','low')
sfit1=survfit(Surv(time, fustat)~risk_group, data=risk_data)
ggsurvplot(sfit1,pval =TRUE, data = risk_data, risk.table = TRUE)

#validation  dataset
risk_data1$risk_group <- ifelse(risk_data1$riskscore > median(risk_data1$riskscore),'high','low')
sfit1=survfit(Surv(time, fustat)~risk_group, data=risk_data1)
ggsurvplot(sfit1,pval =TRUE, data = risk_data1, risk.table = TRUE)


################################################################################################
#ROC analysis 
#training dataset 
train_roc <- roc(risk_data$fustat,risk_data$riskscore,ci=TRUE)
plot.roc(train_roc,print.auc = TRUE,main="training-ROC")
valid_roc <- roc(risk_data1$fustat,risk_data1$riskscore,ci=TRUE)
plot.roc(valid_roc,print.auc = TRUE,main="validation ROC")
roc.test(train_roc,valid_roc)

uni_cox_in_bulk <- function(meth_list,survival_info_df){
  library(survival)
  uni_cox <- function(single_meth){
    formula <- as.formula(paste0('Surv(days_to_last_followup,fustat)~',single_meth))
    surv_uni_cox <- summary(coxph(formula,data = survival_data,method = ))
    ph_hypothesis_p <- cox.zph(coxph(formula,data = survival_data))$table[1,3]
    if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){
      single_cox_report <- data.frame('uni_cox_sig_meth'= single_meth,
                                      'beta'=surv_uni_cox$coefficients[,1],
                                      'Hazard_Ratio'=exp(surv_uni_cox$coefficients[,1]),
                                      'z_pvalue'=surv_uni_cox$coefficients[,5],
                                      'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
                                      'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
      single_cox_report
    }}
  uni_cox_list <- lapply(meth_list,uni_cox)
  do.call(rbind,uni_cox_list)
}
uni_cox_df <- uni_cox_in_bulk(meth_list = meth_list, survival_info_df = survival_data)

