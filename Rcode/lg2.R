#chemotherapy 
save(list = ls(),file="/HDD8T/eunji/proj/lg/lg.RData")
setwd("/home/eunji/R/tmp/lg/data/")

cancer  <- TCGAbiolinks:::getGDCprojects()$project_id
cancer <- str_subset(cancer, "TCGA") 
cancer_type <- sort(cancer) #33  types&proj
cancer <- strsplit(cancer_type,split = "-")
for(i in 1:length(cancer)){
  cancer[i] <- cancer[[i]][2]
}
cancer <- unlist(cancer)

#33 cancer - drug data (LAML remove, no data)
i=14
drugtcga <- list()
for(i in 15:33){
  query <- GDCquery(project = cancer_type[i], 
                    data.category = "Clinical",
                    data.type = "Clinical Supplement", 
                    data.format = "BCR Biotab")
  GDCdownload(query)
  clinical.BCRtab.all <- GDCprepare(query)
  drugtcga[[i]] <- as.data.frame(clinical.BCRtab.all[[as.character(str_subset(names(clinical.BCRtab.all), "drug"))]])
  names(drugtcga)[[i]] <- cancer[i]
}

drugtcga[[14]] <- NULL
drugtcga[[7]] <- NULL
save(drugtcga,file = "/HDD8T/eunji/proj/lg/drugtcga.RData")


for(i in 1:31){
  drugtcga[[i]] <- drugtcga[[i]] %>%
    dplyr::select(bcr_patient_barcode,pharmaceutical_therapy_drug_name,
                  #clinical_trial_drug_classification,
                  pharmaceutical_therapy_type,
                  pharmaceutical_tx_started_days_to,pharmaceutical_tx_ended_days_to,
                  pharmaceutical_tx_ongoing_indicator,treatment_best_response)
  drugtcga[[i]]$cancer <- names(drugtcga)[[i]]
}
drugdata <- Reduce(function(x,y)rbind(x,y),drugtcga)

drugdata <- drugdata %>%
  dplyr::filter(!bcr_patient_barcode %in% c("bcr_patient_barcode","CDE_ID:2003301"))
drugdata2 <- drugdata %>%
  #dplyr::filter(pharmaceutical_therapy_type == "Chemotherapy") %>%
  dplyr::filter(treatment_best_response %in% c("Stable Disease",
                                               "Partial Response",
                                               "Complete Response",
                                               "Clinical Progressive Disease")) %>%
  dplyr::filter(!pharmaceutical_therapy_drug_name %in% c("[Not Available]",
                                                         "[Unknown]"))

drugdata$pharmaceutical_therapy_drug_name <-toupper(drugdata$pharmaceutical_therapy_drug_name) 

save(drugdata,file = "/HDD8T/eunji/proj/lg/drugdata.RData")

drugdata2$pharmaceutical_therapy_drug_name <-toupper(drugdata2$pharmaceutical_therapy_drug_name) 
drugdata2$pharmaceutical_therapy_drug_name <- gsub(' ','-',drugdata2$pharmaceutical_therapy_drug_name,fixed = T)
drugdata2$pharmaceutical_therapy_drug_name <- gsub('/','-',drugdata2$pharmaceutical_therapy_drug_name,fixed = T)
drugdata2$pharmaceutical_therapy_drug_name <- gsub('_','-',drugdata2$pharmaceutical_therapy_drug_name,fixed = T)


drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "5-FLU")] = '5-FU'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "5-FL")] = '5-FU'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "5FU")] = '5-FU'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "XELODA")] = 'XELODA'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "VINORELBIN")] = 'VINORELBINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "VINCRISTINE")] = 'VINCRISTINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TOPOTECAN")] = 'TOPOTECAN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TEMODAL")] = 'TEMOZOLOMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TEMODAR")] = 'TEMOZOLOMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DOXETAXOL")] = 'DOCETAXEL'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TAXOL")] = 'PACLITAXEL'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TAXOTERE")] = 'DOCETAXEL'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ADRIAMY")] = 'ADRIAMYCIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DOXO")] = 'DOXORUBICIN'

drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ADRIAMYCIN")] = 'DOXORUBICIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DOCETAXELUM")] = 'DOCETAXEL'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DOXETAXEL")] = 'DOCETAXEL'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ALIM")] = 'PEMETREXED'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "PEMETREXED")] = 'PEMETREXED'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ALMITA")] = 'PEMETREXED'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "NAVELBINE")] = 'VINORELBINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "AFINITOR")] = 'EVEROLIMUS'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "AVASTIN")] = 'BEVACIZUMAB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "AXITNIB")] = 'AXITINIB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CARMUSTIN")] = 'CARMUSTINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "BCNU")] = 'CARMUSTINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CAMPTOSAR")] = 'IRINOTECAN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CAPECYTABINUM")] = 'CAPECITABINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CARBOPLATINUM")] = 'CARBOPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "PARAPLATIN")] = 'CARBOPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "FLUOROURACIL")] = '5-FU'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CARBOPLATIN-5-AUC")] = 'CARBOPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "LOMUSTIN")] = 'LOMUSTIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CCNU")] = 'LOMUSTIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CDDP")] = 'CISPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ERBITUX")] = 'CETUXIMAB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "FLUORURACIL")] = '5-FU'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CISPLAT")] = 'CISPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "COPOLANG-CAPS")] = 'COPOLANG'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CPT-11")] = 'IRINOTECAN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CYCLOPHOSPHA")] = 'CYCLOPHOSPHAMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "CYTOXAN")] = 'CYCLOPHOSPHAMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DACABARZINE")] = 'DACARBAZINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DACARBAZINE")] = 'DACARBAZINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "DTIC")] = 'DACARBAZINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ELOXATIN")] = 'OXALIPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "EPIRUBIC")] = 'EPIRUBICIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TARCEVA")] = 'ERLOTINIB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ET-743")] = 'TRABECTEDIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "VEPESID")] = 'ETOPOSIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "ETOP")] = 'ETOPOSIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "PLATINOL")] = 'CISPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "GEMCITABINE")] = 'GEMCITABINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "GEMZAR")] = 'GEMCITABINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "FOLINIC-ACID")] = 'LEUCOVORIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "GLEEVEC")] = 'IMATINIB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "GLIADEL")] = 'CARMUSTINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "HALICHONDRIN")] = 'HALICHONDRIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "HERCEPTIN")] = 'TRASTUZUMAB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "IFEX")] = 'IFOSFAMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "IFOSFAMID")] = 'IFOSFAMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "IRINOTECAN-HCL")] = 'IRINOTECAN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "IRINOTECAN-HYDROCHLORIDE")] = 'IRINOTECAN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "JEVTANA")] = 'CABAZITAXEL'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "LEUCOVORIN-CALCIUM")] = 'LEUCOVORIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "MITOMYCIN")] = 'MITOMYCIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "NAXAVAR")] = 'SORAFENIB '
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "NEXAVAR")] = 'SORAFENIB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "OXALIPLATINUM")] = 'OXALIPLATIN'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "SAHA-VORINOSTAT")] = 'VORINOSTAT'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "SORAFENIB")] = 'SORAFENIB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "SUTENT")] = 'SUNITINIB'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "TEMOZOLAMIDE")] = 'TEMOZOLOMIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "VP-16")] = 'ETOPOSIDE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "XELODA")] = 'CAPECITABINE'
drugdata2$pharmaceutical_therapy_drug_name[str_which(drugdata2$pharmaceutical_therapy_drug_name, "YONDELIS")] = 'TRABECTEDIN'
drugname <- as.data.frame(table(drugdata2$pharmaceutical_therapy_drug_name))
drugname2 <- drugname %>%
  dplyr::filter(Freq > 10)
drugname2 <-drugname2[order(drugname2$Freq,decreasing = T),]
drugname3 <- as.character(drugname2$Var1)
drugdata3 <- drugdata2 %>% dplyr::filter(pharmaceutical_therapy_drug_name %in% c(drugname3))
drugdata3$pharmaceutical_therapy_drug_name[str_which(drugdata3$pharmaceutical_therapy_drug_name, "DOXIL")] = 'DOXORUBICIN'

drugsample <- as.data.frame(table(drugdata3$bcr_patient_barcode))

#monotherapy 
a <- drugsample %>%
  dplyr::filter(bcr_patient_barcode %in% c("TCGA-AF-A56N","TCGA-L1-A7W4","TCGA-E1-A7Z2","TCGA-DU-A7TB","TCGA-06-1806","TCGA-19-A6J5","TCGA-DU-A6S6",
                                           "TCGA-DU-A6S7","TCGA-DU-A7TA","TCGA-E1-A7YS","TCGA-E1-A7Z6","TCGA-FG-7638",
                                           "TCGA-FG-8185","TCGA-FG-8189","TCGA-FG-8191","TCGA-FG-A4MW","TCGA-HT-8011",
                                           "TCGA-HT-8113","TCGA-HT-A4DV","TCGA-HW-A5KL","TCGA-S9-A89V") )
drugsample2 <- drugsample %>%
  dplyr::filter(Freq ==1)
drugsample2 <- rbind(a,drugsample2)
names(drugsample2)[1] <-names(drugdata)[1]
drugsample2 <- merge(drugsample2,drugdata3,by="bcr_patient_barcode")
drugsample2$response <- ifelse(drugsample2$treatment_best_response == "Clinical Progressive Disease","NO","YES")
b <- drugsample2[,c(1,2,3,9,10)]
b <- unique(b)
b$N <- rownames(b)
b <- b %>% dplyr::filter(!N %in% c(203,191,2,295,205))
c <- dcast(b, pharmaceutical_therapy_drug_name~cancer)
c <- data.frame(row.names = c$pharmaceutical_therapy_drug_name,c[,2:28])

#cisplatin-CESC(57),HNSC(43), capecitabine-COAD(10),ESCA(9),
#TEMOZOLOMIDE-GBM(8),LGG(96), SORAFENIB-LIHC(9),GEMCITABINE-PAAD(49)
#5-FU-READ(10),STAD(44), DACARBAZINE-SKCM(12)
mono <- readxl::read_excel("/home/eunji/refpaper.xlsx",
                         sheet = "Sheet6",
                         range = "B8:D18", 
                         col_names = F,
                         na="NA")
names(mono) <- c("drug","cancer","Total")

monolist <- list()
for(i in 1:nrow(mono)){
  monolist[[i]] <- b %>% dplyr::filter(cancer == as.character(mono[i,2])
                                       & pharmaceutical_therapy_drug_name==
                                         as.character(mono[i,1]))
  names(monolist)[[i]] <- as.character(mono[i,2])}

for(i in 1:11){
  mono[i,4] <- as.data.frame(table(monolist[[i]][["response"]]))[1,2]
  mono[i,5] <- as.data.frame(table(monolist[[i]][["response"]]))[2,2]
}
names(mono)[c(4,5)] <- c('NO','YES')
mono[9,4] <- NA
mono[9,5] <- 10

mono <- mono %>% dplyr::filter(NO > 2 & YES >2)
monode <- list()
i="CESC"
tcount2 <- data.frame(Barcode=tcount$Barcode,tcount[,303:55611])
for( i in c(mono$cancer)) {
  a<- tcount2 %>%  dplyr::filter (Barcode %in% c(as.character(monolist[[i]]$bcr_patient_barcode)))
  d <- data.frame(Barcode = as.character(monolist[[i]]$bcr_patient_barcode),chemotherapy=monolist[[i]]$response)
  d <- merge(d,a,by="Barcode")
  d <- d[order(d$chemotherapy),]
  monode[[i]][["count"]] <- data.frame(row.names=d$Barcode,d[,-c(1:2)])
  monode[[i]][["anno"]] <- data.frame(row.names=d$Barcode,chemotherapy=d[,2])
}

monode[["GBM"]] <- NULL

for(i in 1:4){
  monode[[i]]$anno$chemotherapy <- factor(monode[[i]]$anno$chemotherapy,levels = c("NO","YES"))
  a <- monode[[i]]$count
  a <- as.data.frame(t(a))
  b <- as.matrix(sapply(a,as.numeric))
  b[is.na(b)] <- 0
  row.names(b) <- rownames(a)
  dds <- DESeqDataSetFromMatrix(countData = b,colData = monode[[i]]$anno,design = ~chemotherapy)
  dds <- DESeq(dds)
  monode[[i]]$dds <- dds
  res <-  DESeq2::results(dds)
  res <- as.data.frame(res)
  res <- res[order(res$padj),]
  res <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
  monode[[i]]$normcount <- as.data.frame(counts(dds,normalized=TRUE))
  monode[[i]]$res <- res}

for(i in 1:5){
  monode[[i]]$rlog <- rlog(monode[[i]]$dds)}

FA <- readxl::read_excel("/HDD8T/eunji/proj/lg/FA.xlsx",
                         sheet = "Sheet1",
                         range = "A1:B21", 
                         col_names = F,
                         na="NA")
names(FA) <- c('gene_name','PATHWAY')

glu <- 

drugsample3 <- dcast(drugdata3, bcr_patient_barcode ~pharmaceutical_therapy_drug_name)
names(drugsample)[1] <-names(drugdata)[1]
drugsample4 <- merge(drugsample,drugsample3,by="bcr_patient_barcode")
drugsample5 <- dcast(drugsample2, pharmaceutical_therapy_drug_name~cancer)

drugsample4$bcr_patient_barcode <- as.character(drugsample4$bcr_patient_barcode)
i=371
h=6
for(i in 1:nrow(drugsample4)){
  for(h in 3:34){
    a <- ifelse(drugsample4[i,h] == drugsample4[i,2] ,drugsample4[i,1],NA)
  }
  drugsample4[i,35] <- a
}


#BRCA-tamoxifen 

tamoxifen <- drugdata[str_which(drugdata$pharmaceutical_therapy_drug_name, "TAMOXIFEN"),]
tamoxifen <- drugdata %>%
  dplyr::filter(treatment_best_response %in% c("Stable Disease",
                                               "Partial Response",
                                               "Complete Response",
                                               "Clinical Progressive Disease")) %>%
  dplyr::filter(cancer == "BRCA")
tamoxifen$response <- ifelse(tamoxifen$treatment_best_response == "Clinical Progressive Disease","NO","YES")
tamoxifen <- tamoxifen[order(tamoxifen$treatment_best_response),]
tamoxifen <- tamoxifen[!duplicated(tamoxifen$bcr_patient_barcode),]
a<- tcount2 %>%  dplyr::filter (Barcode %in% c(as.character(tamoxifen$bcr_patient_barcode)))
d <- data.frame(Barcode = as.character(tamoxifen$bcr_patient_barcode),therapy=tamoxifen$response)
d <- merge(d,a,by="Barcode")
d <- d[order(d$therapy),]
monode[["BRCA"]][["count"]] <- data.frame(row.names=d$Barcode,d[,-c(1:2)])
monode[["BRCA"]][["anno"]] <- data.frame(row.names=d$Barcode,therapy=d[,2])
monode[["BRCA"]][["anno"]]$therapy <- factor(monode[["BRCA"]][["anno"]]$therapy,levels = c("NO","YES"))
i="BRCA"
a <- monode[[i]]$count
a <- as.data.frame(t(a))
b <- as.matrix(sapply(a,as.numeric))
b[is.na(b)] <- 0
row.names(b) <- rownames(a)
dds <- DESeqDataSetFromMatrix(countData = b,colData = monode[[i]]$anno,design = ~therapy)
dds <- DESeq(dds)
monode[[i]]$dds <- dds
res <-  DESeq2::results(dds)
res <- as.data.frame(res)
res <- res[order(res$padj),]
res <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
monode[[i]]$normcount <- as.data.frame(counts(dds,normalized=TRUE))
monode[[i]]$res <- res

wilcox.test(before, after, paired = TRUE)

paad <-  paad <- monode[["PAAD"]][["res"]]

paad2 <- paad %>% dplyr::filter(pvalue < 0.05 & log2FoldChange < 0) 
paad2 <- data.frame(gene=paad2$Row.names)

paad3 <- paad %>% dplyr::filter(pvalue < 0.05 & log2FoldChange > 0) 
paad3 <- data.frame(gene=paad3$Row.names)

cesc2 <- data.frame(gene=cesc2$Row.names)
write.table(cesc2,file = "/home/eunji/cesc.txt",row.names = F)
write.table(paad2,file = "/home/eunji/paad2.txt",row.names = F,quote = F)

brca2 <- brca %>% dplyr::filter(pvalue < 0.05 & log2FoldChange < 0) 
brca2 <- data.frame(gene=brca2$Row.names)
write.table(brca2,file = "/home/eunji/brca.txt",row.names = F,quote = F)

stad2 <- STAD %>% dplyr::filter(padj < 0.05 & log2FoldChange < 0) 
stad2 <- data.frame(gene=stad2$Row.names)
write.table(stad2,file = "/home/eunji/stad.txt",row.names = F,quote = F)

lgg2 <- lgg %>% dplyr::filter(pvalue < 0.05 & log2FoldChange < 0) 
lgg2 <- data.frame(gene=lgg2$Row.names)
write.table(lgg2,file = "/home/eunji/lgg.txt",row.names = F,quote = F)


load("/HDD8T/eunji/proj/multi/tgdata.RData")

lgdata <- list()
for ( i in c(names(tgdata))){
  lgdata[[i]]$tpmcount<- tgdata[[i]]$tpmcount 
  lgdata[[i]]$coldata <- tgdata[[i]]$coldata}
rm(tgdata)
i="ACC"
for(i in c(names(lgdata))){
  n <- lgdata[[i]]$coldata %>% dplyr::filter(group=="T")
  lgdata[[i]]$coldata <- n
  lgdata[[i]]$tpmcount <- lgdata[[i]]$tpmcount[,c(rownames(lgdata[[i]]$coldata))]}
save(lgdata,file = "/HDD8T/eunji/proj/lg/lgdata.RData")
rm(tgdata)
h="ACC"
library(Hmisc)
load("/HDD8T/eunji/proj/2D3D/ID.RData")
pc <- ID %>% dplyr::filter(gene_biotype == "protein_coding")
pc <- pc[!duplicated(pc$gene_name),]
pc <- data.frame(row.names = pc$gene_name,gene_id=pc$gene_id)
for(h in 1:9){
  a <- lgdata5[[h]]$tpmcount
  a <- merge(pc,a,by="row.names")
  a <- data.frame(row.names = a$Row.names,a[,-c(1:2)])
  b <- as.data.frame(t(a))
  res <- rcorr(as.matrix(b))
  flattenCorrMatrix <- function(cormat, pmat){
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  res2 <- flattenCorrMatrix(res$r, res$P)
  lgdata5[[h]][['pccor']] <- res2}
lgdata1 <- lgdata[1:6]
lgdata2 <- lgdata[7:12]
lgdata3 <- lgdata[13:18]
lgdata4 <- lgdata[19:24]
lgdata5 <- lgdata[25:33]

save(lgdata5,file = "/HDD8T/eunji/proj/lg/lgdata5.RData")


b <- res2 %>% dplyr::filter(column == "HIF1A")
i=1
lg1 <- list()
for ( i in 1:6){
  res2 <- lgdata1[[i]][['pccor']]
  b <- res2 %>% dplyr::filter(column == "HIF1A")
  b <- b[grep("^ABC", b$row),]
  lg1[[i]] <- b
}

#brca-pam50 
bcpam<- readxl::read_excel("/HDD8T/eunji/proj/lg/brcapam50.xlsx",
                   sheet = "Sheet1",
                   range = "A1:B1070", 
                   col_names = F,
                   na="NA")
names(bcpam) <- c('sampleID','subtype')

bcpam2<- readxl::read_excel("/HDD8T/eunji/proj/lg/brcapam50.xlsx",
                           sheet = "Sheet2",
                           range = "A2:U827", 
                           col_names = T,
                           na="NA")
bcpam2 <- bcpam2[,c(1,21)]
names(bcpam2) <- c('sampleID','subtype')

bcpam3<- readxl::read_excel("/HDD8T/eunji/proj/lg/brcapam50.xlsx",
                            sheet = "Sheet3",
                            range = "A2:C730", 
                            col_names = F,
                            na="NA")
bcpam3 <- bcpam3[,c(1,3)]
names(bcpam3) <- c('sampleID','subtype')

brca <- lgdata1[['BRCA']]$coldata
brca <- data.frame(sampleID=rownames(brca),brca)
brca <- brca %>% dplyr::filter(group=="T")
brca$sampleID <- gsub('.','-', str_sub(rownames(brca),1,12),fixed = T)
brca <- merge(brca,bcpam,by="sampleID",all.x=T)
brca <- brca[!duplicated(brca$sampleID),]
brca <- merge(brca,bcpam2,by="sampleID",all.x=T)
brca <- merge(brca,bcpam3,by="sampleID",all.x=T)

load("/HDD8T/eunji/proj/lg/brcasurv.RData")
names(survival_data)[1] <- 'sampleID'
brca <- merge(brca,survival_data,by="sampleID",all.x=T)
lumA <- brca %>% dplyr::filter(subtype.x=="BRCA.LumA" & subtype.y=="Luminal A") %>%
  dplyr::filter(subtype=="BRCA_LumA")

count <- as.data.frame(t(lgdata[["BRCA"]]$tpmcount))
count$sampleID <- gsub('.','-', str_sub(rownames(count),1,12),fixed = T)
lumA <- merge(lumA,count,by="sampleID")
lumA$hif1a <- ifelse(lumA$HIF1A > median(lumA$HIF1A),"High","Low")
kmfit <-prodlim(Hist(OS.Time,OS) ~ hif1a, data = lumA)
plot(kmfit,percent=FALSE,logrank=T,digits = 4,axes=TRUE,col=c("#800000FF","#374E55FF"),
     axis1.at=seq(0,kmfit$maxtime+1,1),axis1.lab=seq(0,kmfit$maxtime+1,1),
     marktime=F,atrisk=TRUE,xlab="Years",
     confint=TRUE,confint.citype="shadow",#col=c(4,3),
     legend=TRUE,legend.x=2,
     legend.y=0.4,legend.cex=1,
     legend.title="group\n",
     atrisk.labels=paste(c("High","Low"),": "),
     atrisk.title="")

ggscatter(lumA, x = "HIF1A", y = "FABP6",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
