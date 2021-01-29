#pheno file 만들기 
## DIST_ID를 기준으로 csv file merge 

library(stringr)

working_dir="/home/eunji/proj/GWAS/"
setwd(working_dir)
options(stringsAsFactors = F)
dir <- c("./EXCEL/1")
filelist <- str_subset(list.files(dir),"csv")
filenames <- dir(dir,recursive = T,pattern = ".csv",full.names = T)
data <- lapply(filenames,function(x){
  read.table(x,header = T,sep = ",",stringsAsFactors = F,fill = T,as.is = T,
             encoding="UTF-8",fileEncoding = "CP949" )})
count <- Reduce(function(x,y)merge(x,y,by="DIST_ID"),data)

write.csv(count,file = "./EXCEL/1/merge.csv",row.names = F)

pheno <- count %>%
  dplyr::select(DIST_ID,AS1_SEX,AS1_AGE,AS1_FAMNUM,AS1_JOBB,AS1_EDUA,AS1_INCOME,AS1_DRINK,
                AS1_SMOKEA,AS1_PHYSTB,AS1_PHYSIT,AS1_PHYACTL,AS1_PHYACTM,AS1_PHYACTH,
                AS1_MEDCST,AS1_HTFDCST,AS1_PDHT,AS1_PDDM,AS1_PDMI,AS1_PDCH,
                AS1_PDCD,AS1_PDLP,AS1_PDAS,AS1_PDCL,AS1_PDTOTCA1NA,AS1_PDCV,
                AS1_RGMEALFQA,AS1_BEMEALFQA,AS1_SKPMEAL,AS1_EATOUTFQ,
                AS1_WAIST1,AS1_WAIST2,AS1_WAIST3,AS1_HIP1,AS1_HIP2,AS1_HIP3,
                AS1_HEIGHT,AS1_WEIGHT,AS1_BDFTR,AS1_ABFTR,AS1_BMI,AS1_OBDG)
write.csv(pheno,file = "./EXCEL/1/pheno.csv",row.names = F)
                
