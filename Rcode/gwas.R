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
