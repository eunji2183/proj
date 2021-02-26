### 2021/02/26 ###
dir <- c("./multi_PRS")
filelist <- str_subset(list.files(dir),"txt")
filenames <- dir(dir,recursive = T,pattern = ".txt",full.names = T)
data <- lapply(filenames,function(x){
  read.table(x,header = T,sep = " ",stringsAsFactors = F,fill = T,as.is = T,
             encoding="UTF-8",fileEncoding = "CP949" )})
count <- Reduce(function(x,y)merge(x,y,by="DIST_ID",all=TRUE),data)
rownames(count) <- count[,1]
count <- count[,-1]
head <- read.table("./bfile/con",header = F)
names(count) <- head$V1
count <- data.frame(DIST_ID=rownames(count),count)
rownames(count) <- NULL

catPRS <- read.table("./cat_prsice/merge.txt",sep = "\t",fill = T,header = T)
contPRS <- read.table("./cont_prsice/merge.txt",sep = "\t",fill = T,header = T)

merge <- merge(catPRS,contPRS,by="DIST_ID",all = TRUE)
merge2 <- merge(merge,count,by="DIST_ID",all=TRUE)


write.csv(merge2,file = "./merge.csv",row.names = T)
write.csv(count,file = "./multi.PRS.csv",row.names = T)
