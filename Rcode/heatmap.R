#excel-pheatmap 
setwd("/home/eunji/R/project/200717_heatmap")

library(pheatmap)

normalized <- readxl::read_excel("./data/200324 heatmap (intensity log and RNA copy number) (1).xlsx",
                                   sheet = "Sheet1",
                                   range = "B2:F24", 
                                   col_names = T)
normalized <- as.data.frame(normalized)
rownames(normalized) <- normalized[,1]
normalized <- normalized[,-1]
pheatmap(normalized,scale = "row", clustering_distance_row = "correlation",
                 show_colnames = T,show_rownames = T,
                 cluster_cols = T,cluster_rows = T,
                 color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                 border_color = 'white',
                 display_numbers = F,main = "",key=T)

raw <- readxl::read_excel("./data/200324 heatmap (intensity log and RNA copy number) (1).xlsx",
                                 sheet = "Sheet1",
                                 range = "B2:j24", 
                                 col_names = T)
raw <- as.data.frame(raw)
rownames(raw) <- raw[,1]
raw <- raw[,-1]
raw <- raw[,5:8]
colnames(raw) <- c("control1","control2","UPM1","UPM2")
pheatmap(raw,scale = "row", clustering_distance_row = "correlation",
         show_colnames = T,show_rownames = T,
         cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         border_color = 'white',
         display_numbers = F,main = "",key=T)
save(list = ls(),file="./data/heatmap.RData")
