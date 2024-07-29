#metabolites-protein 
library(clusterProfiler)
mp<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/meta-network.xlsx",
                        sheet = "Sheet1",
                        range = "A1:F31207", 
                        col_names = T,
                        na="NA")
mptable<- as.data.frame(table(mp$a))
mp <- mp %>% dplyr::filter (!a %in% c("Water","Hydrogen Ion","Adenosine triphosphate",
                                      "ADP","Phosphate","Pyrophosphate","Adenosine monophosphate",
                                      "NAD","NADP","NADPH","Oxygen","NADH","Coenzyme A","ubiquitin",
                                      "Protein lysine","O-phosphoprotamine","protamine"))

mprot<- as.data.frame(table(mp$b))
mptable2<- as.data.frame(table(mp$aids))

BiocManager::install("KEGGREST")
BiocManager::install("EnrichmentBrowser")
library("KEGGREST")
library("EnrichmentBrowser")

kegghsa<- keggList("hsa")
#step 3: download the pathways of that organism:
hsapathway <- downloadPathways("hsa")
#step 4: retrieve gene sets for an organism from databases such as GO and KEGG:
hsa <- getGenesets(org = "hsa", db = "kegg", cache = TRUE, return.type="list")

i=1
gs <- data.frame(geneset=names(hsapathway)[1],ENTREZID=hsa[[1]])
for(i in 2:247){
  b <- data.frame(geneset=names(hsapathway)[i],ENTREZID=hsa[[i]])
  gs <- rbind(gs,b)}

gene <- gs$ENTREZID
gene <- bitr(gene,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
gs <- merge(gs,gene,by="ENTREZID")

hsa2 <- read.table("/HDD8T/eunji/proj/lg/luma/hsa2.txt")
hsa3 <- read.table("/HDD8T/eunji/proj/lg/luma/hsa3.txt")
hsa <- cbind(hsa3,hsa2)
names(hsa) <- c('Hsa',"Pathway")
hsa4<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/hsa-met.xlsx",
                        sheet = "Sheet1",
                        range = "A2:C96", 
                        col_names = T,
                        na="NA")
hsa <- merge(hsa,hsa4,by="Hsa")
names(gs)[2] <- 'Hsa'
gs2 <- merge(gs,hsa,by="Hsa")

i=1
for(i in 2:1853){
  h=as.character(mptable2$Var1[i])
  a <- mp %>% dplyr::filter(aids == h)
  b <- a[grep(paste0("^",h,collapse =""),a$direction),]
  b <- b[-c(grep(';',b$direction)),]
  b$PS <- "DownPs"
  c <- a[c(grep(paste0(h,"$",collapse =""),a$direction)),]
  c <- c[-c(grep(';',c$direction)),]
  c$PS <- "UpPs"
  d <- rbind(b,c)
  e <- rbind(e,d)}
etable <- as.data.frame(table(e$a))
mpgene <- as.data.frame(table(e$b))
save(e,file = "/HDD8T/eunji/proj/lg/luma/e.RData")
load("/HDD8T/eunji/proj/lg/luma/tT.RData")
names(e)[2] <- "Gene.symbol"
geo <- merge(e,tT,by="Gene.symbol")
names(e)[2] <- "b"
geo2 <- geo[,c(1,2,5,7,25,26,29)]
geo2 <- geo2 %>% dplyr::filter(P.Value < 0.05)

i=1
mpi <- data.frame()
for(i in 1:710){
  m <- etable[i,1]
  mpi[i,'Metabolite'] <- etable[i,1]
  up <- geo2 %>% dplyr::filter(a == m & PS == "UpPs")
  up <- sum(up$logFC)
  mpi[i,'delta(UpPs)'] <- up
  down <- geo2 %>% dplyr::filter(a == m & PS == "DownPs")
  down <- sum(down$logFC)
  mpi[i,'delta(DownPs)'] <- down
  mpi[i,'deltaM'] <- mpi[i,'delta(UpPs)'] - mpi[i,'delta(DownPs)']
}
mpi <- mpi[order(mpi$deltaM,decreasing = T),]
mpi2 <- mpi %>% dplyr::filter(deltaM > 0.5 | deltaM < -0.5)
write.csv(mpi,file = "/HDD8T/eunji/proj/lg/luma/mpi.csv")
write.csv(mpi2,file = "/HDD8T/eunji/proj/lg/luma/mpi2.csv")


library(igraph)
glut <- geo2 %>% dplyr::filter(a == "Glutathione")
names(glut)[1] <- 'from'
names(glut)[2] <- 'to'

edges <- glut
nrow(unique(edges[,c("from", "to")]))
edges1 <- edges %>% dplyr::filter(PS == "DownPs")
names(edges1)[1] <- 'to'
names(edges1)[2] <- 'from'
edges2 <- edges %>% dplyr::filter(PS=="UpPs")
edges <- rbind(edges2,edges1)
edges$logp <- -log10(edges$adj.P.Val)
names(edges)[7] <- 'weight'
node <- glut[,c(1,4)]
names(node) <- c('id','type.label')
node2 <- data.frame(id=unique(glut$to),type.label="metabolite")
node <- rbind(node,node2)
node$type <- c(rep("protein",nrow(glut)),"metabolite")


nodes <- read.csv("/HDD8T/eunji/proj/lg/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("/HDD8T/eunji/proj/lg/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL

net <- graph_from_data_frame(d=edges, vertices=node, directed=T) 
net
E(net)$width <- E(net)$weight
E(net)$arrow.size <-0.7
plot(net, edge.arrow.size=.5, edge.curved=0,
     vertex.color=c( "#BC3C29FF","#3C5488FF")[1+(V(net)$type=="protein")],
     vertex.frame.color="#555555",
     vertex.label=V(net)$id, vertex.label.color="black",vertex.label.dist=2.5,
     vertex.label.cex=.8) 

mpifa<- readxl::read_excel("/HDD8T/eunji/proj/lg/luma/MPI.xlsx",
                        sheet = "Sheet1",
                        range = "B1:F385", 
                        col_names = T,
                        na="NA")
mpifa <- mpifa[!is.na(mpifa$FA),] 
mpiunfa <- mpifa %>% dplyr::filter(FA=="Unsaturated fatty acid")
mpiunfa <- mpiunfa[order(mpiunfa$deltaM,decreasing = T),]
mpiunfa$Metabolite <- factor(mpiunfa$Metabolite,levels = mpiunfa$Metabolite)

mpiunfa <- within(mpiunfa, 
                   Metabolite <- factor(Metabolite, 
                                      levels=sort(deltaM),decreasing=TRUE))

ggplot(mpiunfa, aes(x = reorder(Metabolite,deltaM), y = deltaM)) +
  geom_bar(stat = "identity",fill="steelblue")+ggtitle("Unsaturated Fatty Acids")+
  theme_bw()+ylab("deltaM=delta(UpPs)-delta(DownPs)")+
  theme(axis.text.x = element_text(size = 12,color = "black",angle = 90,vjust = 0.5, hjust=1),                      # adjusting the position
        axis.title.x = element_text(size = 12,color = "black"),                   # face the x axit title/label
        axis.text.y = element_text(size = 10,color = "black"),                   # face the y axis title/label
        axis.title.y = element_text(hjust = 0.1,face = "bold",size = 13))

#TCGA MPI 
load("/HDD8T/eunji/proj/lg/luma/lumA.RData")
luma <- as.data.frame(t(data.frame(row.names = lumA$sampleID,lumA[,12:49314])))
table(lumA$epas1)
luma$c1 <- apply(luma[,1:216],1,mean)
luma$c2 <- apply(luma[,217:431],1,mean)

luma$fc <- luma$c1/luma$c2
luma$logfc <- log2(luma$c1/luma$c2)

for(h in 1:nrow(luma)){
  x=as.numeric(luma[h, 1:216])
  y=as.numeric(luma[h, 217:431])
  luma$p[h]=t.test(x, y,alternative = "two.sided")$p.value}
resluma <- data.frame(Gene.symbol=rownames(luma),luma[,432:436])


names(e)[2] <- "Gene.symbol"
lumampi <- merge(e,resluma,by="Gene.symbol")
names(e)[2] <- "b"
geo2 <- geo[,c(1,2,5,7,25,29)]
lumampi2 <- lumampi %>% dplyr::filter(p < 0.05)

i=1
lumampi3 <- data.frame()
for(i in 1:710){
  m <- etable[i,1]
  lumampi3[i,'Metabolite'] <- etable[i,1]
  up <- lumampi2 %>% dplyr::filter(a == m & PS == "UpPs")
  up <- sum(up$logfc)
  lumampi3[i,'delta(UpPs)'] <- up
  down <- lumampi2 %>% dplyr::filter(a == m & PS == "DownPs")
  down <- sum(down$logfc)
  lumampi3[i,'delta(DownPs)'] <- down
  lumampi3[i,'deltaM'] <- lumampi3[i,'delta(UpPs)'] - lumampi3[i,'delta(DownPs)']
}
lumampi3 <- lumampi3[order(lumampi3$deltaM,decreasing = T),]
lumampi4 <- lumampi3 %>% dplyr::filter(deltaM > 0.5 | deltaM < -0.5)

library(igraph)
BiocManager::install("networkD3")
library(networkD3)
glut <- lumampi2 %>% dplyr::filter(a == "Hydrogen peroxide")
names(glut)[1] <- 'from'
names(glut)[2] <- 'to'

simpleNetwork(glut, height="100px", width="100px",        
              Source = 1,                 # column number of source
              Target = 2,                 # column number of target
              linkDistance = 10,          # distance between node. Increase this value to have more space between nodes
              charge = -900,                # numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)
              fontSize = 14,               # size of the node names
              fontFamily = "serif",       # font og node names
              linkColour = "#666",        # colour of edges, MUST be a common colour for the whole graph
              nodeColour = "#69b3a2",     # colour of nodes, MUST be a common colour for the whole graph
              opacity = 0.9,              # opacity of nodes. 0=transparent. 1=no transparency
              zoom = T                    # Can you zoom on the figure?
)
simpleNetwork(glut, height="100px", width="100px")

edges <- glut
nrow(unique(edges[,c("from", "to")]))
edges1 <- edges %>% dplyr::filter(PS == "DownPs")
names(edges1)[1] <- 'to'
names(edges1)[2] <- 'from'
edges2 <- edges %>% dplyr::filter(PS=="UpPs")
edges <- rbind(edges2,edges1)
edges$logp <- -log10(edges$p)
names(edges)[13] <- 'weight'
node <- glut[,c(1,7)]
names(node) <- c('id','type.label')
node2 <- data.frame(id=unique(glut$to),type.label="metabolite")
node <- rbind(node,node2)
node$type <- c(rep("protein",nrow(glut)),"metabolite")

link2 <- edges[,c(1,2,13)]


nodes <- read.csv("/HDD8T/eunji/proj/lg/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("/HDD8T/eunji/proj/lg/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
links <- aggregate(links[,3], links[,-3], sum)
links <- links[order(links$from, links$to),]
colnames(links)[4] <- "weight"
rownames(links) <- NULL

net <- graph_from_data_frame(d=edges, vertices=node, directed=T)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net

E(net)$width <- E(net)$weight
E(net)$arrow.size <-1
V(net)$size <- 30
plot(net, edge.arrow.size=.5, edge.curved=0,
     vertex.color=c( "firebrick3","darkgray")[1+(V(net)$type=="protein")],
     vertex.frame.color="#555555",
     vertex.label=V(net)$id, vertex.label.color="black",vertex.label.dist=2.5,
     vertex.label.cex=0.8) 

mpimerge <- merge(mpi,lumampi3,by="Metabolite")
mpimerge2 <- mpimerge %>% dplyr::filter(deltaM.x > 0 & deltaM.y > 0)
mpimerge3 <- mpimerge %>% dplyr::filter(deltaM.x < 0 & deltaM.y < 0)
mpimerge4 <- rbind(mpimerge2,mpimerge3)
mpimerge4 <- mpimerge4[,c(1,4,7)]
mpimerge4 <- data.frame(row.names = mpimerge4$Metabolite,mpimerge4[,c(2,3)])
mpimerge4[mpimerge4$deltaM.x > 2.5] =2.5
mpimerge4[mpimerge4 < -2.5]=-2.5
names(mpimerge4) <- c('MCF7','LumA')
mpimerge4 <- mpimerge4[order(mpimerge4$MCF7,decreasing=T),]
pheatmap(mpimerge4)
pheatmap(p, cluster_cols =T ,cluster_rows = T,
         clustering_distance_row = "correlation",
         fontsize_col=8,fontsize_row=8, cellwidth = 8,cellheight = 11,
         color=myColor,legend = T,breaks = myBreaks,top_annotation=ha,
         show_colnames = T,scale = 'none')

paletteLength <- 50
myColor <- colorRampPalette(c('blue','white','firebrick'))(paletteLength)
myBreaks <- c(seq(min(h),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(h)/paletteLength, max(h), length.out=floor(paletteLength/2)))

library(pheatmap)
pheatmap(mpimerge4, cluster_cols =F ,cluster_rows = F,scale = "none",
         clustering_distance_row = "correlation",
         fontsize_col=10,#fontsize_row=8, cellwidth = 8,cellheight = 11,
         color=myColor,legend = T,#breaks = myBreaks,
         #annotation_row = hypoglu,annotation_col = hgcol,
         show_colnames = T,show_rownames = F)
write.csv(mpimerge2,file = "/HDD8T/eunji/proj/lg/mpimerge2.csv")



