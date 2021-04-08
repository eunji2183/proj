#https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html
#https://www.r-graph-gallery.com/dendrogram.html#:~:text=A%20dendrogram%20(or%20tree%20diagram,to%20build%20one%20with%20R.

#dendrogram 
data <- matrix( sample(seq(1,2000),200), ncol = 10 )
rownames(data) <- paste0("sample_" , seq(1,20))
colnames(data) <- paste0("variable",seq(1,10))

# Euclidean distance
dist <- dist(data[ , c(4:8)] , diag=TRUE)

# Hierarchical Clustering with hclust
hc <- hclust(dist)

# Plot the result
plot(hc)

library(tidyverse)

# Data
head(mtcars)

# Clusterisation using 3 variables
mtcars %>% 
  select(mpg, cyl, disp) %>% 
  dist() %>% 
  hclust() %>% 
  as.dendrogram() -> dend

# Plot
par(mar=c(7,3,1,1))  # Increase bottom margin to have the complete label
plot(dend)

library(dendextend)

# Chart (left)
dend %>% 
  # Custom branches
  set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
  # Custom labels
  set("labels_col", "orange") %>% set("labels_cex", 0.8) %>%
  plot()
par(mar=c(1,1,1,7))
dend %>%
  set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(horiz=TRUE, axes=FALSE)
abline(v = 350, lty = 2)

my_colors <- ifelse(mtcars$am==0, "forestgreen", "green")

# Make the dendrogram
par(mar=c(10,1,1,1))
dend %>%
  set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  set("leaves_pch", 19)  %>% 
  set("nodes_cex", 0.7) %>% 
  plot(axes=FALSE)

# Add the colored bar
colored_bars(colors = my_colors, dend = dend, rowLabels = "am")

# Make 2 dendrograms, using 2 different clustering methods
d1 <- USArrests %>% dist() %>% hclust( method="average" ) %>% as.dendrogram()
d2 <- USArrests %>% dist() %>% hclust( method="complete" ) %>% as.dendrogram()

# Custom these kendo, and place them in a list
dl <- dendlist(
  d1 %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
  d2 %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
)

# Plot them together
tanglegram(dl, 
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=7,
           lwd=2)
