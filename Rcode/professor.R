#20220511 -해부학중간고사채점,histogram,SD,mean 그림 

h <- hist(his$중간고사,
     breaks = 7,probability = T,
     col = "steelblue",
     main = "중간고사점수분포",
     ylab = "Prob",
     xlab = "Score")
lines(density(his$중간고사, bw = 1.199),lwd = 2)
abline(v = mean(his$중간고사),                       # Add line for mean
       col = "red",
       lwd = 1)
text(12.5, 0.13, paste("Mean =", round(mean(his$중간고사), 2), "\n Median =", 
                     round(median(his$중간고사), 2), "\n Std.Dev =", round(sd(his$중간고사), 2)))
text(h$mids,h$density-0.004,labels=h$counts, adj=c(0.5, -0.5))
