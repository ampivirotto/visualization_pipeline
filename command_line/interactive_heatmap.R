#! /usr/bin/Rscript
list.of.packages <- c("gplots", "vioplot", "heatmaply")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, ,repos = "http://cran.us.r-project.org")
#install.packages("vioplot",repos = "http://cran.us.r-project.org")
#install.packages("gplots",repos = "http://cran.us.r-project.org")
library(vioplot)
library(gplots)
library(heatmaply)

args = commandArgs(trailingOnly=TRUE) # https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/

setwd(args[1])

relat <- as.matrix(read.csv(args[2], header=TRUE, row.names=1))

ind <- colnames(relat)


#heatmap.2(relat, trace="none", Rowv=TRUE, Colv=TRUE, cexRow=0.6,cexCol = 0.6, labRow = ind, labCol = ind, col= colorRampPalette(c("lemonchiffon", "lemonchiffon", "lemonchiffon", "yellow", "red", "magenta", "midnightblue", "black"))(70))
fig <- heatmaply(relat)
htmlwidgets::saveWidget(fig, "interactive_heatmap.html")

#https://www.datanovia.com/en/blog/how-to-create-a-beautiful-interactive-heatmap-in-r/
