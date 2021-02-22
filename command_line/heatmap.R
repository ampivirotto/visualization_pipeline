#! /usr/bin/Rscript

install.packages("vioplot",repos = "http://cran.us.r-project.org")
install.packages("gplots",repos = "http://cran.us.r-project.org")
library(vioplot)
library(gplots)

args = commandArgs(trailingOnly=TRUE) # https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/

setwd(args[1])

relat <- as.matrix(read.csv(args[2], header=TRUE, row.names=1))

ind <- colnames(relat)

jpeg(args[3])
heatmap.2(relat, trace="none", Rowv=TRUE, Colv=TRUE, cexRow=0.6,cexCol = 0.6, labRow = ind, labCol = ind, col= colorRampPalette(c("lemonchiffon", "lemonchiffon", "lemonchiffon", "yellow", "red", "magenta", "midnightblue", "black"))(70))
dev.off()
