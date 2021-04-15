
list.of.packages <- c("plotly")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, ,repos = "http://cran.us.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.12")

is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
  } 

  # check if package "SNPRelate" is installed
  if (!is.installed("SNPRelate")){
    BiocManager::install("SNPRelate")
  }
library("gdsfmt")
library("SNPRelate")
library("ggplot2")
library("plotly")

args = commandArgs(trailingOnly=TRUE)

#code from http://corearray.sourceforge.net/tutorials/SNPRelate/#principal-component-analysis-pca
setwd(args[1]) #set directory to where your vcf/gds
snpgdsVCF2GDS(args[2], args[4]) #converts vcf to gds given input.vcf, output.gds
gds <- snpgdsOpen(args[4]) #set the gds to a variable
pca <- snpgdsPCA(gds, num.thread=2,eigen.cnt = 0)
#?snpgdsPCA
pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))
#pc.percent
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    #the thrid eigenvector (if needed for 3D plot)
                  stringsAsFactors = FALSE)
#head(tab)
#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#scatter3D(tab$EV2,tab$EV1,tab$EV3) #3D scatterplot
fig <- plot_ly(tab, x=~EV2, y=~EV1, z=~EV3, name = ~sample.id) #create interactive 3D plot

fig <- fig %>% layout(
    title = "PCA",
    scene = list(
      xaxis = list(title = paste("EV1:", round(pc.percent[1], 2), "%")),
      yaxis = list(title = paste("EV2:", round(pc.percent[2], 2), "%")),
      zaxis = list(title = paste("EV3:", round(pc.percent[3], 2), "%"))))
    
fig #show plot
#htmlwidgets::saveWidget(as_widget(fig), args[3])
htmlwidgets::saveWidget(fig, args[3])
#ggplot(tab, aes(EV2,EV1)) + geom_point() +geom_text(aes(label=sample.id),hjust=0, vjust=0, size = 2) 
#?aes
