list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, ,repos = "http://cran.us.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.12")

is.installed <- function(mypkg){
    is.element(mypkg, installed.packages()[,1])
  } 

  # check if package "ggtree" is installed
  if (!is.installed("ggtree")){
    BiocManager::install("ggtree")
  }

library(tidyverse)
library(ggtree)

args = commandArgs(trailingOnly=TRUE)

setwd(args[2])

png(args[3])
tree <- read.tree(args[1])
labeled_tree <- ggtree(tree) + geom_tiplab()
labeled_tree
dev.off()
