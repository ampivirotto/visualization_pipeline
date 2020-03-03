################  IMPORT STATEMENTS ####################
library(GEOquery)
library(tools)
library(crlmm)
library(ff)
#library(beadarray)

################ FUNCTIONS #########################

## implementing these things in python for now 

# getColumns <- function(filename){
#   filename <- file_path_sans_ext(basename(filename))
#   x <- unlist(strsplit(filename, "_"))
#   return (x)
# }
# 
# makeSampleSheet <- function(listofFiles, wd){
#   wddir <- paste(wd, "\\samplesheet.csv", sep = "")
#   fileConn <- file(wddir)
#   write("Sample_id, Sentrix_barcode, sentrix_position", fileConn)
#   for(i in 1:length(listofFiles)){
#     
#     newlist <- getColumns(listofFiles[i])
#     if(identical(newlist[4], "Red")){
#        line <- c(newlist[1], newlist[2], newlist[3])
#        line <- paste(line, collapse = ", ")
#        write(line, fileConn)
#        #samplesheet <- rbind(samplesheet, newlist)
#     }
#   }
#   close(fileConn)
#   samplesheet = read.csv(wddir)
#   return(samplesheet)
# }


pullGEOdata <- function(geonum){
  # may want to allow user to specifiy directory 
  downloadDir <- "D:\\visualization_pipeline\\data\\beadarray_data"
  getGEOSuppFiles(geonum, makeDirectory = FALSE, baseDir = downloadDir)
  idatFiles <- list.files(path = downloadDir, pattern = ".idat.gz", full.names = TRUE)
  sapply(idatFiles, gunzip)
  idatFiles <- list.files(path = downloadDir, pattern = ".idat", full.names = TRUE)
  #print(idatFiles)
  return(idatFiles)
}

retrieveDataGZ <- function(dirpath){
  idatFiles <- list.files(path = dirpath, pattern = ".idat.gz", full.names = TRUE)
  sapply(idatFiles, gunzip)
  idatFiles <- list.files(path = dirpath, pattern = ".idat", full.names = TRUE)
  return(idatFiles)
}

retrieveData <- function(dirpath){
  downloadDir <- dirpath
  idatFiles <- list.files(path = downloadDir, pattern = ".idat", full.names = TRUE)
  return(idatFiles)
}

##################### MAIN  ############################

################## TESTING #################################
geonum <- "GSE27073"

list_of_idat_files <- pullGEOdata(geonum)

wddir <- "D:\\visualization_pipeline\\data\\canine_snps"
list_of_idat_files <- retrieveData(wddir)

#samplesheet <- makeSampleSheet(list_of_idat_files, wddir)
#############################################################

samplesheet <- read.csv("D:\\visualization_pipeline\\data\\canine_snps\\samplesheet.csv")
arrayNames <- file.path(list_of_idat_files, unique(samplesheet[, "sentrix_position"]))
arrayInfo <- list(barcode=NULL, position="SentrixPosition")
dataset <- genotype.Illumina(sampleSheet=samplesheet,
                             arrayNames= list_of_idat_files,
                             path = wddir, 
                             arrayInfoColNames= list(barcode="sentrix_barcode", position="sentrix_position"), 
                             batch=rep("1", nrow(samplesheet)))
