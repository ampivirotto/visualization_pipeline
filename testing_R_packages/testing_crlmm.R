library(crlmm)
library(ff)
library(GEOquery)

downloadDir <- "D:/visualization_pipeline/data/beadarray_data/GSE27073"
getGEOSuppFiles("GSE27073", makeDirectory = FALSE, baseDir = downloadDir)
idatFiles <- list.files(path = downloadDir, pattern = ".idat.gz", full.names = TRUE)
sapply(idatFiles, gunzip)
idatFiles <- list.files(path = downloadDir, pattern = ".idat", full.names = TRUE)

samplesheet <- read.csv(file = 'D:/visualization_pipeline/data/beadarray_data/GSE27073/samplesheet.csv')
arrayNames <- file.path(idatFiles, unique(samplesheet[, "sentrix_position"]))
arrayInfo <- list(barcode=NULL, position="sentrix_position")
dataset <- genotype.Illumina(sampleSheet=samplesheet,
                             arrayNames=arrayNames,
                             cdfName="human370v1c",
                             path = downloadDir)
