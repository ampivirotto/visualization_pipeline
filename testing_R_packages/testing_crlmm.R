library(crlmm)
library(ff)
library(GEOquery)

downloadDir <- tempdir()
getGEOSuppFiles("GSE27073", makeDirectory = FALSE, baseDir = downloadDir)
idatFiles <- list.files(path = downloadDir, pattern = ".idat.gz", full.names = TRUE)
sapply(idatFiles, gunzip)
idatFiles <- list.files(path = downloadDir, pattern = ".idat", full.names = TRUE)

samplesheet <- read.csv()
arrayNames <- file.path(idatFiles, unique(samplesheet[, "SentrixPosition"]))
arrayInfo <- list(barcode=NULL, position="SentrixPosition")
dataset <- genotype.Illumina(sampleSheet=samplesheet,
                             arrayNames=arrayNames,
                             arrayInfoColNames=arrayInfo,
                             cdfName="human370v1c",
                             batch=rep("1", nrow(samplesheet)))
