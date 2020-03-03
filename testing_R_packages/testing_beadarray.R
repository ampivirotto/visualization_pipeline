library(beadarray)
library(GEOquery)
library(crlmm)

downloadDir <- tempdir()
getGEOSuppFiles("GSE27073", makeDirectory = FALSE, baseDir = downloadDir)
idatFiles <- list.files(path = downloadDir, pattern = ".idat.gz", full.names = TRUE)
sapply(idatFiles, gunzip)
idatFiles <- list.files(path = downloadDir, pattern = ".idat", full.names = TRUE)

#genGSE <- genotype.Illumina(idatFiles, cdfName = "nopackage", call.method = "krlmm")

data <- readIdatFiles(idatFiles)
