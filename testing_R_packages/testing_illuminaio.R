library(illuminaio)
library(argyle)

setwd("D:/visualization_pipeline/canine_snps")
idat <- readIDAT("GSM3502154_3998887025_R06C01_Grn.idat")

names(idat)

idatData <- idat$Quants

head(idatData)
length(idatData)

idatField <- idat$fields

head(idatField)

data(ex)
summary(ex)
print(ex)
head(ex)
markers(ex)
samples(ex)
ex <- run.sample.qc(ex, max.H = 533, max.N = 500)
summarize.filters(ex)

x <- read.table("D:\\visualization_pipeline\\data\\canine_snps\\GSE123368_series_matrix.txt")

genotypes("D:\\visualization_pipeline\\data\\canine_snps\\GSE123368_series_matrix.txt")
