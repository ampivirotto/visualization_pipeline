
list.of.packages <- c("vcfR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, ,repos = "http://cran.us.r-project.org")

library("vcfR")
setwd("/content/data/GSE83225")

vcf <- read.vcfR("test.vcf.gz", verbose = FALSE)
dna <- ape::read.dna("/content/reference/canFam2.fa", format = "fasta")
gff <- read.table("/content/reference/canFam3.1.103.gff3", sep="\t", quote="")

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=TRUE)

#chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
#chrom <- proc.chromR(chrom, verbose = TRUE)

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_DP = 300, max_DP = 700)
chrom <- proc.chromR(chrom, verbose = FALSE)
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)


png("chromo_plot.png")
#chromoqc(chrom, dp.alpha = 22)
#plot(chrom)
heatmap.bp(dp[1001:1500,])
dev.off()
