# data from 
# https://www.encodeproject.org/experiments/ENCSR000EVI/
# https://www.encodeproject.org/experiments/ENCSR000EVB/
# Location hg19 chr1 43824529
# CDC20 gene location
# hg38 chr1 43358955 43363203
# hg19 chr1 43824626 43828874
# chr1:43824626-43828874

## loading packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

peak = readPeakFile("~/biodata/CDC20/ENCFF001VHL.bed.gz")
peak = readPeakFile("~/biodata/CDC20/ENCFF002CRM.bed.gz")
peak

covplot(peak, weightCol = "V5")
covplot(peak, weightCol="V5", chrs=c("chr1"), xlim=c(4.3e7, 4.4e7))
covplot(peak, weightCol="V10", chrs=c("chr1"), xlim=c(43824500, 43828874))

subset(peak, seqnames == "chr1", starts > 4.3e7)

peak[seqnames(peak) == "chr1" & start(peak) > 4.3e7 & end(peak) < 4.4e7]
