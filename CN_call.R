library(QDNAseq)
load("../data/cn/hg38.100kbp.SR50.rda")
hg38.100kbp.SR50@data[hg38.100kbp.SR50@data == "X"] = "23"
hg38.100kbp.SR50@data[hg38.100kbp.SR50@data == "Y"] = "24"

bins<- hg38.100kbp.SR50
readCounts <- binReadCounts(bins, path = "../data/cn/primary/")
readCounts <- applyFilters(readCounts, chromosomes = c("23", "24"))
readCounts <- estimateCorrection(readCounts)
readCounts <- applyFilters(readCounts, chromosomes = NA)
copyNumbers <- correctBins(readCounts)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)#no transformation
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented)
exportBins(copyNumbersCalled, file = "merged_pa_segmented_100kbp.tsv", format = "tsv", type = "segments", logTransform = FALSE)
