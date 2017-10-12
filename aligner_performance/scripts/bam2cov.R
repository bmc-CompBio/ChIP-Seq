rm(list=ls())
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
options(scipen=100)



for (f in list.files(".","*.bam$")) {
	gc()
 	bam.file <- f
	cat("reading", f, "\n")
	gr <- readGAlignments(bam.file)
	grs <- as(gr, "GRanges")
	grsr <- resize(grs, 200)
	covs <- coverage(grsr)
	saveRDS(covs, file=sub(".bam", "_covs.rds",f))
	gd <- as(covs, "RangedData")
	export.bedGraph(gd, sub(".bam",".wig",f))
}
