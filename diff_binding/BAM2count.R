library(IRanges)
library(ShortRead)
library(rtracklayer)
library(GenomicAlignments)

srbam2count <- function(bam.file, extend=200) {
  require(ShortRead)
  gr <- readGAlignments(bam.file)
  grs <- as(gr, "GRanges")
  grsr <- resize(grs, extend)
  return(grsr)
}

args <- commandArgs(trailingOnly = TRUE)
file <- args[[1]]
fragment.length <- as.integer(args[[2]])

counts <- srbam2count(file, extend=fragment.length)
save(counts, file=paste(file,".counts.rda",sep=""))