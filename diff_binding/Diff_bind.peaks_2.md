# Differential Binding at Peaks


# Get sequencing file and align to dm6
Download sra-file and convert to fastq-file.

```r
module load ngs/sratoolkit/2.8.0

prefetch SRR520444 #Input_Kc_ctrl_1
prefetch SRR520445 #CLAMP_Kc_ctrl_1
prefetch SRR520446 #Input_Kc_ctrl_2
prefetch SRR520447 #CLAMP_Kc_ctrl_2
prefetch SRR520448 #Input_S2_msl2_1
prefetch SRR520449 #CLAMP_S2_msl2_1
prefetch SRR520450 #Input_S2_msl2_2
prefetch SRR520451 #CLAMP_S2_msl2_2
prefetch SRR520452 #Input_S2_ctrl_1
prefetch SRR520453 #CLAMP_S2_ctrl_1
prefetch SRR520454 #Input_S2_ctrl_2
prefetch SRR520455 #CLAMP_S2_ctrl_2

fastq-dump SRR520444 #Input_Kc_ctrl_1
fastq-dump SRR520445 #CLAMP_Kc_ctrl_1
fastq-dump SRR520446 #Input_Kc_ctrl_2
fastq-dump SRR520447 #CLAMP_Kc_ctrl_2
fastq-dump SRR520448 #Input_S2_msl2_1
fastq-dump SRR520449 #CLAMP_S2_msl2_1
fastq-dump SRR520450 #Input_S2_msl2_2
fastq-dump SRR520451 #CLAMP_S2_msl2_2
fastq-dump SRR520452 #Input_S2_ctrl_1
fastq-dump SRR520453 #CLAMP_S2_ctrl_1
fastq-dump SRR520454 #Input_S2_ctrl_2
fastq-dump SRR520455 #CLAMP_S2_ctrl_2
```

Align fastq-file with bowtie1 (with -m 1 option) to *Drosophila* genome dm6.
(here, reads mapped with quality score < 10 were removed, not necessary)

```r
module load ngs/bowtie1
module load ngs/samtools

bowtie -S -p 18 -m 1 /work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome ${FILENAME}.fastq | samtools view -bS -q 10 - | samtools sort - | tee ${FILENAME}.sorted.bam | samtools index - ${FILENAME}.sorted.bam.bai
```

# Create count-file from bam-file - BAM2count.R
Read in bam-files, convert to GRanged Object, extend read to 200 bp and save as .rda-file using the BAM2count.R skript.
(Usage: BAM2count.R <bam file> <fragment.length>,
Example: BAM2count.R input_gst_1 200)

```r
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
```

# Compute log2(IP/input) enrichment for peaks
Load required packages and load peak file (bed format), seqnames must be changed (e.g. '2R' > 'chr2R') to match count file seqnames; and convert to GRanged Object.

```r
options(scipen=100)
rm(list=ls())

suppressPackageStartupMessages({
library(rtracklayer)
library(limma)
library(RColorBrewer)
})

HAS <- read.delim("msl2_in_vivo_peaks.bed", header = F)
HAS[,1] <- paste0('chr',HAS[,1])
colnames(HAS) <- c('chr','start','end')
HAS <- as(HAS, 'GRanges')

kc_peaks <- import.bed('CLAMP_Kc_ctrl_robust_peaks.bed')
s2_peaks <- import.bed('CLAMP_S2_ctrl_robust_peaks.bed')

overlaps <- queryHits(findOverlaps(s2_peaks,kc_peaks))
s2_peaks$kc <- 'n'
s2_peaks$kc[overlaps] <- 'y'

overlaps <- queryHits(findOverlaps(s2_peaks,HAS))
s2_peaks$has <- 'n'
s2_peaks$has[overlaps] <- 'y'

clamp_peaks <- subset(s2_peaks, kc == 'y' | has == 'y')
```

Load count files in for loop and count overlapping reads at peaks with minimal overlap of 100 bp (half read length) and store overall number of aligned reads.

```r
files <- system("ls *counts.rda",intern=T)
names <- unlist(strsplit(files,split = "_counts.rda"))[seq(1,length(files)*2,2)]

for(i in seq_along(files)){
  load(files[i])
  if(i==1){
    counts <- countOverlaps(clamp_peaks, counts, minoverlap = 100, ignore.strand=T)
    reads <- length(counts)
  }else{
    counts <- cbind(counts,countOverlaps(clamp_peaks, counts, minoverlap = 100, ignore.strand=T))
    reads <- c(reads,length(counts))
  }
}

colnames(counts) <- names
saveRDS(counts, "clamp_counts.rds")
names(reads) <- names
saveRDS(reads, "clamp_reads.rds")
```

Calculate size factor to acount for different sequencing depth. Size Factor is calcuted as propoartin of number of reads for each sample to the mean read number.

```r
counts <- readRDS("clamp_counts.rds")
head(counts)
```

```
##      CLAMP_Kc_ctrl_1 CLAMP_Kc_ctrl_2 CLAMP_S2_ctrl_1 CLAMP_S2_ctrl_2
## [1,]            4833            3695            2948            1133
## [2,]            4354            3114            2738             868
## [3,]            1981             725             548             214
## [4,]            4834            2208            1578             603
## [5,]            2462            1807             892             412
## [6,]             602             278             608             204
##      CLAMP_S2_msl2_1 CLAMP_S2_msl2_2 Input_Kc_ctrl_1 Input_Kc_ctrl_2
## [1,]            4293            2121              68              46
## [2,]            3369            1613             112              54
## [3,]             757             253             140              43
## [4,]            2413            1134             123              65
## [5,]            1847             907              74              34
## [6,]             643             235             107              49
##      Input_S2_ctrl_1 Input_S2_ctrl_2 Input_S2_msl2_1 Input_S2_msl2_2
## [1,]              31              24             157              65
## [2,]              59              41             123              99
## [3,]              46              33             126              58
## [4,]              60              39             225              96
## [5,]              39              33             140              62
## [6,]              83              42             347              83
```

```r
reads <- readRDS("clamp_reads.rds")
reads
```

```
## CLAMP_Kc_ctrl_1 CLAMP_Kc_ctrl_2 CLAMP_S2_ctrl_1 CLAMP_S2_ctrl_2 
##        53212323        18194728        19577866        15774853 
## CLAMP_S2_msl2_1 CLAMP_S2_msl2_2 Input_Kc_ctrl_1 Input_Kc_ctrl_2 
##        40343206        18738524        47431197        16107592 
## Input_S2_ctrl_1 Input_S2_ctrl_2 Input_S2_msl2_1 Input_S2_msl2_2 
##        18031415        13960771        42648596        23732843
```

```r
sizeFactor <- reads/mean(reads)
sizeFactor
```

```
## CLAMP_Kc_ctrl_1 CLAMP_Kc_ctrl_2 CLAMP_S2_ctrl_1 CLAMP_S2_ctrl_2 
##       1.9482540       0.6661606       0.7168012       0.5775621 
## CLAMP_S2_msl2_1 CLAMP_S2_msl2_2 Input_Kc_ctrl_1 Input_Kc_ctrl_2 
##       1.4770791       0.6860705       1.7365906       0.5897446 
## Input_S2_ctrl_1 Input_S2_ctrl_2 Input_S2_msl2_1 Input_S2_msl2_2 
##       0.6601812       0.5111434       1.5614860       0.8689267
```

```r
counts.norm <- t(t(counts) / sizeFactor)
head(counts.norm)
```

```
##      CLAMP_Kc_ctrl_1 CLAMP_Kc_ctrl_2 CLAMP_S2_ctrl_1 CLAMP_S2_ctrl_2
## [1,]       2480.6827       5546.7107       4112.7164       1961.6938
## [2,]       2234.8215       4674.5486       3819.7482       1502.8687
## [3,]       1016.8079       1088.3262        764.5077        370.5229
## [4,]       2481.1960       3314.5162       2201.4473       1044.0436
## [5,]       1263.6956       2712.5592       1244.4176        713.3432
## [6,]        308.9946        417.3168        848.2129        353.2088
##      CLAMP_S2_msl2_1 CLAMP_S2_msl2_2 Input_Kc_ctrl_1 Input_Kc_ctrl_2
## [1,]       2906.4116       3091.5191        39.15719        77.99987
## [2,]       2280.8527       2351.0704        64.49419        91.56506
## [3,]        512.4979        368.7668        80.61773        72.91292
## [4,]       1633.6295       1652.8914        70.82844       110.21720
## [5,]       1250.4408       1322.0216        42.61223        57.65207
## [6,]        435.3186        342.5304        61.61498        83.08681
##      Input_S2_ctrl_1 Input_S2_ctrl_2 Input_S2_msl2_1 Input_S2_msl2_2
## [1,]        46.95680        46.95355       100.54525        74.80493
## [2,]        89.36940        80.21232        78.77112       113.93367
## [3,]        69.67784        64.56114        80.69237        66.74902
## [4,]        90.88414        76.29953       144.09351       110.48113
## [5,]        59.07469        64.56114        89.65818        71.35240
## [6,]       125.72305        82.16872       222.22421        95.52014
```

```r
zeros <- which(apply(counts.norm,1,function(x){min(x) == 0}))
counts.norm <- counts.norm[-zeros,]
clamp_peaks <- clamp_peaks[-zeros,]
```

Calculate log2(IP/input) and plot log2 fold change *mle* RNAi - *gst* RNAi.

```r
enrichment <- log2(counts.norm[,1:6]/counts.norm[,7:12])
head(enrichment)
```

```
##      CLAMP_Kc_ctrl_1 CLAMP_Kc_ctrl_2 CLAMP_S2_ctrl_1 CLAMP_S2_ctrl_2
## [1,]        5.985316        6.152017        6.452614        5.384722
## [2,]        5.114847        5.673886        5.417553        4.227751
## [3,]        3.656806        3.899793        3.455759        2.520825
## [4,]        5.130563        4.910377        4.598280        3.774364
## [5,]        4.890238        5.556138        4.396787        3.465858
## [6,]        2.326229        2.328452        2.754177        2.103860
##      CLAMP_S2_msl2_1 CLAMP_S2_msl2_2
## [1,]       4.8533222        5.369039
## [2,]       4.8557627        4.367052
## [3,]       2.6670421        2.465890
## [4,]       3.5030035        3.903120
## [5,]       3.8018576        4.211640
## [6,]       0.9700556        1.842355
```

```r
S2.ctrl <- enrichment[,3:4]
S2.msl2 <- enrichment[,5:6]
Kc_ctrl <- enrichment[,1:2]

my_color <- brewer.pal(9, 'Set1')

par(mfrow = c(1,2))
plot(rowMeans(S2.ctrl[clamp_peaks$has == 'n',]),
     rowMeans(Kc_ctrl[clamp_peaks$has == 'n',]),
     main = "CLAMP enrichment at peaks",
     ylab = "Kc", xlab = "S2",
     pch=19, col=paste0(my_color[9],'5F'),
     ylim = c(-2, 8), xlim = c(-2,8))
points(rowMeans(S2.ctrl[clamp_peaks$has == 'y',]),
     rowMeans(Kc_ctrl[clamp_peaks$has == 'y',]),
     pch=19, col=paste0(my_color[1],'5F'))

plot(rowMeans(S2.ctrl[clamp_peaks$has == 'n',]),
     rowMeans(S2.msl2[clamp_peaks$has == 'n',]),
     main = "CLAMP enrichment at peaks",
     ylab = "msl2 RNAi S2", xlab = "ctrl RNAi S2",
     pch=19, col=paste0(my_color[9],'5F'),
     ylim = c(-2, 8), xlim = c(-2,8))
points(rowMeans(S2.ctrl[clamp_peaks$has == 'y',]),
     rowMeans(S2.msl2[clamp_peaks$has == 'y',]),
     pch=19, col=paste0(my_color[1],'5F'))
```

![](Diff_bind.peaks_2_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

# Test for significantly differnt log2(IP/input) enrichment values at peaks
To test for significant differences use the limma package.
Create 'Design Matrix' for fitting a linear model with lmFit() and calculate statistices with eBayes() allow to use intensity-trend for prior variance and allow the prior estimation to be robustified against outlier sample variance. Get ranked table with topTable() using 'fdr method' to adjust for multiple testing.

```r
my_RNAi <- c(rep("gst",2),rep("msl2",2))
my_RNAi <- factor(my_RNAi, levels = c("gst", "msl2"))
design <- model.matrix(~ my_RNAi)
design
```

```
##   (Intercept) my_RNAimsl2
## 1           1           0
## 2           1           0
## 3           1           1
## 4           1           1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$my_RNAi
## [1] "contr.treatment"
```

```r
fit <- lmFit(enrichment[,3:6], design)
fit <- eBayes(fit, trend=T, robust=T)
tt <- topTable(fit, coef=ncol(design), adjust.method = 'fdr', number = length(clamp_peaks))
head(tt)
```

```
##          logFC  AveExpr         t                     P.Value
## 1237 -4.889104 1.747855 -9.069834 0.0000000000000000002098279
## 1771 -4.385270 2.018660 -7.869340 0.0000000000000078035869444
## 1242 -3.811581 1.232733 -7.457879 0.0000000000001036146783542
## 1246 -3.965406 1.749444 -7.359410 0.0000000000002158801767901
## 1343 -3.638236 1.085316 -7.223390 0.0000000000005861248211237
## 1792 -3.669748 1.215142 -7.193055 0.0000000000007306182166773
##                     adj.P.Val        B
## 1237 0.0000000000000005079934 32.93456
## 1771 0.0000000000094462419962 22.90669
## 1242 0.0000000000836170454319 20.59157
## 1246 0.0000000001306614770022 19.88611
## 1343 0.0000000002838016383881 18.92706
## 1792 0.0000000002948044504293 18.71561
```



```r
fdr <- seq(0,1,0.001)
roc_mat <- matrix(0,nrow = length(fdr), ncol = 3)
colnames(roc_mat) <- c("tpr","fpr","fdr")
for(k in seq_along(fdr)){
  clamp_peaks$test <- "y"
  clamp_peaks$test[as.numeric(rownames(subset(tt, adj.P.Val < fdr[k])))] <- 'n'
  roc_mat[k,"tpr"] <- sum(clamp_peaks$kc == 'n' & clamp_peaks$test == 'n')/sum(clamp_peaks$kc == 'n')
    roc_mat[k,"fpr"] <- sum(clamp_peaks$kc == 'y' & clamp_peaks$test == 'n')/sum(clamp_peaks$kc == 'y')
    roc_mat[k,"fdr"] <- fdr[k]
}
roc_mat <- as.data.frame(roc_mat)
plot(tpr ~ fpr, roc_mat, type="l",lwd=3,col=my_color[1],
     ylab = 'true positive rate', xlab = 'false positive rate', main = 'ROC')
abline(0,1, lwd=2)
points(tpr ~ fpr, subset(roc_mat, fdr %in% c(0.01, 0.05, 0.1, 0.5)), pch = 16)
```

![](Diff_bind.peaks_2_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

Plot log2 fold change *msl2* RNAi - *ctrl* RNAi and highlight significantly different peaks (FDR < 0.1 and < 0.5).

```r
fdr0.5 <- clamp_peaks[as.numeric(rownames(subset(tt, adj.P.Val < 0.5)))]
fdr0.1 <- clamp_peaks[as.numeric(rownames(subset(tt, adj.P.Val < 0.1)))]

par(mfrow = c(1,2))
plot(rowMeans(S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val >= 0.5)))],
     rowMeans(S2.msl2-S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val >= 0.5)))],
     main = "CLAMP enrichment at peaks",
     ylab = "lfc (msl2 RNAi - ctrl RNAi)", xlab = "CLAMP enrichment",
     pch=19, col=paste0(my_color[9],'8F'), 
     ylim = c(-6, 6), xlim = c(0,8))
points(rowMeans(S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val < 0.5)))],
       rowMeans(S2.msl2-S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val < 0.5)))],
       pch=19, col = paste0(my_color[1],'8F'))
legend('bottomleft', bty = 'n', legend = paste0('fdr < 0.5, n=',nrow(subset(tt, adj.P.Val < 0.5))),
       pch = 19, col = my_color[1])

plot(rowMeans(S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val >= 0.1)))],
     rowMeans(S2.msl2-S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val >= 0.1)))],
     main = "CLAMP enrichment at peaks",
     ylab = "lfc (msl2 RNAi - ctrl RNAi)", xlab = "CLAMP enrichment",
     pch=19, col=paste0(my_color[9],'8F'), 
     ylim = c(-6, 6), xlim = c(0,8))
points(rowMeans(S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val < 0.1)))],
       rowMeans(S2.msl2-S2.ctrl)[as.numeric(rownames(subset(tt, adj.P.Val < 0.1)))],
       pch=19, col = paste0(my_color[1],'8F'))
legend('bottomleft', bty = 'n', legend = paste0('fdr < 0.1, n=',nrow(subset(tt, adj.P.Val < 0.1))), pch = 19, col = my_color[1])
```

![](Diff_bind.peaks_2_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


Convert count files to coverage vector and store in List.

```r
suppressPackageStartupMessages({
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
})

files <- system("ls CLAMP*counts.rda",intern=T)
names <- unlist(strsplit(files,split = ".fastq"))[seq(1,length(files)*2,2)]

covList <- list()
for(i in seq_along(files)){
  load(files[i])
  covList[[i]] <- coverage(counts.ip)/length(counts.ip)*10^6}
names(covList) <- names
saveRDS(covList, "clamp_covList.rds")
```

Plot Browser profile at region of interest using "tsTools"

```r
suppressPackageStartupMessages({
  library(ShortRead)
  library(rtracklayer)
  library(GenomicRanges)
  library(gdata)
  library(GenomicFeatures)
  library(tsTools)
  library(RColorBrewer)
})

covList <- readRDS("clamp_covList.rds")

txdb <- makeTxDbFromGFF("~/Desktop/mount/work/data/genomes/fly/Drosophila_melanogaster/UCSC/dm6/Annotation/Archives/archive-2015-07-24-09-25-49/Genes/genes.gtf")
```

```
## Import genomic features from the file as a GRanges object ... OK
## Prepare the 'metadata' data frame ... OK
## Make the TxDb object ... OK
```

```r
cols <- c(brewer.pal(4, "Greens")[3:4],brewer.pal(4, "Blues")[3:4],brewer.pal(4, "Reds")[3:4])
lim <- c(0,100)
lims <- rep(list(lim),6)

chr <- "chrX"
start <- 10*10^6
end <- 11*10^6

plotProfiles(start,end,chr,profs=covList,txdb=txdb,collapse=T, plot.labels=F, cols=cols, grid=T, ylims=lims,ylab="coverage (rpm)")
```

```
## Loading required package: grid
## Loading required package: HilbertVis
## Loading required package: lattice
```

![](Diff_bind.peaks_2_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


```r
peakList <- list()
peakList[[1]] <- coverage(subset(clamp_peaks, has == "y"))
peakList[[2]] <- coverage(subset(clamp_peaks, has == "y" & kc == "n"))
peakList[[3]] <- coverage(subset(fdr0.5, has == "y"))
peakList[[4]] <- coverage(subset(fdr0.1, has == "y"))
names(peakList) <- c('HAS','non-Kc','FDR < 0.5','FDR < 0.1')

cols <- brewer.pal(4, "Set1")
lim <- c(0,1)
lims <- rep(list(lim),5)

chr <- "chrX"
start <- 10*10^6
end <- 11*10^6

plotProfiles(start,end,chr,profs=peakList,txdb=txdb,collapse=T, plot.labels=F, cols=cols, grid=T, ylims=lims,ylab="peaks")
```

![](Diff_bind.peaks_2_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

