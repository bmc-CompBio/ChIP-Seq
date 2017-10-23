# Homer




## Peak Finding


```bash
#!/bin/sh
# homer.sh


# make pooled tag directory for bam files

makeTagDirectory E03pool_K36.dir ../BAM/E03A_K36_TTAGGCAT_1.bam ../BAM/E03B_K36_ACAGTGAT_1.bam
makeTagDirectory E11pool_K36.dir ../BAM/E11A_K36_GCCAATAT_1.bam ../BAM/E11B_K36_GATCAGAT_1.bam
makeTagDirectory E03pool_INP.dir ../BAM/E03A_INP_ATCACGAT_1.bam ../BAM/E03B_INP_CGATGTAT_1.bam
makeTagDirectory E11pool_INP.dir ../BAM/E11A_INP_TTAGGCAT_1.bam ../BAM/E11B_INP_ACAGTGAT_1.bam


# peak finding against input

findPeaks E03pool_K36.dir/  -i E03pool_INP.dir/ -style histone -C 0 -F 2 -o E03pool_K36.histone.F2.txt
findPeaks E11pool_K36.dir/  -i E11pool_INP.dir/ -style histone -C 0 -F 2 -o E11pool_K36.histone.F2.txt

findPeaks E03pool_K36.dir/  -i E03pool_INP.dir/ -style histone -C 0 -F 3 -o E03pool_K36.histone.F3.txt
findPeaks E11pool_K36.dir/  -i E11pool_INP.dir/ -style histone -C 0 -F 3 -o E11pool_K36.histone.F3.txt

findPeaks E03pool_K36.dir/  -i E03pool_INP.dir/ -region -size 500 -minDist 100 -C 0 -F 2 -o E03pool_K36.size500.minDist100.F2.txt
findPeaks E11pool_K36.dir/  -i E11pool_INP.dir/ -region -size 200 -minDist 100 -C 0 -F 2 -o E11pool_K36.size500.minDist100.F2.txt

findPeaks E03pool_K36.dir/  -i E03pool_INP.dir/ -region -size 500 -minDist 500 -C 0 -F 2 -o E03pool_K36.size500.minDist500.F2.txt
findPeaks E11pool_K36.dir/  -i E11pool_INP.dir/ -region -size 200 -minDist 500 -C 0 -F 2 -o E11pool_K36.size500.minDist500.F2.txt



# convert txt to BED

cat E03pool_K36.histone.F2.txt | grep -v "#" | cut -f 2,3,4 > E03pool_K36.histoneF2.bed
cat E11pool_K36.histone.F2.txt | grep -v "#" | cut -f 2,3,4 > E11pool_K36.histoneF2.bed

cat E03pool_K36.histone.F3.txt | grep -v "#" | cut -f 2,3,4 > E03pool_K36.histoneF3.bed
cat E11pool_K36.histone.F3.txt | grep -v "#" | cut -f 2,3,4 > E11pool_K36.histoneF3.bed

cat E03pool_K36.size500.minDist100.F2.txt | grep -v "#" | cut -f 2,3,4 > E03pool_K36.size500.minDist100.F2.bed
cat E11pool_K36.size500.minDist100.F2.txt | grep -v "#" | cut -f 2,3,4 > E11pool_K36.size500.minDist100.F2.bed

cat E03pool_K36.size500.minDist500.F2.txt | grep -v "#" | cut -f 2,3,4 > E03pool_K36.size500.minDist500.F2.bed
cat E11pool_K36.size500.minDist500.F2.txt | grep -v "#" | cut -f 2,3,4 > E11pool_K36.size500.minDist500.F2.bed



# all pooled for -F 2

cat E03pool_K36.histoneF2.bed E11pool_K36.histoneF2.bed >  Eallpool_K36.histoneF2.bed
sort -k1,1 -k2,2n Eallpool_K36.histoneF2.bed > Eallpool_K36.histoneF2.sorted.bed
mergeBed -i Eallpool_K36.histoneF2.sorted.bed > Eallpool_K36.histoneF2.merged.bed

cat E03pool_K36.size500.minDist100.F2.bed E11pool_K36.size500.minDist100.F2.bed >  Eallpool_K36.size500.minDist100.F2.bed
sort -k1,1 -k2,2n Eallpool_K36.size500.minDist100.F2.bed > Eallpool_K36.size500.minDist100.F2.sorted.bed
mergeBed -i Eallpool_K36.size500.minDist100.F2.sorted.bed > Eallpool_K36.size500.minDist100.F2.merged.bed

cat E03pool_K36.size500.minDist500.F2.bed E11pool_K36.size500.minDist500.F2.bed >  Eallpool_K36.size500.minDist500.F2.bed
sort -k1,1 -k2,2n Eallpool_K36.size500.minDist500.F2.bed > Eallpool_K36.size500.minDist500.F2.sorted.bed
mergeBed -i Eallpool_K36.size500.minDist500.F2.sorted.bed > Eallpool_K36.size500.minDist500.F2.merged.bed
```



## BAM Files


```r
library(GenomicAlignments)
library(GenomicRanges)
library(csaw)
library(edgeR)

my_ChIP <- "K36"
my_control <- "INP"

bam.files <- list.files("../BAM/", pattern = ".bam$")
bam.files <- bam.files[grep(paste(my_ChIP,my_control,sep="|"), bam.files)]
bam.files.path <- file.path("../BAM/", bam.files)

my_samples <- gsub("_[G,A,T,C].*","",bam.files)

print(bam.files)
```

```
## [1] "E03A_INP_ATCACGAT_1.bam" "E03A_K36_TTAGGCAT_1.bam"
## [3] "E03B_INP_CGATGTAT_1.bam" "E03B_K36_ACAGTGAT_1.bam"
## [5] "E11A_INP_TTAGGCAT_1.bam" "E11A_K36_GCCAATAT_1.bam"
## [7] "E11B_INP_ACAGTGAT_1.bam" "E11B_K36_GATCAGAT_1.bam"
```

```r
print(my_samples)
```

```
## [1] "E03A_INP" "E03A_K36" "E03B_INP" "E03B_K36" "E11A_INP" "E11A_K36"
## [7] "E11B_INP" "E11B_K36"
```



## Counting





```r
my_peaks <- read.delim("Eallpool_K36.size500.minDist500.F2.merged.bed", header = FALSE)
colnames(my_peaks) <- c("chr","start","end")

my_peak_ranges <- makeGRangesFromDataFrame(my_peaks)

my_peak_counts <- my_peak_ranges

for(i in seq_along( bam.files.path)){

    my_peak_counts <- bam2bins(bam.file =  bam.files.path[i],
                                 bins = my_peak_counts,
                                 my_sample = my_samples[i]
                                )
}


saveRDS(my_peak_counts, paste(my_ChIP,"counts.rds", sep="."))

my_peak_counts
```

```
## GRanges object with 6059 ranges and 8 metadata columns:
##          seqnames               ranges strand |  E03A_INP  E03A_K36
##             <Rle>            <IRanges>  <Rle> | <integer> <integer>
##      [1]    chr2L      [ 7606,   8106]      * |        92       197
##      [2]    chr2L      [ 8503,  16061]      * |      1401      6868
##      [3]    chr2L      [67842,  72637]      * |       897      3353
##      [4]    chr2L      [83112,  83612]      * |       115       215
##      [5]    chr2L      [94968, 102319]      * |      1450      5684
##      ...      ...                  ...    ... .       ...       ...
##   [6055]     chrX [23087925, 23090444]      * |       344       685
##   [6056]     chrX [23097487, 23102483]      * |       686      2921
##   [6057]     chrX [23105988, 23108518]      * |       372      1586
##   [6058]     chrY [ 3258929,  3259429]      * |       143       254
##   [6059]     chrY [ 3606553,  3606803]      * |        66       106
##           E03B_INP  E03B_K36  E11A_INP  E11A_K36  E11B_INP  E11B_K36
##          <integer> <integer> <integer> <integer> <integer> <integer>
##      [1]       106       353       102       193       166       164
##      [2]      1408     14053      1539      4059      1881      4153
##      [3]       900      6139       955      2315      1169      2395
##      [4]       106       397       130       129       144       103
##      [5]      1344     11025      1466      4172      1896      4181
##      ...       ...       ...       ...       ...       ...       ...
##   [6055]       341      1342       389       739       482       690
##   [6056]       657      5809       753      2112       936      2059
##   [6057]       347      3276       354      1096       481      1042
##   [6058]       139       484       159       144       185       131
##   [6059]        78       170        76       151        87       107
##   -------
##   seqinfo: 7 sequences from an unspecified genome; no seqlengths
```



## Trended Bias


```r
mcols(my_peak_counts) <- mcols(my_peak_counts)[!(grepl("INP", colnames(mcols(my_peak_counts) )))]


par(mfrow=c(1,2), oma=c(3,0,0,0),mar=c(5,4,4,1), mgp = c(2.5,1,0),
    cex=1, cex.axis=1, cex.lab=1.25, cex.main=1.5)

win.ab <- aveLogCPM(DGEList(as.matrix(mcols(my_peak_counts))))
adjc <- log2(as.matrix(mcols(my_peak_counts))+0.5)
logfc <- adjc[,1] - adjc[,4]
smoothScatter(win.ab, logfc, ylim=c(-3, 3), xlim=c(2, 12),
              xlab="Average abundance", ylab="Log-fold change")


offsets <- normOffsets(as.matrix(mcols(my_peak_counts)), type="loess")
norm.adjc <- adjc - offsets/log(2)
norm.fc <- norm.adjc[,1]-norm.adjc[,4]
smoothScatter(win.ab, norm.fc, ylim=c(-3, 3), xlim=c(2, 12),
              xlab="Average abundance", ylab="Log-fold change")
```

<img src="K36_homer_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />



## Estimate Variability


```r
stage <- factor(gsub("[A-B]_.*","", bam.files[!(grepl("INP", bam.files))]))
design <- model.matrix(~0+stage)
colnames(design) <- levels(stage)
design
```

```
##   E03 E11
## 1   1   0
## 2   1   0
## 3   0   1
## 4   0   1
## attr(,"assign")
## [1] 1 1
## attr(,"contrasts")
## attr(,"contrasts")$stage
## [1] "contr.treatment"
```

```r
y <- DGEList(as.matrix(mcols(my_peak_counts)))
y <- scaleOffset(y, offsets)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0005321 0.0007485 0.0010344 0.0017329 0.0025833 0.0042589
```

```r
fit <- glmQLFit(y, design, robust=TRUE)

par(mfrow=c(1,2), oma=c(3,0,0,0),mar=c(5,4,4,1), mgp = c(2.5,1,0),
    cex=1, cex.axis=1, cex.lab=1.25, cex.main=1.5)

plotBCV(y)

plotQLDisp(fit)
```

<img src="K36_homer_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

```r
summary(fit$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.3938 41.7571 41.7654 40.2899 41.7654 41.7654
```



## Differential Test


```r
contrast <- makeContrasts(E11-E03, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

toptags <- topTags(res, n = nrow(res$table), sort.by = "none" )

toptags[1:5,]
```

```
## Coefficient:  -1*E03 1*E11 
##         logFC   logCPM          F       PValue          FDR
## 1 -0.30649641 4.178522  5.4289999 2.448140e-02 3.768618e-02
## 2 -0.18322706 9.105798 26.4985295 5.984431e-06 1.773635e-05
## 3  0.02995419 8.128933  0.3939302 5.349322e-01 5.985028e-01
## 4 -1.13785355 4.052948 64.3916752 3.863933e-10 2.393821e-09
## 5  0.13871712 8.943304 13.9811493 5.356486e-04 1.135582e-03
```


## Region Annotation


```r
out.ranges <- granges(my_peak_counts, use.mcols=FALSE)

mcols(out.ranges) <- cbind(mcols(out.ranges), toptags$table)

out.ranges$direction <- ifelse(out.ranges$logFC > 0, "up", "down")

library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

my_genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
my_genes <- my_genes[seqnames(my_genes) %in% unique(seqnames(out.ranges))]

out.ranges <- out.ranges[seqnames(out.ranges) %in% unique(seqnames(my_genes))]

mid.ranges <- makeGRangesFromDataFrame(data.frame(chr=seqnames(out.ranges),
                                                   start=mid(ranges(out.ranges)),
                                                   end=mid(ranges(out.ranges))))

my_nearest <- nearest(mid.ranges,my_genes)
out.ranges$gene_id <- my_genes[my_nearest]$gene_id

out.ranges$gene_name <- mapIds(org.Dm.eg.db, out.ranges$gene_id, "SYMBOL", keytype="ENSEMBL", multiVals="first")

saveRDS(out.ranges, paste(my_ChIP,"results.rds", sep="."))

out.ranges
```

```
## GRanges object with 6059 ranges and 8 metadata columns:
##          seqnames               ranges strand |       logFC    logCPM
##             <Rle>            <IRanges>  <Rle> |   <numeric> <numeric>
##      [1]    chr2L      [ 7606,   8106]      * | -0.30649641  4.178522
##      [2]    chr2L      [ 8503,  16061]      * | -0.18322706  9.105798
##      [3]    chr2L      [67842,  72637]      * |  0.02995419  8.128933
##      [4]    chr2L      [83112,  83612]      * | -1.13785355  4.052948
##      [5]    chr2L      [94968, 102319]      * |  0.13871712  8.943304
##      ...      ...                  ...    ... .         ...       ...
##   [6055]     chrX [23087925, 23090444]      * |  0.22274553  6.124255
##   [6056]     chrX [23097487, 23102483]      * | -0.03018335  7.974859
##   [6057]     chrX [23105988, 23108518]      * | -0.28327741  7.079496
##   [6058]     chrY [ 3258929,  3259429]      * | -1.09069305  4.310906
##   [6059]     chrY [ 3606553,  3606803]      * |  0.01689416  3.430219
##                    F       PValue          FDR   direction     gene_id
##            <numeric>    <numeric>    <numeric> <character> <character>
##      [1]   5.4289999 2.448140e-02 3.768618e-02        down FBgn0031208
##      [2]  26.4985295 5.984431e-06 1.773635e-05        down FBgn0002121
##      [3]   0.3939302 5.349322e-01 5.985028e-01          up FBgn0067779
##      [4]  64.3916752 3.863933e-10 2.393821e-09        down FBgn0002931
##      [5]  13.9811493 5.356486e-04 1.135582e-03          up FBgn0031216
##      ...         ...          ...          ...         ...         ...
##   [6055]  9.47889342 3.582397e-03 6.514328e-03          up FBgn0039945
##   [6056]  0.46933592 4.969054e-01 5.633786e-01        down FBgn0003559
##   [6057] 27.56446613 4.257563e-06 1.285971e-05        down FBgn0039946
##   [6058] 69.93974749 1.273559e-10 8.535943e-10        down FBgn0085644
##   [6059]  0.01134323 9.156691e-01 9.311916e-01          up FBgn0046323
##            gene_name
##          <character>
##      [1]     CG11023
##      [2]      l(2)gl
##      [3]         dbr
##      [4]         net
##      [5]         Zir
##      ...         ...
##   [6055]     CG17159
##   [6056]       su(f)
##   [6057]        ATbp
##   [6058]        <NA>
##   [6059]         ORY
##   -------
##   seqinfo: 7 sequences from an unspecified genome; no seqlengths
```


