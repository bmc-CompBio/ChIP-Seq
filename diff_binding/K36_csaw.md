# csaw



## BAM Files



```r
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
param <- readParam(minq=0, BPPARAM = MulticoreParam(4))

win.data <- windowCounts(bam.files.path, param=param, width=250, ext=150)
bin.data <- windowCounts(bam.files.path, bin=TRUE, param=param, width=10000)

saveRDS(win.data, file = paste(my_ChIP,"win.data.rds", sep="."))
saveRDS(bin.data, file = paste(my_ChIP,"bin.data.rds", sep="."))

win.data
```

```
## class: RangedSummarizedExperiment 
## dim: 2515498 8 
## metadata(6): spacing width ... param final.ext
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames: NULL
## colData names(4): bam.files totals ext rlen
```

```r
metadata(win.data)
```

```
## $spacing
## [1] 50
## 
## $width
## [1] 250
## 
## $shift
## [1] 0
## 
## $bin
## [1] FALSE
## 
## $param
##     Extracting reads in single-end mode
##     Duplicate removal is turned off 
##     Minimum allowed mapping score is 0 
##     Reads are extracted from both strands
##     No restrictions are placed on read extraction
##     No regions are specified to discard reads
##     Using MulticoreParam with 4 workers
## $final.ext
## [1] NA
```

```r
colData(win.data)
```

```
## DataFrame with 8 rows and 4 columns
##                         bam.files    totals       ext      rlen
##                       <character> <integer> <integer> <integer>
## 1 ../BAM//E03A_INP_ATCACGAT_1.bam  20069531       150        50
## 2 ../BAM//E03A_K36_TTAGGCAT_1.bam  19728462       150        50
## 3 ../BAM//E03B_INP_CGATGTAT_1.bam  19478697       150        50
## 4 ../BAM//E03B_K36_ACAGTGAT_1.bam  36394416       150        50
## 5 ../BAM//E11A_INP_TTAGGCAT_1.bam  20731254       150        50
## 6 ../BAM//E11A_K36_GCCAATAT_1.bam  19856620       150        50
## 7 ../BAM//E11B_INP_ACAGTGAT_1.bam  25354820       150        50
## 8 ../BAM//E11B_K36_GATCAGAT_1.bam  18435605       150        50
```

```r
rowRanges(win.data)
```

```
## GRanges object with 2515498 ranges and 0 metadata columns:
##             seqnames             ranges strand
##                <Rle>          <IRanges>  <Rle>
##         [1]    chr2L       [4801, 5050]      *
##         [2]    chr2L       [4851, 5100]      *
##         [3]    chr2L       [4901, 5150]      *
##         [4]    chr2L       [4951, 5200]      *
##         [5]    chr2L       [5001, 5250]      *
##         ...      ...                ...    ...
##   [2515494]     chrY [3667151, 3667352]      *
##   [2515495]     chrY [3667201, 3667352]      *
##   [2515496]     chrY [3667251, 3667352]      *
##   [2515497]     chrY [3667301, 3667352]      *
##   [2515498]     chrY [3667351, 3667352]      *
##   -------
##   seqinfo: 8 sequences from an unspecified genome
```

```r
head(assay(win.data))
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
## [1,]    5    6    5    6    5    9    7    3
## [2,]    9   13   11   19    8   13   10    6
## [3,]   17   19   20   35   19   21   22   14
## [4,]   29   30   29   52   27   36   28   23
## [5,]   39   41   37   74   35   47   35   36
## [6,]   51   48   47   87   50   56   46   42
```


```r
ChIP.win.data <- win.data[,!(grepl("INP", bam.files))]
INP.win.data <- win.data[,grep("INP", bam.files)]

ChIP.bin.data <- bin.data[,!(grepl("INP", bam.files))]
INP.bin.data <- bin.data[,grep("INP", bam.files)]
```

## Filtering



```r
filter.stat <- filterWindows(data =       ChIP.win.data, 
                             background = INP.win.data, 
                             type = "control", 
                             prior.count = 5, 
                             norm.fac = list(ChIP.bin.data, INP.bin.data))
min.fc <- 2
keep <- filter.stat$filter > log2(min.fc)

ChIP.filt.data <- ChIP.win.data[keep,]

saveRDS(ChIP.filt.data, paste(my_ChIP,"filt.data.rds", sep="."))

summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical 2125893  389605
```

```r
library(RColorBrewer)
my_colors <- brewer.pal(9,"Set1")

par(mfrow=c(1,2),  oma=c(3,0,0,0),mar=c(5,4,4,1), mgp = c(2.5,1,0),
    cex=1, cex.axis=1, cex.lab=1.25, cex.main=1.5)


plot(density(filter.stat$back.abundances), xlim=c(-1,5), 
     type="l", col=my_colors[9], lwd=3,
     xlab="log2 CPM", main="Pre-Filtering")
lines(density(filter.stat$abundances), col=my_colors[5], lwd=3)

plot(density(filter.stat$back.abundances[keep]), xlim=c(-1,5), 
     type="l", col=my_colors[9], lwd=3,
     xlab="log2 CPM", main="Post-Filtering")
lines(density(filter.stat$abundances[keep]), col=my_colors[5], lwd=3)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)

plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-0.33,-0.75, legend =  c("ChIP","Input"),
       xpd = TRUE, horiz = TRUE, inset = c(2, 0), bty = "n", pch = 19, col = my_colors[c(5,9)], cex = 1.25)
```

<img src="K36_csaw_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

```r
par(mfrow=c(1,1), oma=c(3,0,0,0),mar=c(5,4,4,1), mgp = c(2.5,1,0),
    cex=1, cex.axis=1, cex.lab=1.25, cex.main=1.5)


my_region <- 1:2000
keep_region <- keep[my_region]

ylim <- ceiling(max(filter.stat$abundances[my_region]))*1.2

plot(start(rowRanges(win.data))[my_region], 
     filter.stat$abundances[my_region],
     type="l", col=my_colors[5], lwd=2, ylim=c(1,ylim),
     xlab="Coordinates", ylab="log2 CPM", main = "Example Region")

xx <- c(start(rowRanges(win.data))[my_region], rev(start(rowRanges(win.data))[my_region]))
yy <- c(rep(0,length(my_region)), rev(filter.stat$abundances[my_region]))
polygon(xx, yy, col=my_colors[5], lty = 0)
lines(filter.stat$back.abundances[my_region],type="l", col=my_colors[9], lwd=2)

xx <- c(start(rowRanges(win.data))[my_region], rev(start(rowRanges(win.data))[my_region]))
yy <- c(rep(0,length(my_region)), rev(filter.stat$back.abundances[my_region]))
polygon(xx, yy, col=my_colors[9], lty = 0)

xx <- c(start(rowRanges(win.data))[my_region], rev(start(rowRanges(win.data))[my_region]))
yy <- c(rep(0,length(my_region)), rev(keep_region*(ylim/20)))+(ylim*0.95)
polygon(xx, yy, col=my_colors[5], lty = 0)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)

plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-0.33,-0.75, legend =  c("ChIP","Input"),
       xpd = TRUE, horiz = TRUE, inset = c(2, 0), bty = "n", pch = 19, col = my_colors[c(5,9)], cex = 1.25)
```

<img src="K36_csaw_files/figure-html/unnamed-chunk-4-2.png" style="display: block; margin: auto;" />

## Trended Bias


```r
par(mfrow=c(1,2), oma=c(3,0,0,0),mar=c(5,4,4,1), mgp = c(2.5,1,0),
    cex=1, cex.axis=1, cex.lab=1.25, cex.main=1.5)

win.ab <- aveLogCPM(asDGEList(ChIP.filt.data))
adjc <- log2(assay(ChIP.filt.data)+0.5)
logfc <- adjc[,1] - adjc[,4]
smoothScatter(win.ab, logfc, ylim=c(-4, 4), xlim=c(1, 5),
    xlab="Average abundance", ylab="Log-fold change")


offsets <- normOffsets(ChIP.filt.data, type="loess")
norm.adjc <- adjc - offsets/log(2)
norm.fc <- norm.adjc[,1]-norm.adjc[,4]
smoothScatter(win.ab, norm.fc, ylim=c(-4, 4), xlim=c(1, 5),
    xlab="Average abundance", ylab="Log-fold change")
```

<img src="K36_csaw_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

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
y <- asDGEList(ChIP.filt.data)
y <- scaleOffset(y, offsets)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0008006 0.0010608 0.0013904 0.0018307 0.0024503 0.0050647
```

```r
fit <- glmQLFit(y, design, robust=TRUE)

par(mfrow=c(1,2), oma=c(3,0,0,0),mar=c(5,4,4,1), mgp = c(2.5,1,0),
    cex=1, cex.axis=1, cex.lab=1.25, cex.main=1.5)

plotBCV(y)

plotQLDisp(fit)
```

<img src="K36_csaw_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

```r
summary(fit$df.prior)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.3218 111.8252 111.8252 111.5266 111.8252 111.8252
```



## Differential Test



```r
contrast <- makeContrasts(E11-E03, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
head(res$table)
```

```
##          logFC   logCPM           F    PValue
## 1 -0.043737440 2.134307 0.063850851 0.8009657
## 2  0.009903156 2.119649 0.003290543 0.9543564
## 3 -0.026324947 2.241569 0.024879801 0.8749464
## 4  0.031531013 2.280505 0.037485498 0.8468250
## 5  0.072039269 2.347003 0.206854684 0.6501102
## 6  0.126838549 2.399358 0.667434045 0.4156543
```

```r
merged <- mergeWindows(rowRanges(ChIP.filt.data), tol=100, max.width=5000)
tabcom <- combineTests(merged$id, res$table)

head(tabcom)
```

```
##   nWindows logFC.up logFC.down       PValue          FDR direction
## 1        2        0          0 9.543564e-01 9.688800e-01     mixed
## 2       88        0          1 2.213927e-04 4.728013e-04      down
## 3       87        0          8 3.337728e-06 1.049756e-05      down
## 4       91        0          9 7.381236e-06 2.179870e-05      down
## 5       15        0         12 7.120959e-09 3.812280e-08      down
## 6       77        0          0 5.103529e-02 6.782729e-02        up
```

```r
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)
```

```
##    Mode   FALSE    TRUE 
## logical    2311    6228
```

```r
table(tabcom$direction[is.sig])
```

```
## 
##  down mixed    up 
##  2715    94  3419
```

```r
tabbest <- getBestTest(merged$id, res$table)
head(tabbest)
```

```
##   best       logFC   logCPM           F       PValue          FDR
## 1    1 -0.04373744 2.134307  0.06385085 1.000000e+00 1.000000e+00
## 2   84 -0.50000849 3.762860 23.62088696 3.330082e-04 7.019810e-04
## 3  166 -0.75140749 3.084744 34.79488062 3.337728e-06 1.065055e-05
## 4  263 -0.78573331 2.765574 30.81011524 1.713317e-05 4.792013e-05
## 5  269 -1.17775991 2.282769 46.32112601 7.423893e-09 4.017277e-08
## 6  292  0.30633815 4.071351 12.25749842 5.103529e-02 6.987179e-02
```




## Region Annotation


```r
out.ranges <- merged$region

elementMetadata(out.ranges) <- data.frame(tabcom,
    best.pos=mid(ranges(rowRanges(ChIP.filt.data[tabbest$best]))),
    best.logFC=tabbest$logFC,
    ave.logFC=sapply(1:nrow(tabcom), function(i){mean(res$table$logFC[merged$id == i]) })
    )


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
## GRanges object with 8539 ranges and 11 metadata columns:
##          seqnames             ranges strand |  nWindows  logFC.up
##             <Rle>          <IRanges>  <Rle> | <integer> <integer>
##      [1]    chr2L     [ 6301,  6600]      * |         2         0
##      [2]    chr2L     [ 7351, 11950]      * |        88         0
##      [3]    chr2L     [11751, 16300]      * |        87         0
##      [4]    chr2L     [67801, 72550]      * |        91         0
##      [5]    chr2L     [83051, 84000]      * |        15         0
##      ...      ...                ...    ... .       ...       ...
##   [8535]     chrY [3576051, 3576950]      * |        11         0
##   [8536]     chrY [3577701, 3579300]      * |        22         0
##   [8537]     chrY [3579801, 3580350]      * |         7         0
##   [8538]     chrY [3606351, 3607050]      * |        10         3
##   [8539]     chrY [3625201, 3625450]      * |         1         0
##          logFC.down       PValue          FDR   direction  best.pos
##           <integer>    <numeric>    <numeric> <character> <integer>
##      [1]          0 9.543564e-01 9.688800e-01       mixed      6425
##      [2]          1 2.213927e-04 4.728013e-04        down     11525
##      [3]          8 3.337728e-06 1.049756e-05        down     15625
##      [4]          9 7.381236e-06 2.179870e-05        down     72175
##      [5]         12 7.120959e-09 3.812280e-08        down     83175
##      ...        ...          ...          ...         ...       ...
##   [8535]         11 0.0395606003 0.0538768685        down   3576825
##   [8536]         17 0.0059037740 0.0094725091        down   3578725
##   [8537]          7 0.0003341019 0.0006887726        down   3579975
##   [8538]          0 0.0028749554 0.0049236350          up   3606725
##   [8539]          1 0.0555435987 0.0733020679        down   3625325
##           best.logFC   ave.logFC     gene_id   gene_name
##            <numeric>   <numeric> <character> <character>
##      [1] -0.04373744 -0.01691714 FBgn0031208     CG11023
##      [2] -0.50000849 -0.14047097 FBgn0031208     CG11023
##      [3] -0.75140749 -0.25243921 FBgn0002121      l(2)gl
##      [4] -0.78573331 -0.03072580 FBgn0067779         dbr
##      [5] -1.17775991 -0.69850134 FBgn0002931         net
##      ...         ...         ...         ...         ...
##   [8535]  -1.0534218  -0.8798012 FBgn0046323         ORY
##   [8536]  -1.2167192  -0.8425203 FBgn0046323         ORY
##   [8537]  -1.8172350  -1.5552453 FBgn0046323         ORY
##   [8538]   0.5124549   0.3765305 FBgn0046323         ORY
##   [8539]  -0.8085177  -0.8085177 FBgn0046323         ORY
##   -------
##   seqinfo: 8 sequences from an unspecified genome
```



