
load("/Users/cwon2/Project/estimateGenomicOverlaps/R/grl.rda")
# which one has more then one transcripts
# 39
# 7
## case grl[["7"]]
#source("poissonMLE.R")
gr <- grl[["7"]]
values(gr)$exon_id <- seq.int(length(gr))
values(gr)$tx_id <- values(gr)$tx
C <- estimateGenomicOverlaps:::.getGM4gr(gr)
X <- values(gr)$V1
kb <- width(gr) / 1e3
rpm <- 1
poissonMLE(X, C, kb, rpm, return.read=TRUE)

## case grl[["39"]]
gr <- grl[["39"]]
values(gr)$exon_id <- seq.int(length(gr))
values(gr)$tx_id <- values(gr)$tx
C <- estimateGenomicOverlaps:::.getGM4gr(gr)
X <- values(gr)$V1
kb <- width(gr) / 1e3
rpm <- 1
poissonMLE(X, C, kb, rpm, return.read=TRUE)

## trivial case
C <- matrix(c(1,1,0,1), nrow=2)
rownames(C) <- c("A", "B")
X <- c(20, 10)
poissonMLE(X=X, C=C, return.read=TRUE)

## yeast case
#source("poissonMLE.R")
txdbFile <- "/Users/cwon2/Project/cwon2/TransExpression/extdata/sacCer2_sgdGene.sqlite"
txdb <- loadFeatures(txdbFile)
bamFile <- "/Users/cwon2/Project/cwon2/TransExpression/extdata/SRR002051.chrI-V.bam"
query <- DataFrame(files=c(bamFile, bamFile))
which <- GRanges(c("chrI", "chrII", "chrIII"), IRanges(1, 1e+07))
rename=NULL
res <- estimateGenomicOverlaps_poisson(query[1], subject=txdb, rename=NULL, which=which)

### human case
bamFile <- "/Users/cwon2/Project/accepted_hits.bam"
txdb <-    "/Users/cwon2/Project/hg18txdbrefGene.sqlite"
which <- GRanges(c("chr3", "chr4", "chr5"), IRanges(1, 1e+7))
rename=NULL

query <- DataFrame(files=bamFile)

res <- estimateGenomicOverlaps_poisson(query, subject=txdb,
                                       rename=NULL, which=which)
