## 
## unit tests for estimateGenomicOverlaps 
## 

library(GenomicFeatures)
library(MASS)
source("estimateGenomicOverlaps_glm-methods.R")
source("estimateGenomicOverlaps-methods.R")
source("glmEstimateTx.R")
source("cleanSeq.R")


## A : 2 exon, 1 tx, same length
## B : 4 exon, 3 tx, same length, 2 exons hit same tx
## C : 2 exon, 2 tx, different length *not sure
## D1 : 2 exon, 2 tx, each hit one tx, same length, same counts
## D2 : 2 exon, 2 tx, each hit one tx, same length, different counts
## D3 : 2 exon, 2 tx, each hit one tx, different length, same counts
## D4 : 2 exon, 2 tx, each hit one tx, different length, different counts
## E1 : 2 exon, 2 tx, one exon hits both tx, same length, same counts
## E2 : 2 exon, 2 tx, one exon hits both tx, same length, different counts
## E3 : 2 exon, 2 tx, one exon hits both tx, different length, same counts
## E4 : 2 exon, 2 tx, one exon hits both tx, different length, different counts
rng1 <- function(s, w)
GRanges(seq="chr1", IRanges(s, width=w), strand="+")

gr <- c(
    ## A, B, C
    rng1(c(1000, 2000), c(10, 10)),
    rng1(c(3000, 3050, 4000, 4050), rep.int(50, 4)),
    rng1(c(5000, 6000), c(10, 200)),
    ## D's
    rng1(c(7000, 7100), c(50, 50)),
    rng1(c(7200, 7300), c(50, 50)),
    rng1(c(7400, 7500), c(50, 100)),
    rng1(c(7700, 7800), c(50, 100)),
    ## E's
    rng1(c(8000, 8100), c(50, 50)),
    rng1(c(8200, 8300), c(50, 50)),
    rng1(c(8400, 8600), c(150, 50)),
    rng1(c(8700, 8900), c(150, 50)))

counts <- c(
    ## A, B, C
    100, 100, 100, 100, 100, 100, 100, 200, 
    ## D's 
    100, 100, 100, 200, 100, 100, 100, 200,
    ## E's 
    100, 100, 200, 100, 100, 100, 200, 100)

metaData <- DataFrame(
    exons = seq_len(length(gr)),
    tx_id = IntegerList(list(
            ## A, B, C
            83, 84, 85, 86, 87, 87, 88, 89, 
            ## D's
            90, 91, 92, 93, 94, 95, 96, 97,
            ## E's 
            c(98, 99), 99, c(100, 101), 101, 
            c(102, 103), 103, c(104, 105), 105)),
    genes = c(rep.int("A", 2), rep.int("B", 4), rep.int("C", 2),
              rep(c("D1", "D2", "D3", "D4"), each=2), 
              rep(c("E1", "E2", "E3", "E4"), each=2)), 
    samp1 = counts, 
    samp2 = counts + runif(length(counts), 1, 10),
    samp3 = counts + runif(length(counts), -10, 0),
    samp4 = counts + runif(length(counts), -15, -5))
values(gr) <- metaData

DF <- list(
    #groups = c(0 , 0, 1, 1),
    groups = 1,
    # counts = c("samp1", "samp2", "samp3", "samp4"),
    counts = c("samp1"),
    genes = "genes",
    tx = "tx_id")


test_glmEstimateTx <- function()
{
    DF <- list(groups = 1, counts = "samp1", genes = "genes", tx = "tx_id")
    res <- glmEstimateTx(gr, DF)
 
   # ## Case 1 : 2 exons, 2 tx, exons hit different tx 
   #     ## Gene A : same length (10), same counts (100)
   #     A <- res[["A"]][,1] 
   #     checkIdentical(c(10, 10), A)
   #     ## Gene D1 : same length (50), same counts (100)
   #     D1 <- res[["D1"]][,1] 
   #     checkIdentical(c(2, 2), D1)
   #     ## Gene D2 : same length (50), different counts (100, 200)
   #     D2 <- res[["D2"]][,1] 
   #     checkIdentical(c(2, 4), D2)
   #     ## Gene D3 : different length (50, 100), same counts (100)
   #     D3 <- res[["D3"]][,1] 
   #     checkIdentical(c(2, 1), D3)
   #     ## Gene D4 : different length (50, 100), different counts (100, 200)
   #     D4 <- res[["D4"]][,1] 
   #     checkIdentical(c(2, 1), D4)
 
   # ## Case 2 : 2 exons, 2 tx, 1 exon hits both tx 
   #     ## Gene E1 : exons same length(50), same counts (100)
   #     E1 <- res[["E1"]][,1] 
   #     checkIdentical(c(1, 2), E1)
   #     ## Gene E2 : exons same length (50), different counts (200, 100)
   #     E2 <- res[["E2"]][,1] 
   #     checkIdentical(c(2, 2), E2)
   #     ## Gene E3 : exons different length (150, 50), same counts (100)
   #     E3 <- res[["E3"]][,1] 
   #     checkIdentical(c(0.5, 2), E3)
   #     ## Gene E4 : exons diff length (150, 50), diff counts (200, 100) 
   #     E4 <- res[["E4"]][,1] 
   #     checkIdentical(c(1, 2), E4)
}


