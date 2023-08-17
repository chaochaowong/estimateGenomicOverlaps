poissonMLE <- function(X, C, kb = rep(1, length(X)), rpm=1, return.read=FALSE)
{
    ## X: observed exons hits
    ## C: gene model - an n-by-m binary matrix where n is teh number of
    ##    transcrpts and m is the numbe of exons. C_ij = 1L if transcript i
    ##    containing exon j, or 0L otherwise.
    ## kb: kilobase
    ## rpm: read per million
    ## C <- t(model)
    if (return.read) { ## ensure that the MLE is not scaled by kb and rpm
        rpm <-1
        kb <- rep(1, length(kb))
    }
    if (is.null(rownames(C))) rownames(C) <- as.character(seq.int(nrow(C)))
        
    #C <- t(model)
    FLAG <- (ncol(C)==1  | nrow(C)==1)
    
    if (FLAG) {
        if (ncol(C) == 1) # one exon (with one or multiple transcripts)
            Theta <- X / (rpm * kb)
        else # one transcript with multiple exons
            Theta <- sum(X) / sum(rpm * kb)
        if (return.read) Theta <- Theta * sum(kb)
        
        names(Theta) <- as.character(rownames(C))
        return(Theta)
    } else {
        ## simple model: transcripts do not share exon
        iniAns <- as.vector(C %*% X )/ (rpm * C %*% kb))
        names(iniAns) <- as.character(rownames(C))
        Theta <- iniAns
        hasValues <- iniAns > 0
        if (sum(hasValues) > 1) { # if more than one tx has value
            if (any(as.logical(colSums(C)-1L))) {
                ## keep the exons that are contained by the transcripts
                keep <- colSums(C[hasValues,, drop=FALSE]) > 0L
                mC <- C[hasValues, keep, drop=FALSE]
                
                op  <- optim(iniAns[hasValues], fn=.mle, gr=.grr,
                         lower=1e-24, method="L-BFGS-B",
                         control=list(fnscale=-1), w=rpm,
                         L=kb[keep], C=mC, X=X[keep])
                if (op$convergence==0)
                    tmpTheta <- op$par
                else
                    tmpTheta <- rep(NA, length(op$par))

                 #names(tmpTheta) <- rownames(mC)
                 Theta[names(tmpTheta)] <- tmpTheta
            }
        }
        if (return.read) Theta <- Theta * rowSums(C)
        return(Theta)

    }

}



 
## mle is the log maximum likelihood used to estimates transcript expression
##   theta = transcript expression estimator
##   w = total number of mapped reads (per million reads)
##   L = length of each exon (per kilobase)
##   C = an n-by-m binary matrix where n is the numbe of transcripts
##       and m is the number of exons. C_ij = TRUE if transcript i containing
##       exon j, or FALSE otherwise.
##   X = observed exon hits
.mle <- function(theta, w=1, L, C, X)
{ 
    thetaC <- theta %*% C

    -w * (thetaC) %*% L + (log(thetaC) + log(L * w)) %*% X -
      sum(lfactorial(X)) 
}

##  the constent term of the mle can be dropped
.mle2 <- function(theta, w=1, L, C, X)
{ 
    thetaC <- theta %*% C

    -w * (thetaC) %*% L + (log(thetaC) + log(L * w)) %*% X
    #- sum(lfactorial(X))

}

## grr is the gradient of the mle function
.grr <- function(theta, w=1, L, C, X)
{ 
    - w * C %*% L + C %*% t(X/theta %*% C)
}



convTxHits2RPKM <- function(grl, gM, hits, rpm)
{   ## hits: a list of hits
    ## gM: gene model -> only for genes with one transcripts
    ## grl: GRagnesList
    x <- unlist(hits[names(grl)])
    kb <- width(grl@unlistData) / 1e3
    ans <- (x %*% gM) / (rpm * (kb %*% gM))
 
}

## .getGeneModel is just used for testing when the input is a GRange instance..
## Build gene model, C, an n-by-m matrix where n is the number of transcripts
## and m is the number of exons. C_ij=1L if transcript i containing exon j,
## or 0L o/w
.getGM4gr <- function(gr)
{
    ## check existance of exon_id and tx_id
    
    ex_id <- values(gr)[["exon_id"]]
    uex <- unique(ex_id)
    
    tx <- values(gr)[["tx_id"]]
    utx <- unique(unlist(tx))
    ex <- rep(ex_id, sapply(tx, length))
    C <- matrix(0L, length(utx), length(uex), dimnames=list(utx, uex))
    C[cbind(match(unlist(tx), utx), match(ex, ex_id))] <- 1L
    C
}
