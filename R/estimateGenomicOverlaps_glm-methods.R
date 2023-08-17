setGeneric("estimateGenomicOverlaps_glm", signature = c("query", "subject"),
            function(query, subject, rename = NULL, 
                     which = NULL, ignore.strand = FALSE,
                     control=glm.control(), ...)
            standardGeneric("estimateGenomicOverlaps_glm")
)

setMethod("estimateGenomicOverlaps_glm",  c("DataFrame", "character"),
            function(query, subject, rename = NULL, 
                     which = NULL, ignore.strand = FALSE,
                     control=glm.control(), ...)
{
    subj <- loadFeatures(subject)
    callGeneric(query, subj, ...)
})

setMethod("estimateGenomicOverlaps_glm", c("DataFrame", "TranscriptDb"),
    function(query, subject, rename = NULL,
             which=NULL, ignore.strand = FALSE, control=glm.control(), ...)
{
   .checkValidQuery(query, checkColName=c("files", "group"),
                    checkFile.exist=TRUE)
   
   preps <- .internal_preprocessExons(query, subject, rename=rename,
       which=which, ignore.strand=ignore.strand, ...)
   Exons <- preps$Exons
   metadata <- preps$metadata
    
  ## need Exons and metadata
    if (!"counts" %in% names(metadata))
        stop("metadata must have a column named 'counts'")
    if (length(metadata[["groups"]]) != length(metadata[["counts"]]))
        stop("metadata group length must equal the number of count data",
             "columns")
    geneLabel <- metadata[["genes"]]
    txLabel <- metadata[["tx"]]
    countLabel <- metadata[["counts"]]

    ## filter low counts, add length covariate
    ## FIXME : should counts be normalized?
    counts <- as.data.frame(values(Exons)[,names(values(Exons)) %in%
        eval(countLabel)])
    if (any(is.na(counts))) stop("Not currently supporting NAs in counts")
    keep <- rowSums(as.data.frame(counts)) > 5
    Exons <- Exons[keep, ]
    counts <- ceiling(counts[keep, , drop=FALSE])
    exonLength <- width(Exons)/1e3
    values(Exons) <- append(values(Exons),
        DataFrame(exonLength=exonLength))

    ## gene models
    ## assumptions : no tx is in > 1 gene group 
    ##               no exon is in > 1 gene group
    geneModels <- .makeGeneModels(Exons, metadata)
    geneGroup <- rep.int(unlist(values(Exons)[, eval(geneLabel)]),
        elementLengths(values(Exons)[, eval(txLabel)]))
    txNames <- unlist(values(Exons)[ , eval(txLabel)])
    mapTxGeneGroup <- DataFrame(geneGroup, txNames)
    mapTxGeneGroup <- mapTxGeneGroup[!duplicated(mapTxGeneGroup[["txNames"]]), ]

    ## glm
    ## FIXME : need to pass both counts and Exons? 
    ## FIXME : loop through multiple file but not combine (ie, no group?)
    ## FIXME : if too slow, investigate analytical soln for 1 exon or 1 tx models 
    lst <- lapply(unique(mapTxGeneGroup[["geneGroup"]]), .glmTx,
                  geneModels, metadata, mapTxGeneGroup, counts, Exons)
    names(lst) <- unique(mapTxGeneGroup[["geneGroup"]])
    do.call(rbind, lst)
})


.estTheta <- function(data, formula)
{
    nG <- length(levels(factor(data$group)))
    nTx <- ncol(data$tx)

    ## estimate mu
    if (nG <=1 & nTx <= 1)
        mu <- mean(data$counts)
    else
        mu <- fitted(glm(formula, data, family=poisson))
      
    theta <- try(theta.ml(data$counts, mu, limit=20),
                 silent=TRUE)
    
}


.glmTx <- function(i, geneModels, metadata, mapTxGeneGroup, counts, Exons,
                   control=glm.control(), ...)
{
    geneLabel <- metadata[["genes"]]
    nsamples <- length(metadata[["groups"]])
    geneIdx <- unlist(values(Exons)[ , eval(geneLabel)] == i)
    if (is.null(dim(counts))) {
        y <- counts[geneIdx]
    } else {
        y <- unname(unlist(counts[geneIdx, ]))
    }
    length <- unlist(values(Exons)[["exonLength"]])[geneIdx]
    model <- geneModels[geneIdx, ,drop=FALSE]
    minimalx <- model[, colSums(model) != 0, drop=FALSE]
    txn <- mapTxGeneGroup[["txNames"]][
        mapTxGeneGroup[["geneGroup"]] == i]
    dimnames(minimalx) <- list(NULL, txn)
    x <- minimalx[rep(seq_len(nrow(minimalx)), nsamples), ,drop=FALSE]
    group <- as.factor(metadata[["groups"]])
    data <- data.frame(tx=I(x), counts=y, group=group,
        length=rep.int(length, nsamples))
   
    ## define formula, result row and column names 
    if (length(levels(as.factor(metadata[["groups"]]))) > 1) {
        formula <- counts ~ tx + as.factor(group) -1
        rn <- c(paste("tx", mapTxGeneGroup[["txNames"]][
            mapTxGeneGroup[["geneGroup"]] == i],
            sep=""), paste("group", levels(group), sep=""))
    } else {
        formula <- counts ~ tx -1 
        rn <- c(paste("tx", mapTxGeneGroup[["txNames"]][
            mapTxGeneGroup[["geneGroup"]] == i],
            sep=""))
    }
    cn <- c("Estimate", "Std.Error", "z value", "Pr(>|z|)")
  
    ## single exon case
    if (nrow(x) == 1){
        res <- matrix(NA, length(rn), 4L, dimnames = list(rn, cn))
    } else {
        ## estimate theta 
        theta <- .estTheta(data, formula)
        ## glm.fit might not converge when using neg. binomial
        if (inherits(theta, "try-error"))
            fit <- try(glm(formula, family=poisson,
                       data=data, control=control), silent=TRUE)
        else
            fit <- try(glm(formula, family=negative.binomial(theta),
                      data=data, control=control), silent=TRUE)
        if (!inherits(fit, "try-error")) {
                if (!fit$converged)
                    fit <- try(glm(formula, family=poisson,
                        data=data, control=control), silent=TRUE)
        }

        ## results
        if (inherits(fit, "try-error")) {
            res <- matrix(NA, length(rn), 4L, dimnames = list(rn, cn))
        } else {
            res <- matrix(NA, length(rn), 4L, dimnames = list(rn, cn))
            nonNAresults <-
                unlist(dimnames(summary(fit)[["coefficients"]])[1])
            if (ncol(minimalx) == 1)
                nonNAresults[1] <- paste("tx", colnames(minimalx), sep="")
            res[rn %in% nonNAresults] <- summary(fit)[["coefficients"]]
        }
    }
    round(res, 4)
}



