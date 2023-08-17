setGeneric("estimateGenomicOverlaps_poisson", signature = c("query", "subject"),
    function(query, subject, rename=NULL, keep=NULL,
             which=NULL, ignore.strand = FALSE, ...)
    standardGeneric("estimateGenomicOverlaps_poisson")
)

setMethod("estimateGenomicOverlaps_poisson",  c("DataFrame", "character"),
    function(query, subject,
             rename=NULL, keep=NULL, which=NULL,
             ignore.strand = FALSE, ...)
          
{
    if (!file.exists(subject)) 
        stop("The TranscriptDb file does not exist")
    else {
      subj <- loadFeatures(subject)
      callGeneric(query, subj, ...)
    }
})

setMethod("estimateGenomicOverlaps_poisson",
          c("DataFrame", "TranscriptDb"),
    function(query, subject,
             rename=NULL, keep=NULL, which=NULL,
             ignore.strand = FALSE, ...)
{
    ## FIXME : Eventually want choice of count resolution and type
    .checkValidQuery(query, checkColName="files", checkFile.exist=TRUE)

    preps <- .internal_preprocessExons(query, subject, rename=rename,
           keep=keep, ignore.strand=ignore.strand, which=which, ...)
    Exons <- preps$Exons
    metadata <- preps$metadata
    
    if (!"counts" %in% names(metadata))
        stop("metadata must have a column named 'counts'")

    geneLabel <- metadata[["genes"]]
    txLabel <- metadata[["tx"]]
    countLabel <- metadata[["counts"]]
    counts <- as.data.frame(values(Exons)[,names(values(Exons)) %in%
        countLabel])
    colnames(counts) <- countLabel
    ## gene models
    ## assumptions : no tx is in > 1 gene group 
    ##               no exon is in > 1 gene group
    gOneTx <- .getMapsTxGene(Exons, geneLabel, txLabel, return.gOneTx=TRUE)

    ## subsetting into two GRanges
    geneModels <- .glmGeneModels(Exons, metadata) 

    .getMapsTxGene <- function(Exons, geneLabel, exLabel, return.gOneTx=FALSE) {
    geneGroup <- rep.int(unlist(values(Exons)[, eval(geneLabel)]),
        elementLengths(values(Exons)[, eval(txLabel)]))
    txNames <- unlist(values(Exons)[ , eval(txLabel)])
    mapTxGeneGroup <- DataFrame(geneGroup, txNames)
    mapTxGeneGroup <- mapTxGeneGroup[!duplicated(mapTxGeneGroup[["txNames"]]), ]

     ## identify single tx models 
        singleTxRle <- Rle(mapTxGeneGroup[["geneGroup"]])
        singleGroups <- runValue(singleTxRle)[runLength(singleTxRle) == 1]
        singleIdx <- mapTxGeneGroup[["geneGroup"]] %in% singleGroups
    mapTxGeneGroup[["singleTx"]] <- FALSE
    mapTxGeneGroup[["singleTx"]][singleIdx] <- TRUE
    
    if (return.gOneTx)
       return(singleGroups)
    else
       return(mapTxGeneGroup=mapTxGeneGroup) 
    }  
    
    ## kbLength 

    if(any(is.na(counts))) stop("Not currently supporting NAs in counts")
    kbLength <- width(Exons)/1e3
    values(Exons) <- append(values(Exons), DataFrame(kbLength=kbLength))

    ## poisson MLE 
    ## for each file 
    lapply(countLabel, 
        function(c, mapTxGenegroup, geneModels, metadata,
                 counts, Exons) {
            subcounts <- counts[, names(counts) == eval(c), drop=FALSE]
            rpm <- (metadata[["totalReads"]]/1e6)[metadata[["counts"]] == eval(c)]
            ## FIXME : - insert analytical solns for cases of 1 exon and 1 tx
            ##         - subset counts and/or Exons on 'singleTx' 

            ## for each gene model
            lst <- lapply(unique(mapTxGeneGroup[["geneGroup"]]),
                function(i, geneModels, metadata, mapTxGeneGroup, 
                         subcounts, rpm, Exons) {
                    geneIdx <- unlist(values(Exons)[,eval(geneLabel)] == i)
                    model <- geneModels[geneIdx, ,drop=FALSE]
                    x <- model[, colSums(model) != 0, drop=FALSE]
                    txn <- mapTxGeneGroup[["txNames"]][
                        mapTxGeneGroup[["geneGroup"]] == i]
                    dimnames(x) <- list(NULL, txn)
                    hits <- subcounts[geneIdx, names(counts) == eval(c)]
                    kbL <- values(Exons)[["kbLength"]][geneIdx]
                    expression <- poissonMLE(hits, x, kbL, rpm) 
                  }, geneModels, metadata, mapTxGeneGroup, subcounts,
                  rpm, Exons)

        round(unlist(lst), 4)
        }, mapTxGeneGroup, geneModels, metadata, counts, Exons)

} )


.getMapsTxGene <- function(Exons, geneLabel, exLabel, return.gOneTx=FALSE) {
    geneGroup <- rep.int(unlist(values(Exons)[, eval(geneLabel)]),
        elementLengths(values(Exons)[, eval(txLabel)]))
    txNames <- unlist(values(Exons)[ , eval(txLabel)])
    mapTxGeneGroup <- DataFrame(geneGroup, txNames)
    mapTxGeneGroup <- mapTxGeneGroup[!duplicated(mapTxGeneGroup[["txNames"]]), ]

     ## identify single tx models 
    singleTxRle <- Rle(mapTxGeneGroup[["geneGroup"]])
    singleGroups <- runValue(singleTxRle)[runLength(singleTxRle) == 1]
    singleIdx <- mapTxGeneGroup[["geneGroup"]] %in% singleGroups
    mapTxGeneGroup[["singleTx"]] <- FALSE
    mapTxGeneGroup[["singleTx"]][singleIdx] <- TRUE
    
    if (return.gOneTx)
       return(singleGroups)
    else
       return(mapTxGeneGroup=mapTxGeneGroup) 
}
    
