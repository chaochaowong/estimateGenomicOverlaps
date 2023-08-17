setGeneric("estimateGenomicOverlaps", signature = c("query", "subject"),
            function(query, subject, 
                     method = c("poisson_mle", "glm"), rename = NULL,
                     which = NULL, ignore.strand = FALSE, ...)
            standardGeneric("estimateGenomicOverlaps")
)

setMethod("estimateGenomicOverlaps",  c("DataFrame", "character"),
            function(query, subject, 
                     method = c("poisson_mle", "glm"), rename = NULL, 
                     which = NULL, ignore.strand = FALSE, ...)
{
    txdb <- loadFeatures(subject)
    callGeneric(query, txdb, ...)
})


setMethod(estimateGenomicOverlaps, c("DataFrame", "TranscriptDb"),
    function(query, subject,
             method=c("poisson_mle", "glm"), rename = NULL, 
             which = NULL, ignore.strand = FALSE, ...)
{
    method <- match.arg(method)
    dotargs <- list(...)
    if (length(dotargs) != 0L && is.null(names(dotargs)))
        stop("extra arguments must be named")
 
    FUN <- switch(method,
                 poisson_mle=estimateGenomicOverlaps_poisson,
                 glm=estimateGenomicOverlaps_glm)
 
    FUN(query, subject, rename=rename, which=which,
        ignore.strand=ignore.strand, ...)
})

